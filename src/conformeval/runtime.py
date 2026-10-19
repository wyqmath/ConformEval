"""Run provenance and safe paths shared by inference and evaluation."""
from __future__ import annotations
import hashlib
import json
import platform
try:
    import resource
except ImportError:
    resource = None
import time
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def json_write(path: Path, obj) -> None:
    Path(path).write_bytes((json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n').encode('utf-8'))


def safe_path(root: Path, name: str) -> Path:
    rel = Path(name)
    path = (root / rel).resolve()
    if rel.is_absolute() or '..' in rel.parts or not path.is_relative_to(root.resolve()):
        raise ValueError(f'unsafe archive path: {name}')
    return path


def versions() -> dict:
    result = {'python': platform.python_version()}
    for name in ['numpy', 'scipy', 'scikit-learn', 'torch', 'transformers', 'sentencepiece', 'huggingface-hub', 'tzdata']:
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            pass
    return result


class RunRecord:
    """Save actual inputs, configuration and outcomes, including failed runs."""
    def __init__(self, output: Path, configuration: dict, code_files=None):
        self.output = output
        if output.exists():
            raise ValueError(f'refusing to overwrite existing output directory: {output}')
        output.mkdir(parents=True)
        self.clock = time.monotonic()
        self.code_files = [Path(__file__)] + list(code_files or [])
        self.data = {
            'project': 'ConformEval', 'status': 'RUNNING', 'started_at': self.stamp(),
            'configuration': configuration, 'dependencies': versions(),
            'hardware': {'platform': platform.platform(), 'machine': platform.machine()},
            'input_sha256': {}, 'warnings': [],
            'code_sha256': {p.name: digest(p) for p in self.code_files},
        }
        self.save()

    @staticmethod
    def stamp():
        return datetime.now().astimezone().isoformat(timespec='seconds')

    def track(self, path: Path):
        self.data['input_sha256'][str(path.resolve())] = digest(path)

    def save(self):
        json_write(self.output / 'run.json', self.data)

    def finish(self, error=None):
        if error is None:
            for name, expected in self.data['input_sha256'].items():
                if digest(Path(name)) != expected:
                    error = ValueError(f'input changed during run: {name}')
                    break
            current = {p.name: digest(p) for p in self.code_files}
            if current != self.data['code_sha256']:
                error = ValueError('evaluation or extraction code changed during run')
        self.data.update(status='FAILED' if error else 'PASS', ended_at=self.stamp(),
                         seconds=time.monotonic() - self.clock,
                         peak_rss_bytes=(int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss *
                                            (1 if platform.system() == 'Darwin' else 1024)) if resource else None))
        if error:
            self.data.update(error_type=type(error).__name__, error=str(error))
        self.data['output_sha256'] = {p.relative_to(self.output).as_posix(): digest(p)
                                     for p in sorted(self.output.rglob('*'))
                                     if p.is_file() and p.name != 'run.json'}
        self.save()
        return error
