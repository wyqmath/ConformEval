"""UTF-8/BOM, newline and real scoring checks without encoder dependencies."""
from __future__ import annotations
import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from conformeval.adapters import load_config
from conformeval.data_validation import read_tsv as read_data_tsv
from conformeval.evaluate_predictions import read_table, write_gzip_tsv
from conformeval.readers import parse_fasta
from conformeval.residue_scoring_core import read_tsv, write_tsv
from conformeval.runtime import json_write

ROOT = Path(__file__).resolve().parents[1]


def read_data_table(path):
    return read_data_tsv(path, compressed=path.suffix=='.gz')


class TextIOTests(unittest.TestCase):
    def test_tsv_bom_crlf_and_unicode_plain_or_gzip(self):
        rows=[{'id':'protein_α','note':'assay\r\n25 °C'}]
        for ending in ['\n','\r\n']:
            for bom in [b'',b'\xef\xbb\xbf']:
                stream=io.StringIO(newline='')
                writer=csv.DictWriter(stream,['id','note'],delimiter='\t',lineterminator=ending)
                writer.writeheader();writer.writerows(rows)
                raw=bom+stream.getvalue().encode('utf-8')
                with tempfile.TemporaryDirectory() as temp:
                    for suffix in ['.tsv','.tsv.gz']:
                        path=Path(temp)/('labels'+suffix)
                        path.write_bytes(gzip.compress(raw) if suffix.endswith('.gz') else raw)
                        for reader in [read_table,read_data_table,read_tsv]:
                            with self.subTest(ending=repr(ending),bom=bool(bom),suffix=suffix,reader=reader.__module__):
                                self.assertEqual(reader(path),rows)

    def test_json_fasta_config_and_output_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)/'Unicode α path';root.mkdir()
            config=json.loads((ROOT/'configs/models/esm2_8m.json').read_text(encoding='utf-8'))
            config['note']='assay temperature: 25 °C'
            path=root/'model_α.json';path.write_bytes(b'\xef\xbb\xbf'+json.dumps(config,ensure_ascii=False).encode('utf-8'))
            self.assertEqual(load_config(path)['note'],'assay temperature: 25 °C')
            fasta=root/'sequences_α.fasta';fasta.write_bytes('\ufeff>protein_α\r\nACDE\r\n'.encode('utf-8'))
            self.assertEqual(parse_fasta(fasta),{'protein_α':'ACDE'})
            target=root/'output_α.json';json_write(target,{'temperature_°C':'25 °C','value':1})
            self.assertEqual(json.loads(target.read_bytes().decode('utf-8')),{'temperature_°C':'25 °C','value':1})
            self.assertTrue(target.read_bytes().endswith(b'\n'));self.assertNotIn(b'\r\n',target.read_bytes())
            rows=[{'id':'protein_α','value':'1'}]
            for suffix in ['.tsv','.tsv.gz']:
                path=root/('table'+suffix);write_tsv(path,['id','value'],rows)
                raw=gzip.decompress(path.read_bytes()) if suffix.endswith('.gz') else path.read_bytes()
                self.assertEqual(raw,'id\tvalue\nprotein_α\t1\n'.encode('utf-8'))
            path=root/'scores.tsv.gz';write_gzip_tsv(path,rows,['id','value'])
            self.assertEqual(gzip.decompress(path.read_bytes()),'id\tvalue\nprotein_α\t1\n'.encode('utf-8'))

    @unittest.skipUnless(shutil.which('git'), 'Git is required to test checkout rules')
    def test_git_autocrlf_preserves_scientific_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            repo=Path(temp)/'repo';repo.mkdir()
            subprocess.run(['git','init','-q',str(repo)],check=True,capture_output=True)
            for key,value in [('core.autocrlf','true'),('core.eol','crlf')]:
                subprocess.run(['git','-C',str(repo),'config',key,value],check=True)
            shutil.copy2(ROOT/'.gitattributes',repo/'.gitattributes')
            contents={'data/frozen.tsv':b'id\tlabel\r\nx\t1\r\n',
                      'results/frozen.tsv':b'id\tscore\r\nx\t0.5\r\n',
                      'configs/frozen.json':b'{"value":1}\r\n',
                      'examples/mini/frozen.tsv':b'id\tlabel\r\nx\t1\r\n',
                      'src/probe.py':b'value = 1\n'}
            for name,raw in contents.items():
                p=repo/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
            subprocess.run(['git','-C',str(repo),'add','.'],check=True,capture_output=True)
            checkout=Path(temp)/'checkout';checkout.mkdir()
            subprocess.run(['git','-C',str(repo),'checkout-index','--all','--prefix='+str(checkout)+'/'],check=True,capture_output=True)
            for name,raw in contents.items():
                with self.subTest(name=name):
                    self.assertEqual(subprocess.check_output(['git','-C',str(repo),'show',':'+name]),raw)
                    self.assertEqual((checkout/name).read_bytes(),raw)

    def test_invalid_utf8_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'invalid.tsv';path.write_bytes(b'id\n\xff\n')
            for reader in [read_table,read_data_table,read_tsv]:
                with self.subTest(reader=reader.__module__),self.assertRaises(UnicodeDecodeError):reader(path)

    def test_scoring_with_bom_crlf_and_utf8_mode_disabled(self):
        with tempfile.TemporaryDirectory(prefix='conformeval-io-') as temp:
            root=Path(temp)
            metadata=json.loads((ROOT/'examples/predictions/metadata.json').read_text(encoding='utf-8'))
            metadata['prediction_method']='external prediction at 25 °C'
            meta=root/'model.json';meta.write_bytes(b'\xef\xbb\xbf'+json.dumps(metadata,ensure_ascii=False).encode('utf-8'))
            env={**os.environ,'PYTHONUTF8':'0','PYTHONCOERCECLOCALE':'0','PYTHONIOENCODING':'ascii:backslashreplace','PYTHONWARNDEFAULTENCODING':'1','LC_ALL':'C'}
            for kind,filename in [('knot','knot.tsv'),('residue','residues.tsv')]:
                pred=root/filename;raw=(ROOT/'examples/predictions'/filename).read_bytes().replace(b'\r\n',b'\n').replace(b'\n',b'\r\n');pred.write_bytes(b'\xef\xbb\xbf'+raw)
                out=root/(kind+'-out')
                if kind=='knot':
                    args=['-m','conformeval.evaluate_predictions','--task','KNOT-PRESENCE','--data-dir',str(ROOT/'examples/mini/base'),'--protocol',str(ROOT/'configs/scoring_protocol.json')]
                else:
                    args=['-m','conformeval.evaluate_residue_predictions','--task-dir',str(ROOT/'examples/mini/residue'),'--example']
                command=[sys.executable,'-X','utf8=0',*args,'--predictions',str(pred),'--model-metadata',str(meta),'--output-dir',str(out)]
                run=subprocess.run(command,env=env,capture_output=True)
                self.assertEqual(run.returncode,0,run.stderr.decode('utf-8',errors='replace'))
                self.assertNotRegex(run.stderr.decode('utf-8',errors='replace'),r'[/\\]conformeval[/\\][^\r\n]*EncodingWarning')
                metrics=json.loads((out/'metrics.json').read_bytes().decode('utf-8'))
                if kind=='knot':self.assertEqual(metrics['metrics']['KNOT-PRESENCE']['auroc'],1.0)
                else:self.assertEqual(metrics[metadata['model_id']]['mean_per_protein_AP'],1.0)
                for name in ['metrics.json','run.json']:
                    payload=(out/name).read_bytes();payload.decode('utf-8');self.assertNotIn(b'\r\n',payload)


if __name__=='__main__':unittest.main()
