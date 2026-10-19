import argparse
import json
from pathlib import Path
from conformeval.adapters import HFAdapter, download_snapshot, load_config


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model-config',required=True,type=Path)
    args=p.parse_args();config=load_config(args.model_config)
    adapter=HFAdapter(config,download_snapshot(config))
    sequence='ACDEFGHIKLMNPQRSTVWY'
    pooled,residues,windows=adapter.encode(sequence)
    assert residues.shape[1]==len(sequence)
    print(json.dumps({'model':config['model_id'],'revision':config['revision'],
                      'sequence_shape':list(pooled.shape),'residue_shape':list(residues.shape),
                      'windows_0based_half_open':windows}))


if __name__=='__main__':main()
