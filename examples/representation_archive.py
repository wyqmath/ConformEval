import argparse
import json
from pathlib import Path
import numpy as np
from conformeval.datasets import parse_fasta
from conformeval.runtime import json_write


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-dir',required=True,type=Path)
    args=p.parse_args();output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=False)
    layers=[1,2];dimension=8;revision='example';rng=np.random.default_rng(2026);archives=[]
    source=Path(__file__).parent/'sequences.fasta'
    for index,(sequence_id,sequence) in enumerate(parse_fasta(source).items(),start=1):
        values=rng.normal(size=(len(layers),len(sequence),dimension)).astype(np.float32)
        meta={'sequence_id':sequence_id,'len_full':len(sequence),'len_used':len(sequence),'truncated':False,
              'model':'synthetic-schema-example','revision':revision,'layer_set':'mean_1_2+resid_1_2','pooling':'mean'}
        path=output/f'sequence-{index:06d}.npz'
        # Real adapters supply encoder vectors at the same full-sequence coordinates.
        np.savez_compressed(path,meta=np.asarray(json.dumps(meta)),sequence_layers=values.mean(axis=1),
                            resid_layers=values,resid_indices=np.arange(1,len(sequence)+1))
        archives.append({'sequence_id':sequence_id,'length':len(sequence),'file':path.name})
    json_write(output/'representation_manifest.json',{'schema_version':1,'model_id':'synthetic-schema-example',
               'model_revision':revision,'candidate_layers':layers,'output_dimension':dimension,'pooling':'mean',
               'archives':archives})
    print(json.dumps({'manifest':str(output/'representation_manifest.json'),'unique_sequences':len(archives),'dimension':dimension}))


if __name__=='__main__':main()
