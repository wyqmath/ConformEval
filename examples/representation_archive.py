"""Write complete deterministic synthetic archives for the mini task fixture."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from confomeval.readers import parse_fasta
from confomeval.runtime import digest, json_write


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True,type=Path)
    args=p.parse_args();output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=False)
    layers=[1,2];dimension=8;revision='0'*40;rng=np.random.default_rng(2026);archives=[]
    source=Path(__file__).parent/'mini/sequences.fasta'
    for sequence in parse_fasta(source).values():
        sha=hashlib.sha256(sequence.encode('ascii')).hexdigest()
        values=rng.normal(size=(len(layers),len(sequence),dimension)).astype(np.float32)
        meta={'seq_sha256':sha,'len_full':len(sequence),'len_used':len(sequence),'truncated':False,
              'model':'synthetic-schema-example','revision':revision,'layer_set':'mean_1_2+resid_1_2','pooling':'mean'}
        path=output/(sha+'.npz')
        # Real adapters supply encoder vectors at the same full-sequence coordinates.
        np.savez_compressed(path,meta=np.asarray(json.dumps(meta)),sequence_layers=values.mean(axis=1),
                            resid_layers=values,resid_indices=np.arange(1,len(sequence)+1))
        archives.append({'sequence_sha256':sha,'length':len(sequence),'file':path.name,'archive_sha256':digest(path)})
    json_write(output/'representation_manifest.json',{'schema_version':1,'model_id':'synthetic-schema-example',
               'model_revision':revision,'candidate_layers':layers,'output_dimension':dimension,'pooling':'mean',
               'extractor_code_sha256':digest(Path(__file__)),'archives':archives})
    print(json.dumps({'manifest':str(output/'representation_manifest.json'),'unique_sequences':len(archives),'dimension':dimension}))


if __name__=='__main__':main()
