"""Adapter projection and window checks without downloading model weights."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from confomeval.adapters import HFAdapter, load_config, token_residue_map, windows

ROOT=Path(__file__).resolve().parents[1]


class Tokenizer:
    def convert_ids_to_tokens(self, ids):
        return [{0:'<pad>',1:'<cls>',2:'<eos>',3:'A',4:'▁C',5:'D',6:'[MASK]'}[i] for i in ids]


class AdapterTests(unittest.TestCase):
    def test_special_tokens_and_padding_never_become_residue_features(self):
        tokenizer=Tokenizer()
        got=token_residue_map([1,3,4,5,2,0],[1,1,1,1,1,0],[1,0,0,0,1,1],tokenizer,'ACD')
        np.testing.assert_array_equal(got,[1,2,3])
        with self.assertRaisesRegex(ValueError,'identity mismatch'):
            token_residue_map([1,3,4,2],[1]*4,[1,0,0,1],tokenizer,'ACD')
        with self.assertRaisesRegex(ValueError,'non-residue'):
            token_residue_map([6],[1],[0],tokenizer,'A')

    def test_window_stitching_preserves_all_coordinates_and_pooling(self):
        adapter=object.__new__(HFAdapter)
        adapter.config={'layers':[1,6],'pooling':'mean','long_sequence':{'policy':'window','max_residues':4,'overlap':2}}
        adapter.dimension=2
        def window(sequence):
            a=np.asarray([[ord(c),ord(c)+1] for c in sequence],dtype=np.float32)
            return np.stack([a,a+10])
        adapter.encode_window=window
        sequence='ACDEFGHIK'
        mean,residues,segments=adapter.encode(sequence)
        np.testing.assert_allclose(residues,window(sequence))
        np.testing.assert_allclose(mean,window(sequence).mean(axis=1))
        self.assertEqual(segments,[(0,4),(2,6),(4,8),(5,9)])
        adapter.config['pooling']='max'
        maximum,_,_=adapter.encode(sequence)
        np.testing.assert_allclose(maximum,window(sequence).max(axis=1))
        with self.assertRaisesRegex(ValueError,'exceeds'):
            windows(9,4,2,'reject')

    def test_model_configuration_accepts_6_24_and_33_layer_models(self):
        for name,layers in [('esm2_8m',[3,6]),('esm2_650m',[17,33]),('prott5',[12,24])]:
            config=load_config(ROOT/'configs/models'/f'{name}.json')
            self.assertEqual(config['layers'],layers)
        original=json.loads((ROOT/'configs/models/esm2_8m.json').read_text(encoding='utf-8-sig'))
        for field,value in [('layers',[0]),('revision','main'),('pooling','include_special_tokens')]:
            with self.subTest(field=field),tempfile.TemporaryDirectory() as temp:
                invalid={**original,field:value};path=Path(temp)/'model.json';path.write_text(json.dumps(invalid), encoding='utf-8', newline='\n')
                with self.assertRaises(ValueError):load_config(path)


if __name__=='__main__':unittest.main()
