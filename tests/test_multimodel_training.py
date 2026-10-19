"""End-to-end grouped training with small, explicitly synthetic data."""
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
from conformeval.readers import parse_fasta
from conformeval.runtime import digest, json_write
from conformeval.train import run, load_archives, reader_config, fit_reader, select_candidates, choose_candidate
from conformeval.evaluate_predictions import evaluate, TASK_KNOT, TASK_FS
from conformeval.evaluate_residue_predictions import TASK as RESIDUE_TASK, load_task, read_scores, calculate

ROOT=Path(__file__).resolve().parents[1]


def feature_fixture(root,dimension=4):
    root.mkdir()
    rng=np.random.default_rng(2026);rows=[]
    for sequence in parse_fasta(ROOT/'examples/sequences.fasta').values():
        sha=hashlib.sha256(sequence.encode('ascii')).hexdigest()
        residues=rng.normal(size=(2,len(sequence),dimension)).astype(np.float32)
        meta={'seq_sha256':sha,'len_full':len(sequence),'len_used':len(sequence),'truncated':False,
              'model':'synthetic-interface-test','revision':'a'*40,'layer_set':'mean_1_6+resid_1_6','pooling':'mean'}
        path=root/(sha+'.npz')
        np.savez_compressed(path,meta=np.asarray(json.dumps(meta)),sequence_layers=residues.mean(axis=1),
                            resid_layers=residues,resid_indices=np.arange(1,len(sequence)+1))
        rows.append({'sequence_sha256':sha,'length':len(sequence),'file':path.name,'archive_sha256':digest(path)})
    manifest={'model_id':'synthetic-interface-test','model_revision':'a'*40,'pooling':'mean',
              'output_dimension':dimension,'candidate_layers':[1,6],'extractor_code_sha256':'b'*64,'archives':rows}
    path=root/'representation_manifest.json';json_write(path,manifest)
    return path


class GenericTrainingTests(unittest.TestCase):
    def test_all_three_tasks_train_with_layers_below_seventeen(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);manifest=feature_fixture(root/'features')
            for task,source in [(TASK_KNOT,'base'),(TASK_FS,'base'),(RESIDUE_TASK,'residue')]:
                out=root/task
                args=SimpleNamespace(task=task,data_dir=ROOT/'examples'/source,feature_manifest=manifest,
                                     config=ROOT/'configs/readers/linear.json',protocol=ROOT/'configs/scoring.json',
                                     output_dir=out,example=True)
                actual=run(args)
                self.assertEqual(actual['status'],'PASS')
                self.assertEqual(actual['candidate_layers'],[1,6])
                self.assertFalse(actual['outer_labels_used_for_selection'])
                self.assertEqual(actual['scope'],'software_example')
                self.assertTrue((out/'predictions.tsv.gz').is_file())
                self.assertTrue((out/'parameters/PLM-nested_outer00.npz').is_file())
                self.assertTrue((out/'inner_selection.tsv').is_file())
                # Stored out-of-fold predictions must cover each held row exactly once.
                self.assertEqual(actual['scored_rows_per_arm'],554 if task==RESIDUE_TASK else 20)

    def test_model_dimensions_and_configuration_must_match_archives(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for dimension in [4,7]:
                path=feature_fixture(root/str(dimension),dimension)
                m=json.loads(path.read_text(encoding='utf-8-sig'));first=m['archives'][0];required={first['sequence_sha256']:first['length']}
                manifest,arrays,_=load_archives(path,required)
                self.assertEqual(arrays[first['sequence_sha256']].shape,(2,dimension))
                m['pooling']='max';json_write(path,m)
                with self.assertRaisesRegex(ValueError,'pooling differs'):load_archives(path,required)

    def test_reader_config_selects_ridge_and_scaling_uses_training_only(self):
        config=reader_config(ROOT/'configs/readers/linear.json')
        X=np.asarray([[0.,2.],[2.,4.],[1000.,-1000.]])
        y=np.asarray([0,1,0])
        for predictor in ['logistic','ridge']:
            scores,parameters=fit_reader(X,y,np.asarray([0,1]),np.asarray([2]),predictor,1.0,config,{},[])
            np.testing.assert_allclose(parameters['scaler_mean'],[1.,3.])
            self.assertTrue(np.isfinite(scores).all())
        for bad in [0,-1,float('nan')]:
            with tempfile.TemporaryDirectory() as temp:
                p=Path(temp)/'reader.json';json_write(p,{**config,'C_grid':[bad]}) if np.isfinite(bad) else p.write_bytes((json.dumps({**config,'C_grid':[bad]})).encode('utf-8'))
                with self.assertRaises(ValueError):reader_config(p)

    def test_external_prediction_examples_are_complete_and_score_one(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'score'
            evaluate(ROOT/'examples/base',ROOT/'examples/predictions/knot.tsv',ROOT/'examples/predictions/metadata.json',ROOT/'configs/scoring.json',out,{TASK_KNOT})
            self.assertEqual(json.loads((out/'metrics.json').read_text(encoding='utf-8-sig'))['metrics'][TASK_KNOT]['auroc'],1.0)
        samples,expected,folds,_=load_task(ROOT/'examples/residue',example=True)
        scores=read_scores(ROOT/'examples/predictions/residues.tsv',samples,expected,folds,{'example'},'example')
        _,metrics=calculate(samples,expected,folds,scores)
        self.assertEqual(metrics['example']['mean_per_protein_AP'],1.0)
        with self.assertRaisesRegex(ValueError,'unsupported disorder'):
            load_task(ROOT/'examples/residue')

    def test_outer_labels_cannot_change_inner_candidate_selection(self):
        rng=np.random.default_rng(0);X={1:rng.normal(size=(12,3)),6:rng.normal(size=(12,3))}
        y=np.tile([0,1],6);inners=[('0',np.asarray([0,1,2,3]),np.asarray([4,5,6,7]))]
        config=reader_config(ROOT/'configs/readers/linear.json')
        data={'task':TASK_KNOT,'residue':False,'y':y.copy()}
        first=select_candidates(data,X,inners,config,0,[],[])
        data['y'][8:]=1-data['y'][8:]
        second=select_candidates(data,X,inners,config,0,[],[])
        self.assertEqual(first,second)

    def test_residue_inner_aggregation_weights_proteins_and_honors_ties(self):
        config=reader_config(ROOT/'configs/readers/linear.json')
        config={**config,'C_grid':[0.1]}
        data={'task':RESIDUE_TASK,'residue':True,'y':np.zeros(10),
              'sample_index':np.asarray([0,0,0,0,1,1,2,2,3,3])}
        inners=[('0',np.arange(4,10),np.arange(4)),('1',np.arange(4),np.arange(4,10))]
        with patch('conformeval.train.fit_reader',return_value=(np.zeros(10),{})), \
             patch('conformeval.train.selection_value',side_effect=[0.1,0.9]):
            rows=[];candidates=select_candidates(data,{1:np.zeros((10,2))},inners,config,0,rows,[])
        self.assertAlmostEqual(candidates[0][3],0.7)
        self.assertEqual([r['aggregation_weight'] for r in rows],[1,3])
        config={**config,'C_grid':[0.1,1.0]}
        with patch('conformeval.train.fit_reader',return_value=(np.zeros(10),{})), \
             patch('conformeval.train.selection_value',side_effect=[0.8,0.8+5e-13]):
            candidates=select_candidates(data,{1:np.zeros((10,2))},inners[:1],config,0,[],[])
        self.assertEqual(choose_candidate(candidates,RESIDUE_TASK,1e-12)[2],0.1)


if __name__=='__main__':unittest.main()


class CandidateTieTests(unittest.TestCase):
    def test_ties_are_anchored_to_exact_best_not_chained(self):
        candidates=[(1,'logistic',0.1,0.8),(1,'logistic',1.0,0.8+0.75e-12),(1,'logistic',10.0,0.8+1.5e-12)]
        self.assertEqual(choose_candidate(candidates,RESIDUE_TASK,1e-12)[2],1.0)
        self.assertEqual(choose_candidate(candidates,TASK_FS,1e-12)[2],0.1)
