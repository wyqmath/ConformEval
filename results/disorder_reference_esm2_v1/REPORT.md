# DISORDER-STRUCTURE-REFERENCE-CANDIDATE nested residue evaluation

Completed: 2026-10-15T22:50:18+08:00 (Asia/Shanghai).

This is a separately versioned development evaluation. Label 1 is direct DisProt disorder; label 0 combines direct DisProt order and amino-acid-matched PDB-observed reference positions. It is not strict direct ordered/disordered truth and does not establish family-independent generalization.

The frozen task contains 1279 proteins in 1029 binding groups. It scores 391,718 labeled protein-residue rows over 1,279 proteins; unknown positions remain excluded. All outer/inner roles are consumed from the saved grouped fold tables.

| Arm | Mean protein AP (primary) | Mean protein AUROC | Pooled residue AP (secondary) |
|---|---:|---:|---:|
| ESM2-L11 | 0.722277 | 0.897205 | 0.815593 |
| ESM2-L23 | 0.732330 | 0.900249 | 0.822946 |
| ESM2-L33 | 0.789856 | 0.921073 | 0.849966 |
| traditional-AA-context | 0.585906 | 0.822784 | 0.624381 |
| MLP-h33-seed13 | 0.802670 | 0.926454 | 0.861817 |
| MLP-h33-seed42 | 0.802767 | 0.926744 | 0.861738 |
| MLP-h33-seed2026 | 0.803926 | 0.926312 | 0.861194 |

95% intervals are binding-group bootstrap intervals conditional on the fixed outer predictions. MLP seeds are all reported; no seed is selected by outer performance.

The traditional baseline uses amino-acid one-hot, relative position, and centered 15-residue composition. All scaling is fit on training residues only. The MLP capacity diagnostic uses hidden layer 33 and inner-fold early stopping.
