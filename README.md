# NEURCOMP

Analysis code of the article:

Flexible compression of sensory evidence makes human decisions more resilient to noisy neural computations.
Drevet, J., Benjamin, L., Rondot, L., Drugowitsch, J. & Wyart, V.  

The code reproduces the statistics and figures of the article from the behavioral data (in this repository), the
preprocessed MEG data (Zenodo: https://doi.org/10.5281/zenodo.23266592) and recurrent neural networks trained with the provided code.

## Content

| Folder | Content | Figures |
|---|---|---|
| `data/` | behavioral data of the 37 participants (MATLAB tables): one row per sequence (`data_table_seq.mat`), and the presentations of each sequence in all conditions and repetitions aligned (`data_table_cond.mat`) | |
| `BEH_analyses/` | behavioral analyses, model fitting, model selection, model and parameter recovery, simulations; the fits are in `results/` (see its `README.md`) | 2, 3, Supp. 2–5 |
| `MEG_analyses/` | epoching, sensor- and source-level decoding of the MEG data (see its `README.md`) | 4–6, Supp. 6–10 |
| `RNN_analyses/` | training and analyses of the recurrent neural networks in Python, statistics and figures in MATLAB (see its `README.md`) | 7, Supp. 11a, 12, 13 |
| `functions/` | functions shared by the analyses (repeated-measures ANOVA, cluster-based permutation test, ridge decoding, jackknife, ...) | |
| `helper_plots/` | plotting functions | |
| `config_paths.m` | the only place where paths are set | |

## Requirements

- MATLAB R2022b, with the Statistics and Machine Learning, Signal Processing, Image Processing, Curve Fitting and
  Parallel Computing toolboxes.
- External toolboxes, not included in this repository (set their paths in `config_paths.m`):
  - [bayesFactor](https://doi.org/10.5281/zenodo.4394422) 1.1 (Bayes factors),
  - [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Bayesian model selection),
  - [VBMC](https://github.com/acerbilab/vbmc) 1.0.12 (model fitting),
  - [FieldTrip](https://www.fieldtriptoolbox.org) 20240115 (epoching),
  - [Brainstorm](https://neuroimage.usc.edu/brainstorm) (source-level maps).
- Python 3.11 with PyTorch 2.5.1 (and NumPy, SciPy, Matplotlib, Joblib, Jupyter) for the recurrent neural networks:
  conda environment `RNN_analyses/python/environment.yml` (PyTorch runs on the CPU, or on a CUDA GPU when available).

## Setup

1. Set the paths of the external toolboxes and of the MEG data, either in `config_paths.m` or in a file
   `config_paths_local.m` (not shared, see the header of `config_paths.m`).
2. For the MEG analyses, download the preprocessed MEG data of each participant from Zenodo and unzip the folders
   `SXX` into `cfg_path.megcont`. The source-level decoding results (`zdec_source.zip`, Zenodo) are unzipped
   into `MEG_analyses/results`.

## Reproducing the figures

Run each script from its folder. The scripts print the statistics reported in the article; set `saveornot = true`
to save the figures (PDF) in the folder `figs`.

| Figure | Script |
|---|---|
| 2 | `BEH_analyses/analyses_figure_2_behavior.m` |
| 3 | `BEH_analyses/analyses_figure_3_model.m` |
| 4, Supp. 6 | `MEG_analyses/analyses_figure_4.m` |
| 5, Supp. 7, 8 | `MEG_analyses/analyses_figure_5.m` |
| 6, Supp. 9, 10 | `MEG_analyses/analyses_figure_6.m` |
| 7 | `RNN_analyses/analyses_figure_7.m` |
| Supp. 2 | `BEH_analyses/analyses_supp_RT.m` |
| Supp. 3–5 | `BEH_analyses/analyses_supp.m` |
| Supp. 11a, 12, 13 | `RNN_analyses/analyses_supp_figures_RNN.m` |

Figure 1 and Supplementary Figures 1 and 11b are schematics.

The model fits and simulations (`BEH_analyses/results`) and the sensor-level decoding results
(`MEG_analyses/results/zdec_sensor`) are included, and can be recomputed with the `run_*.m` scripts of each folder
(long). The networks are not included: train and analyze them first with the Python notebooks (see
`RNN_analyses/README.md`).

## Licence

The code is distributed under the BSD 3-Clause License (see `LICENSE`), except `MEG_analyses/default_cortex.mat`
(below).

The code was cleaned up and documented for release with the help of Claude (Anthropic); all analyses were checked by
the authors.

`MEG_analyses/default_cortex.mat` is the cortical surface of the ICBM152 template (Brainstorm default anatomy), with
the Desikan-Killiany atlas (Desikan, R. S. et al. NeuroImage 31, 968–980, 2006). ICBM152 template:

> Copyright (C) 1993–2004 Louis Collins, McConnell Brain Imaging Centre, Montreal Neurological Institute, McGill
> University. Permission to use, copy, modify, and distribute this software and its documentation for any purpose
> and without fee is hereby granted, provided that the above copyright notice appear in all copies. The authors and
> McGill University make no representations about the suitability of this software for any purpose. It is provided
> "as is" without express or implied warranty. The authors are not responsible for any data loss, equipment damage,
> property loss, or injury to subjects or patients resulting from the use or misuse of this software package.
