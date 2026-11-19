# NEURCOMP

Analysis code of the article:

**Flexible compression of sensory evidence makes human decisions more resilient to noisy neural computations.**\
Drevet, J., Benjamin, L., Rondot, L., Drugowitsch, J. & Wyart, V.

The code reproduces all analyses of the article. The preprocessed MEG data are available on Zenodo:
https://doi.org/10.5281/zenodo.23266592.

## Content

| Folder | Content | Figures |
|---|---|---|
| `data/` | behavioral data of the 37 participants | |
| `BEH_analyses/` | behavioral analyses and computational model | 2, 3, Supp. 2–5 |
| `MEG_analyses/` | MEG decoding | 4–6, Supp. 6–10 |
| `RNN_analyses/` | recurrent neural networks | 7, Supp. 11a, 12, 13 |
| `functions/` | shared analysis functions | |
| `helper_plots/` | plotting functions | |
| `config_paths.m` | paths to the data and toolboxes | |

## System requirements

- Tested on macOS 26.4 (Apple M3 Max, 48 GB RAM).
- The decoding of the MEG data needs several GB of memory per parallel worker (see `MEG_analyses/README.md`).
  PyTorch runs on the CPU, or on a CUDA GPU when available.
- MATLAB R2022b, with the Statistics and Machine Learning, Signal Processing, Image Processing, Curve Fitting and
  Parallel Computing toolboxes.
- External toolboxes, not included in this repository:
  - [bayesFactor](https://doi.org/10.5281/zenodo.4394422) 3.0 (Bayes factors),
  - [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) release 7487 (Bayesian model selection),
  - [VBMC](https://github.com/acerbilab/vbmc) 1.0.12 (model fitting),
  - [FieldTrip](https://www.fieldtriptoolbox.org) 20240115 (epoching),
  - [Brainstorm](https://neuroimage.usc.edu/brainstorm) 3.250603 (source-level maps).
- Python 3.11 with PyTorch 2.5.1 (and NumPy, SciPy, Matplotlib, Joblib, Jupyter) for the recurrent neural networks
  (exact versions in `RNN_analyses/python/environment.yml`).

## Installation

1. Clone the repository: `git clone https://github.com/juliedrevet/NEURCOMP.git`
2. Download the external toolboxes, then set their paths and the path of the MEG data, either in `config_paths.m`
   or in a file `config_paths_local.m` (not shared, see the header of `config_paths.m`).
3. For the recurrent neural networks, create the Python environment:
   `conda env create -f RNN_analyses/python/environment.yml`
4. To re-run the MEG decoding, download the preprocessed MEG data of each participant from Zenodo and unzip the
   folders `SXX` into `cfg_path.megcont`. To plot the source-level maps, unzip the source-level decoding results
   (`zdec_source.zip`, Zenodo) into `MEG_analyses/results`.

Typical install time: about 20 minutes (excluding MATLAB), plus the download of the MEG data from Zenodo
(optional, 38 GB, about 1 hour at 100 Mbit/s).

## Demo and reproduction of the figures

All figures of the article can be reproduced from the results included in the repository. Run each script below
from its folder: it plots the figures and prints the statistics of the article (under 1 minute each).

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

Figure 1 and Supplementary Figures 1 and 11b are schematics. The source-level maps (`plotsource = true` in
`analyses_figure_4.m` and `analyses_figure_5.m`) need the source-level decoding results from Zenodo.

## Instructions for use

The model fits and simulations (`BEH_analyses/results`), the sensor-level decoding results
(`MEG_analyses/results/zdec_sensor`) and the results of the recurrent neural networks
(`RNN_analyses/python/rnn_simulations`) are included. To recompute them
from the data, follow the README of `BEH_analyses`, then `MEG_analyses`, then `RNN_analyses`. The longest steps take
up to about 1.5 hours per participant.

## Licence

The code is distributed under the BSD 3-Clause License (see `LICENSE`), except `MEG_analyses/default_cortex.mat`
(below).

The code was cleaned up and documented for release with the help of Claude Opus 5.5 (Anthropic). All analyses were
checked by Julie Drevet.

`MEG_analyses/default_cortex.mat` is the cortical surface of the ICBM152 template (Brainstorm default anatomy), with
the Desikan-Killiany atlas (Desikan, R. S. et al. NeuroImage 31, 968–980, 2006). ICBM152 template:

> Copyright (C) 1993–2004 Louis Collins, McConnell Brain Imaging Centre, Montreal Neurological Institute, McGill
> University. Permission to use, copy, modify, and distribute this software and its documentation for any purpose
> and without fee is hereby granted, provided that the above copyright notice appear in all copies. The authors and
> McGill University make no representations about the suitability of this software for any purpose. It is provided
> "as is" without express or implied warranty. The authors are not responsible for any data loss, equipment damage,
> property loss, or injury to subjects or patients resulting from the use or misuse of this software package.
