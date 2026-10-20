# RNN analyses

Recurrent neural networks trained on the PRE and POST conditions of the task (Figure 7 and Supplementary
Figures 12 and 13). The networks are trained, tested and decoded in Python (`python/`); the statistics and figures
are computed in MATLAB from the results exported as `.mat` files. The results are not in the repository: run the
steps below in order.

## 0. Python environment (once)

```
conda env create -f python/environment.yml
conda activate NEURCOMP
```

## 1. Training: `python/scripts/training_PRE_POST.ipynb`

Trains the 100 networks (each one seeded with its own child of `SeedSequence(seed_base)`, so training is
reproducible) and saves them, every 20 epochs, in `python/rnn_simulations/RNN_train_PRE_POST/sigcat_15`.

## 2. Analyses: `python/scripts/analyses_PRE_POST.ipynb` and `analyses_PRE_POST_SUPP.ipynb`

Test the trained networks, decode the variables from their hidden activity, and save the results in the same folder:

| Notebook | Section | Figure | Exported file |
|---|---|---|---|
| `analyses_PRE_POST` | 1. along training | 7b | `test_epochs.mat` |
| `analyses_PRE_POST` | 2. resilience to computation noise | 7c | `test_noise.mat` |
| `analyses_PRE_POST` | 3. causal perturbation | 7d, Supp. 13 | `test_perturbation.mat` |
| `analyses_PRE_POST_SUPP` | S1. resilience to input noise | Supp. 12 | `test_input_noise.mat` |

## 3. Statistics and figures (MATLAB)

From this folder (`RNN_analyses`), in MATLAB R2022b, with the paths set in `config_paths.m` (in the folder above):
`analyses_figure_7.m` (Figure 7; also reads `data/data_table_seq.mat` for the human accuracy) and
`analyses_supp_figures_RNN.m` (Supplementary Figures 11a, 12 and 13; Supp. Fig. 11a simulates the computational
model with the behavioral fits `BEH_analyses/results/results_fit.mat` and also reads `data/data_table_seq.mat`).
