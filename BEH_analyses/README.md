# Behavioral analyses and computational model

The figure scripts (`analyses_*.m`, see the root `README.md`) read the model fits and simulations saved in
`results`, which are included in the repository. This page lists how to recompute them. All scripts are run from
this folder (`BEH_analyses`), in MATLAB R2022b, with the paths set in `config_paths.m` (in the folder above); the
fits use VBMC (`cfg_path.vbmc`) and take long (participants fitted in parallel with `parfor`).

## Fits and simulations

| Script | Result in `results` | Used by |
|---|---|---|
| `run_fit.m` | `results_fit.mat`: winning model (inference imbalance and inference noise in each condition, common lapse) | Figure 3, Supp. Figures 3–5, MEG decoding (belief updates) and Figure 6, Supp. Figure 11a |
| `run_fit_cond.m` (one participant, called and gathered by `analyses_figure_3_model.m`) | `results_fit_cond.mat`: the 8 model variants of the factorized Bayesian model selection | Figure 3b |
| `run_fit_choice_variability.m` (one participant, called and gathered by `analyses_figure_2_behavior.m`) | `results_fit_choice_variability.mat`: model with inference noise only | Figure 2c |
| `run_bias_variance.m` | `results_sim_rsq.mat`: simulations of the bias-variance decomposition (reads `results_fit_choice_variability.mat`) | Figure 2e |
| `run_fit_sensory.m` | `results_fit_sensory.mat`: model accumulating evidence in sensory space | Supp. Figure 3 |
| `run_rec_cond.m`, `fit_rec_cond.m`, `collect_rec_cond.m` | `results_rec_cond.mat`: model and parameter recovery | Supp. Figure 4a |

The figure scripts run the missing fits themselves (`results_fit.mat`, `results_fit_cond.mat`,
`results_fit_choice_variability.mat`, `results_sim_rsq.mat`); the per-participant files are saved in `results/tmp`.

## Model and parameter recovery

1. `run_rec_cond.m` simulates the 8 model variants for each participant (parameters drawn from the prior
   distributions, saved in `results/tmp/recovery`), then fits each simulated dataset with each model variant: one
   fit per job, `fit_rec_cond(irec,irep)`, run in parallel by the script or as a job array on a cluster
   (37 × 8 × 8 fits).
2. `collect_rec_cond.m` gathers the fits into `results_rec_cond.mat`, read by `analyses_supp.m`.

## Functions

| Function | Role |
|---|---|
| `fit_model.m` | fit of the model to the choices of a participant (particle filter and VBMC) |
| `fit_model_sensory.m` | same for the model accumulating evidence in sensory space |
| `sim_model.m` | simulation of the model |
| `sim_model_rsq.m` | simulation with inference noise shared between the two presentations of a sequence |
