# Epoching and decoding of the MEG data

How to recompute the decoding results read by the figure scripts (`results/zdec_sensor` and `results/zdec_source`)
from the preprocessed MEG data. All scripts are run from this folder (`MEG_analyses`), in MATLAB R2022b, with the
paths set in `config_paths.m` (in the folder above). The epochs are saved with the MEG data (`cfg_path.megepo`),
the decoding results in `MEG_analyses/results`.

## 0. Data and paths (once)

- Download the preprocessed MEG data from Zenodo and unzip the participant folders `SXX` (continuous data,
  events and source model of each participant) into one folder: `cfg_path.megcont` in `config_paths.m`
  (by default the folder `continuous` in `cfg_path.megdata`).
- Set the paths of FieldTrip (epoching) and Brainstorm in `config_paths.m` (or in `config_paths_local.m`, see
  `config_paths.m`).

## 1. Epochs: `run_epoching.m`

Cuts the epochs of all participants (`subjlist`) from the continuous data, low-pass filtered at 16 Hz (stimulus
orientation) and at 8 Hz (other variables), with the Hilbert transform (real and imaginary parts), and saves them as
binary files in `cfg_path.megepo/SXX/epochs_16Hz_bin` and `epochs_08Hz_bin` (one pair of files per condition, with
`cfg_condC.mat`). For one participant: `epoching_meg_from_continuous(isubj,freq)`.

## 2. Sensor-level decoding: `run_decoding_sensor.m`

Checks first that the epochs of all participants exist (otherwise lists the missing files and stops), then decodes,
for each participant, in each condition and at each time sample (-0.2 to 0.8 s around the stimulus onset, 200 Hz), the variables below
from the real and imaginary parts of the planar gradiometers (ridge decoding, `decode_ridg`, 10 folds):

| `varnam` | Epochs | Decoded | Result files |
|---|---|---|---|
| `'seqang'` | 16 Hz | stimulus orientation (cos and sin of the doubled angle) | `seqang` |
| `'absllr'` | 8 Hz | evidence strength \|LLR\| | `absllr` |
| `'conllr'` | 8 Hz | belief updates along the true axis and along the null axis (particle filter with the best-fitting model parameters, `BEH_analyses/results/results_fit.mat`) | `conllr_particle`, `nulllr_particle` |

For one participant: `decoding_sensor(isubj,varnam,freq,nblock)` (`freq` optional: 16 Hz for `seqang`, 8 Hz
otherwise; `nblock` optional, see Memory below).

Results, in `results/zdec_sensor`:
- `tmp/SXX_decoding_<name>_condC_FFHz.mat`: decoding precision (Fisher z of the correlation) of each participant and
  condition, `zdec` (1 x time sample, x cos/sin for `seqang`), and `_optsub.mat` (same x optimal/suboptimal trials);
- `zdec_<name>_FFHz.mat` and `_optsub.mat`: group files (participant x condition), made by
  `collect_decoding_sensor(name,freq,subjlist)` at the end of each section of the script, and read by the figure
  scripts (already in the repository).

## 3. Source-level decoding: `run_decoding_source.m`

Checks first that the source models and the epochs of all participants exist, then decodes the same variables
(`decoding_source(isubj,varnam,freq,nblock)`) at each vertex of the participant's cortex: all MEG channels projected on the
sources with the inverse kernel, and decoding from a searchlight (the vertex and its nearest neighbours on the
inflated cortex, 10 vertices; `decode_source_searchlight.m`). This step is long (about 1 hour per participant for
`conllr` with 12 parallel workers): set `subjlist` in the script to split it into several runs, and `parpool(n)` to
fewer workers if the memory is short.

Results, in `results/zdec_source`:
- `tmp/SXX_decoding_<name>_condC_FFHz.mat` and `_optsub.mat`: `zdec` (vertex x time sample, on the participant's
  cortex), about 13 GB for all participants and variables;
- `zdec_<name>_FFHz.mat` and `_optsub.mat`: group files, projected on the template cortex (`default_cortex.mat`) and
  gathered by `collect_decoding_source(name,freq,0:36)` at the end of each section of the script (once all
  participants are done), and read by the figure scripts (several GB, available on Zenodo).

## Memory

By default, the epochs of a condition are read whole (`read_epochs.m`: two binary files of about 0.7 GB each, real
and imaginary parts). Count about 3 GB of memory for the sensor-level decoding plus 1 to 2 GB per parallel worker,
and about 4 GB plus 3.5 GB per worker for the source level. With less memory:
- start fewer workers (`parpool(n)` in the scripts, before the decoding);
- and/or read and decode the epochs by blocks of time samples (`nblock` in the scripts, e.g. 20: about 10 times less
  memory for the epochs, a bit slower, same results).

## Functions

| Function | Role |
|---|---|
| `read_events.m` | stimulus events of a participant (`SXX_events.tsv`) |
| `read_epochs.m` | epochs of a condition (binary files of step 1) |
| `get_decoding_variables.m` | variables decoded (`absllr`, `conllr`, `seqang`), including the particle filter |
| `decode_ridg.m` (folder `functions`) | ridge decoding with cross-validation |
| `load_source_model.m` | inverse kernel and searchlight neighbours of a participant (`SXX_source_model.mat`) |
| `project_source_default.m` | projection of a participant's map on the template cortex |
