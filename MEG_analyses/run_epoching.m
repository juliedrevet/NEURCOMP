%% run_epoching
% epochs of all participants from the preprocessed continuous MEG data (SXX_meg_continuous.mat, Zenodo),
% needed by the decoding (run_decoding_sensor.m and run_decoding_source.m): low-pass filtered at 16 Hz
% (orientation, seqang) and at 8 Hz (other variables), saved in cfg_path.megepo/SXX/epochs_16Hz_bin and
% epochs_08Hz_bin (epoching_meg_from_continuous.m)

clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

subjlist = 0:36;

for isubj = subjlist
    epoching_meg_from_continuous(isubj,16)
    epoching_meg_from_continuous(isubj,8)
end
