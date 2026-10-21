function cfg_path = config_paths()
% config_paths  paths of the NEURCOMP analyses: the ONLY place where paths are set
%   cfg_path = config_paths adds the code and the toolboxes used by all the analyses to the MATLAB path
%   and returns a structure with the paths of:
%       root         this folder (NEURCOMP)
%       functions    helper functions shared by the analyses (functions/)
%       data         behavioral data tables (data/data_table_seq.mat and data/data_table_cond.mat)
%       results      model fits and simulations (BEH_analyses/results)
%       rnnsim       results of the RNN simulations saved by the Python notebooks (RNN_analyses/python/rnn_simulations)
%       megdata      root folder of the MEG data: continuous/ and epochs/
%       megcont      preprocessed MEG data downloaded from Zenodo: one folder per participant (SXX) with
%                    the continuous data, the events and the source model
%       megepo       epochs computed from the preprocessed MEG data (binary files for the decoding)
%       bayesfactor  bayesFactor toolbox (Bayes factors)
%       spm12        SPM12 (Bayesian model selection: spm_BMS)
%       vbmc         VBMC toolbox (model fitting)
%       fieldtrip    FieldTrip toolbox (epoching)
%       brainstorm   Brainstorm (source-level brain maps)
%   Set the paths marked [local] to your own installation, either below or in a file config_paths_local.m
%   in this folder (not shared: see .gitignore), which then replaces them, e.g.:
%       function cfg_path = config_paths_local(cfg_path)
%       cfg_path.megdata     = '/my/MEG/data';
%       cfg_path.bayesfactor = '/my/toolboxes/bayesFactor';
%       end

cfg_path = [];
cfg_path.root      = fileparts(mfilename('fullpath'));
cfg_path.functions = fullfile(cfg_path.root,'functions');
cfg_path.data      = fullfile(cfg_path.root,'data');
cfg_path.results   = fullfile(cfg_path.root,'BEH_analyses','results');
cfg_path.rnnsim    = fullfile(cfg_path.root,'RNN_analyses','python','rnn_simulations'); % RNN simulations (Python notebooks)

% [local] MEG data: unzip the folders SXX downloaded from Zenodo in cfg_path.megcont (continuous/)
cfg_path.megdata = '/path/to/MEG_data';

% [local] external toolboxes
cfg_path.bayesfactor = '/path/to/bayesFactor';        % bayesFactor 1.1 (MATLAB File Exchange)
cfg_path.spm12       = '/path/to/spm12';              % SPM12
cfg_path.vbmc        = '/path/to/vbmc1.0.12';         % VBMC 1.0.12 (added by the fitting scripts)
cfg_path.fieldtrip   = '/path/to/fieldtrip-20240115'; % FieldTrip 20240115 (added by the epoching script)
cfg_path.brainstorm  = '/path/to/brainstorm3';        % Brainstorm (added by the brain map sections)

% [local] paths of this computer, if config_paths_local.m exists
if exist(fullfile(cfg_path.root,'config_paths_local.m'),'file')
    cfg_path = config_paths_local(cfg_path);
end

% MEG data subfolders
cfg_path.megcont = fullfile(cfg_path.megdata,'continuous');
cfg_path.megepo  = fullfile(cfg_path.megdata,'epochs');

% code and toolboxes used by all the analyses
addpath(cfg_path.root)
addpath(cfg_path.functions)
addpath(fullfile(cfg_path.root,'helper_plots'))
addpath(cfg_path.bayesfactor)
addpath(cfg_path.spm12)

end
