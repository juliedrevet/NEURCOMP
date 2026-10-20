function out_fit = run_fit_choice_variability(idata)
% model with only the inference noise siginf as free 
% siginf1, siginf2, siginf3: one per condition
% all the other parameters set to 0: gamma, sigsel, epsi

addpath('..'); cfg_path = config_paths;
addpath(cfg_path.vbmc)

savepath = './results/tmp';
[~,~] = mkdir(savepath);

% load data table
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = setdiff(unique(dat.participant),[]);
ndata = numel(participants);

fprintf('Processing participant %02d/%02d ...\n',idata,ndata)

fname = sprintf('results_fit_choice_variability_S%02d.mat',participants(idata));

ifilt = dat.participant == participants(idata);
cfg.seqcat  = dat.seqcat(ifilt);
cfg.tltllr  = dat.LLR(ifilt);
cfg.resp    = dat.response(ifilt); % response => 1 or 2
cfg.condtn  = dat.condtn(ifilt);
cfg.verbose = false; % or true if subject loop not run in parallel
cfg.seed    = idata*100;

%% MODEL SIGINF ONLY (NO EPSI)
fprintf('MODEL SIGINF\n')
cfg_fit = cfg;
cfg_fit.gamma1  = 0; cfg_fit.gamma2  = 0; cfg_fit.gamma3  = 0;
cfg_fit.sigsel1 = 0; cfg_fit.sigsel2 = 0; cfg_fit.sigsel3 = 0;
cfg_fit.epsi = 0;

out_fit = fit_model(cfg_fit);
save(fullfile(savepath,fname),'out_fit');
fprintf('done.\n')




