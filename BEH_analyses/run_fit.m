%% run_fit
% run winning model fit to human data
% FREE parameters:
%   gamma + siginf inside each condition
%   + epsi common to all conditions in all models
clearvars
close all
clc

addpath('..'); cfg_path = config_paths;
addpath(cfg_path.vbmc)

savepath = './results/';

% load data table
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = setdiff(unique(dat.participant),[]);
ndata = numel(participants);

out_fit = cell(ndata,1);

parfor idata = 1:ndata
    fprintf('Processing participant %02d/%02d ...\n',idata,ndata)
    cfg_fit = [];

    ifilt = dat.participant == participants(idata);

    cfg_fit.seqcat  = dat.seqcat(ifilt);
    cfg_fit.tltllr  = dat.LLR(ifilt);
    cfg_fit.resp    = dat.response(ifilt); % response => 1 or 2
    cfg_fit.condtn  = dat.condtn(ifilt);
    cfg_fit.verbose = false;
    cfg_fit.seed    = idata*200;

    % MODEL GAMMA+SIGINF
    fprintf('MODEL GAMMA+SIGINF\n')

    cfg_fit.sigsel1 = 0; cfg_fit.sigsel2 = 0; cfg_fit.sigsel3 = 0;


    out_fit{idata} = fit_model(cfg_fit);
end
save(fullfile(savepath,'results_fit.mat'),'out_fit');