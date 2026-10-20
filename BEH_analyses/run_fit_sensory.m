function out_fit = run_fit_sensory 
% fit the sensory-based model with GAMMA+SIGINF+EPSI

addpath('..'); cfg_path = config_paths; 
addpath(cfg_path.vbmc)

savepath = './results/';

% load data table
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = setdiff(unique(dat.participant),[]);
ndata = numel(participants);

out_fit = cell(ndata,1);
fname = 'results_fit_sensory.mat'; 

parfor idata = 1:ndata
    fprintf('Processing participant %02d/%02d ...\n',idata,ndata)
    cfg = [];
    ifilt = dat.participant == participants(idata);
    cfg.seqcat = dat.seqcat(ifilt);
    cfg.catang = dat.catang(ifilt);
    cfg.seqang = dat.seqang(ifilt);
    cfg.resp   = dat.response(ifilt); % response => 1 or 2
    cfg.condtn = dat.condtn(ifilt);
    cfg.verbose = false;
    cfg.seed    = idata*200;

    %% MODEL GAMMA+SIGINF
    fprintf('MODEL GAMMA+SIGINF\n')
    cfg_fit = cfg;
    cfg_fit.sigsel1 = 0; cfg_fit.sigsel2 = 0; cfg_fit.sigsel3 = 0;

    out_fit{idata} = fit_model_sensory(cfg_fit);

    fprintf('done.\n')
end
save(fullfile(savepath,fname),'out_fit');