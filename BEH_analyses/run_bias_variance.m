function out_sim = run_bias_variance
% simulation of a model with inference noise, varying the fraction of the 
% noise variance shared between two presentations of the same sequence 
% (rsq, from 0 = independent to 1 = identical noise)

addpath('..'); cfg_path = config_paths;

savepath = './results';

% load data table (sequence aligned)
load(fullfile(cfg_path.data,'data_table_cond.mat'),'dat');

participants = setdiff(unique(dat.participant),[]);
ndata = numel(participants);

% load output from model with only SIGINF
load(fullfile('./results','results_fit_choice_variability.mat'),'out_fit');
siginf = [cellfun(@(s)s.siginf1,out_fit),...
    cellfun(@(s)s.siginf2,out_fit),...
    cellfun(@(s)s.siginf3,out_fit)];

rsq = 0:.01:1;
nrsq = numel(rsq);

out_sim = cell(ndata,nrsq);

fprintf('Computing consistency on repeated trial increasing inference noise correlation\n')
for idata = 1:ndata
    fprintf('\tparticipant %02d/%02d...\n',idata,ndata)
    ifilt = dat.participant == participants(idata);
    cfg_sim = [];
    cfg_sim.seqllr = dat.seqllr(ifilt);
    cfg_sim.seqlen = dat.seqlen(ifilt);
    cfg_sim.siginf = siginf(idata,:);
    cfg_sim.seed   = idata * 100;
    for irsq = 1:nrsq
        cfg = cfg_sim;
        cfg.rsq = rsq(irsq);
        out_sim{idata,irsq} = sim_model_rsq(cfg);
    end
end
fprintf('done.\n')

fname = 'results_sim_rsq.mat';
save(fullfile(savepath,fname),'out_sim');

end