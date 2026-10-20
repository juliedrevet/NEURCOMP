%% collect_rec_cond
% collects the model and parameter recovery (run_rec_cond.m and fit_rec_cond.m: one file per fit in
% ./results/tmp/recovery, about 110 MB in all) into one light file, read by analyses_supp.m:
% ./results/results_rec_cond.mat with
%   elbo          participant x simulated model x fitted model: ELBO of each fit (model recovery)
%   par_sim       participant x parameter: parameters of the simulated full model (model 1), named in parnam_sim
%   xnam          names of the fitted parameters of the full model
%   xavg, xstd    participant x fitted parameter: posterior mean and s.d., full model simulated and fitted
%                 (parameter recovery)
%   participants  participant numbers, irep: repetition of the recovery (simulated dataset)

clearvars
addpath('..'); cfg_path = config_paths;

respath = './results/tmp/recovery';
irep    = 1; % repetition of the recovery
nmodel  = 8; % models of run_rec_cond.m

% participant numbers (used in the names of the files)
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');
participants = unique(dat.participant);
ndata = numel(participants);

parnam_sim = {'gamma1','siginf1','sigsel1','epsi'}; % condition 1 parameters (switched on/off) and lapse
par_sim = nan(ndata,numel(parnam_sim));
elbo    = nan(ndata,nmodel,nmodel);
for i = 1:ndata
    fprintf('S%02d\n',participants(i))
    sim = load(fullfile(respath,sprintf('results_sim_cond_S%02d.mat',participants(i))),'out_sim');
    par_sim(i,:) = cellfun(@(p)sim.out_sim{1}.(p),parnam_sim);
    for imodsim = 1:nmodel
        for imodfit = 1:nmodel
            fit = load(fullfile(respath,sprintf('results_rec_cond_S%02d_imodsim%d_imodfit%d_rep%d.mat',participants(i),imodsim,imodfit,irep)),'out_fit');
            elbo(i,imodsim,imodfit) = fit.out_fit.elbo;
            if imodsim == 1 && imodfit == 1 % full model simulated and fitted
                if i == 1
                    xnam = fit.out_fit.xnam;
                    xavg = nan(ndata,numel(xnam));
                    xstd = nan(ndata,numel(xnam));
                end
                xavg(i,:) = fit.out_fit.xavg;
                xstd(i,:) = fit.out_fit.xstd;
            end
        end
    end
end

fname = './results/results_rec_cond.mat';
save(fname,'elbo','par_sim','parnam_sim','xnam','xavg','xstd','participants','irep')
fprintf('saved in %s\n',fname)
