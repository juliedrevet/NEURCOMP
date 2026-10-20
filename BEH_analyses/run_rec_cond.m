%% run_rec_cond
% model and parameter recovery for the factorized Bayesian model selection
% (Supplementary Figure 3a), same models as run_fit_cond:
%   1. SIMULATION (all participants): one parameter set drawn from the prior distributions,
%      used to simulate each of the 8 models below on the sequences of the participant,
%      with nrep simulated datasets per model (one per repetition of the recovery)
%   2. RECOVERY: each simulated dataset fitted back with each of the 8 models
% Only the parameters of condition 1 are switched on/off (the other conditions are
% symmetrical): the parameters of conditions 2 and 3 stay free, as in run_fit_cond.
%
% 8 models with following FREE parameters (in condition 1):
%   1. gamma + siginf + sigsel + epsi
%   2. gamma + siginf + epsi
%   3. gamma + sigsel + epsi
%   4. siginf + sigsel + epsi
%   5. gamma + epsi
%   6. siginf + epsi
%   7. sigsel + epsi
%   8. epsi
%
% Each fit is one job, fit_rec_cond(irec,irep), run in parallel here (or as one job of a job array on a cluster):
% irec = 1:ndata*nmodel*nmodel -> (idata, imodsim, imodfit), participant varying fastest, so that
% irec = 1:ndata are the full model simulated and fitted (parameter recovery).
%
% results saved in ./results/tmp/recovery (analysed in analyses_supp.m, Supplementary Figure 3a):
%   results_sim_cond_SXX.mat                       out_sim: 1 x simulated model (nrep simulated datasets in out_sim{.}.r)
%   results_rec_cond_SXX_imodsimY_imodfitK_repZ.mat out_fit: fit of the dataset irep of model Y with model K

clearvars
clc

addpath('..'); cfg_path = config_paths; % script in BEH_analyses, code and config_paths at the root
addpath(cfg_path.vbmc)

savepath = fullfile(cfg_path.root,'BEH_analyses','results','tmp','recovery'); % full path (also for the parallel workers)
[~,~] = mkdir(savepath);

nrep     = 5;    % simulated datasets per model (repetitions of the recovery)
irep     = 1;    % repetition to recover (simulated dataset)
irec2run = 112:296; % fits to run 

% load data table
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = unique(dat.participant);
ndata = numel(participants);

% parameters FREE in each model (columns: gamma, siginf, sigsel), condition 1 only
isfree = logical([ ...
    1 1 1;   % 1. all
    1 1 0;   % 2. gamma + siginf
    1 0 1;   % 3. gamma + sigsel
    0 1 1;   % 4. siginf + sigsel
    1 0 0;   % 5. gamma
    0 1 0;   % 6. siginf
    0 0 1;   % 7. sigsel
    0 0 0]); % 8. epsi
nmodel = size(isfree,1);

%% 1. SIMULATION (all participants)
for idata = 1:ndata
    fprintf('Simulating participant %02d/%02d ...\n',idata,ndata)
    ifilt = dat.participant == participants(idata);
    cfg = [];
    cfg.seqcat  = dat.seqcat(ifilt);
    cfg.tltllr  = dat.LLR(ifilt);
    cfg.condtn  = dat.condtn(ifilt);
    cfg.seed    = idata*200; % seed of the simulation (and of the fits)
    cfg.nsmp    = nrep;      % simulated datasets

    % simulation parameters, drawn from the prior distributions (same set for all simulated models)
    sd = RandStream('mt19937ar','Seed',idata*100);
    RandStream.setGlobalStream(sd);
    % 1/ inference imbalance
    gamma1 = normrnd(0,0.2); gamma2 = normrnd(0,0.2); gamma3 = normrnd(0,0.2);
    % 2/ inference noise
    siginf1 = exprnd(0.5); siginf2 = exprnd(0.5); siginf3 = exprnd(0.5);
    % 3/ selection noise
    sigsel1 = exprnd(1); sigsel2 = exprnd(1); sigsel3 = exprnd(1);
    % 4/ response lapses ~ beta(1,49)
    epsi = betarnd(1,49);

    % simulate each model (parameters of condition 1 switched off if not free)
    out_sim = cell(1,nmodel);
    for imodsim = 1:nmodel
        cfg_sim = cfg;
        cfg_sim.gamma1  = gamma1*isfree(imodsim,1);  cfg_sim.gamma2  = gamma2;  cfg_sim.gamma3  = gamma3;
        cfg_sim.siginf1 = siginf1*isfree(imodsim,2); cfg_sim.siginf2 = siginf2; cfg_sim.siginf3 = siginf3;
        cfg_sim.sigsel1 = sigsel1*isfree(imodsim,3); cfg_sim.sigsel2 = sigsel2; cfg_sim.sigsel3 = sigsel3;
        cfg_sim.epsi    = epsi;
        out_sim{imodsim} = sim_model(cfg_sim);
    end
    save(fullfile(savepath,sprintf('results_sim_cond_S%02d.mat',participants(idata))),'out_sim');
end

%% 2. RECOVERY: one fit per irec (participant varying fastest, then simulated model, then fitted model)
parfor i = 1:numel(irec2run)
    fit_rec_cond(irec2run(i),irep);
end
fprintf('done.\n')
