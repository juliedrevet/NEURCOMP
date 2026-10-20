function fit_rec_cond(irec,irep)
% one fit of the model and parameter recovery (run_rec_cond), as one job on a cluster:
% irec = 1:ndata*nmodel*nmodel -> (idata, imodsim, imodfit), participant varying fastest
% (irec = 1:ndata: full model simulated and fitted, i.e. parameter recovery).
% Fits the simulated dataset irep of model imodsim (participant idata, file results_sim_cond_SXX.mat)
% with model imodfit (parameters of condition 1 fixed to 0 if not free), and saves the fit in
% results_rec_cond_SXX_imodsimY_imodfitK_repZ.mat (out_fit, without the variational posterior vp)
%   irec  index of the fit
%   irep  simulated dataset (repetition of the recovery)

addpath('..'); cfg_path = config_paths;
addpath(cfg_path.vbmc)

savepath = fullfile(cfg_path.root,'BEH_analyses','results','tmp','recovery');

% participant numbers (used in the names of the files)
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');
participants = unique(dat.participant);
ndata = numel(participants);

% parameters FREE in each model (columns: gamma, siginf, sigsel), condition 1 only (as in run_rec_cond)
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

[idata,imodsim,imodfit] = ind2sub([ndata nmodel nmodel],irec);
subj = participants(idata);

fprintf('Fitting participant %02d - simulated model %d - fitted model %d - repetition %d ...\n',subj,imodsim,imodfit,irep)

sim = load(fullfile(savepath,sprintf('results_sim_cond_S%02d.mat',subj)),'out_sim');
out_sim = sim.out_sim{imodsim};

cfg_fit = [];
cfg_fit.seqcat  = out_sim.seqcat;
cfg_fit.tltllr  = out_sim.tltllr;
cfg_fit.condtn  = out_sim.condtn;
cfg_fit.seed    = out_sim.seed;
cfg_fit.verbose = false;
cfg_fit.resp    = out_sim.r(:,irep); % simulated dataset of the repetition
if ~isfree(imodfit,1), cfg_fit.gamma1  = 0; end
if ~isfree(imodfit,2), cfg_fit.siginf1 = 0; end
if ~isfree(imodfit,3), cfg_fit.sigsel1 = 0; end

out_fit = fit_model(cfg_fit);
out_fit = rmfield(out_fit,'vp'); % lighter files: without the variational posterior of VBMC
save(fullfile(savepath,sprintf('results_rec_cond_S%02d_imodsim%d_imodfit%d_rep%d.mat',subj,imodsim,imodfit,irep)),'out_fit');
end
