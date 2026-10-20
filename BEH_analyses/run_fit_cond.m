function out_fit = run_fit_cond(idata)
% this scripts launches all fits to enable subsequent model comparison and
% assess which parameters to keep
%
% for each condition we have 8 possible models with following FREE parameters:
%   1. gamma + siginf + sigsel + epsi (fitted once, shared across conditions)
%   2. gamma + siginf + epsi
%   3. gamma + sigsel + epsi
%   4. siginf + sigsel + epsi
%   5. gamma + epsi
%   6. siginf + epsi
%   7. sigsel + epsi
%   8. epsi

addpath('..'); cfg_path = config_paths;
addpath(cfg_path.vbmc)

savepath = './results/tmp';
[~,~] = mkdir(savepath);

% load data table
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = setdiff(unique(dat.participant),[]);
ndata = numel(participants);

fprintf('Processing participant %02d/%02d ...\n',idata,ndata)

ifilt = dat.participant == participants(idata);
cfg.seqcat  = dat.seqcat(ifilt);
cfg.tltllr  = dat.LLR(ifilt);
cfg.resp    = dat.response(ifilt); % response => 1 or 2
cfg.condtn  = dat.condtn(ifilt);
cfg.verbose = true;
cfg.seed    = idata*200;


out_fit = cell(3,8);
fname = sprintf('results_fit_cond_S%02d.mat',participants(idata));

% the parameters of the condition that are not free in the model are fixed to 0
isfree = logical([ ... % columns: gamma, siginf, sigsel
    1 1 1;   % 1. all
    1 1 0;   % 2. gamma + siginf
    1 0 1;   % 3. gamma + sigsel
    0 1 1;   % 4. siginf + sigsel
    1 0 0;   % 5. gamma
    0 1 0;   % 6. siginf
    0 0 1;   % 7. sigsel
    0 0 0]); % 8. epsi
for icond = 1:3
    for imodel = 1:8
        if imodel == 1 && icond > 1 % the full model is the same for all conditions: fitted once
            continue
        end
        cfg_fit = cfg;
        if ~isfree(imodel,1), cfg_fit.(sprintf('gamma%d',icond))  = 0; end
        if ~isfree(imodel,2), cfg_fit.(sprintf('siginf%d',icond)) = 0; end
        if ~isfree(imodel,3), cfg_fit.(sprintf('sigsel%d',icond)) = 0; end
        out_fit{icond,imodel} = fit_model(cfg_fit);
        save(fullfile(savepath,fname),'out_fit');
    end
end
