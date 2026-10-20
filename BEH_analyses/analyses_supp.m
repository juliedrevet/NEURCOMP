%% analyses_supp
% supplementary figures of the behavioral model
%   * best-fitting parameters of the sensory model versus the category model (Supplementary Figure 3)
%   * model and parameter recovery (Supplementary Figure 4a; needs all the fits of run_rec_cond.m)
%   * best-fitting response lapse and correlation between gamma and siginf (Supplementary Figure 4b-c)
%   * model simulations versus human data: accuracy, fraction optimal, evidence weights after a median
%     split on gamma or siginf (Supplementary Figure 5)


clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

cfg_fig = [];
cfg_fig.saveornot = false;
cfg_fig.figpath = './figs';
[~,~] = mkdir(cfg_fig.figpath);

rgb = lines;
rgb(3,:) = [.7 .7 .7];

condnam = {'PRE','POST','REF'};
ncond = numel(condnam);


%% Supplementary Figure 3: best-fitting parameters: sensory model versus category model
% fit using sensory space (stimulus orientations)
load('./results/results_fit_sensory.mat','out_fit')

xavg_sen = cell2mat(cellfun(@(x)x.xavg,out_fit,'UniformOutput',false));
xstd_sen = cell2mat(cellfun(@(x)x.xstd,out_fit,'UniformOutput',false));

gamma_sen  = xavg_sen(:,contains(out_fit{1}.xnam,'gamma'));
siginf_sen = xavg_sen(:,contains(out_fit{1}.xnam,'siginf'));
epsi_sen   = xavg_sen(:,contains(out_fit{1}.xnam,'epsi'));
gamma_sen_std  = xstd_sen(:,contains(out_fit{1}.xnam,'gamma'));
siginf_sen_std = xstd_sen(:,contains(out_fit{1}.xnam,'siginf'));
epsi_sen_std   = xstd_sen(:,contains(out_fit{1}.xnam,'epsi'));

% fit using categories (logLR)
load('./results/results_fit.mat','out_fit')

xavg_cat = cell2mat(cellfun(@(x)x.xavg,out_fit,'UniformOutput',false));
xstd_cat = cell2mat(cellfun(@(x)x.xstd,out_fit,'UniformOutput',false));

gamma_cat  = xavg_cat(:,contains(out_fit{1}.xnam,'gamma'));
siginf_cat = xavg_cat(:,contains(out_fit{1}.xnam,'siginf'));
epsi_cat   = xavg_cat(:,contains(out_fit{1}.xnam,'epsi'));
gamma_cat_std  = xstd_cat(:,contains(out_fit{1}.xnam,'gamma'));
siginf_cat_std = xstd_cat(:,contains(out_fit{1}.xnam,'siginf'));
epsi_cat_std   = xstd_cat(:,contains(out_fit{1}.xnam,'epsi'));

% correlation GAMMA
fprintf('<strong>Best-fitting inference imbalance gamma - correlation between sensory and category models</strong>\n')
X = gamma_sen;
Y = gamma_cat;
X_std = gamma_sen_std;
Y_std = gamma_cat_std;
acc = 0.2;
maxY = ceil(max([X(:);Y(:)])/acc)*acc;  % same axes for all conditions, sensory and category
minY = floor(min([X(:);Y(:)])/acc)*acc;
xt = minY:acc:maxY;
xtl = cell(1,numel(xt));
for it = 1:numel(xt)
    xtl{it} = sprintf('%3.1f',xt(it));
end
for icond = 1:ncond
    cfg = cfg_fig;
    cfg.figname = sprintf('CORR_sencat_best_fitting_gamma_%s',condnam{icond});
    cfg.xavg = X(:,icond);
    cfg.yavg = Y(:,icond);
    cfg.xstd = X_std(:,icond);
    cfg.ystd = Y_std(:,icond);

    cfg.xlabel = 'sensory gamma';
    cfg.ylabel = 'category gamma';

    cfg.xlim = [minY maxY]; cfg.ylim = cfg.xlim;
    cfg.xtick = xt; cfg.xticklabel = xtl;
    cfg.ytick = cfg.xtick; cfg.yticklabel = cfg.xticklabel;
    cfg.rgb = rgb(icond,:);

    cfg.verbose = false;
    [R,P,RLO,RUP] = corrcoef(cfg.xavg,cfg.yavg);
    fprintf('\t%s: gamma %-4s versus %-4s : RHO = %.5f ; PVAL = %.5f ; RHO^2 = %.5f ; CI = [%+.5f %+.5f]\n',condnam{icond},'sensory','category',R(1,2),P(1,2),R(1,2)^2,RLO(1,2),RUP(1,2));

    helper_plot_corr(cfg)
end

% correlation SIGINF
fprintf('<strong>Best-fitting inference noise siginf - correlation between sensory and category models</strong>\n')
X = siginf_sen;
Y = siginf_cat;
X_std = siginf_sen_std;
Y_std = siginf_cat_std;
acc = 0.5;
maxY = ceil(max([X(:);Y(:)])/acc)*acc;
minY = floor(min([X(:);Y(:)])/acc)*acc;
xt = minY:acc:maxY;
xtl = cell(1,numel(xt));
for it = 1:numel(xt)
    xtl{it} = sprintf('%3.1f',xt(it));
end
for icond = 1:ncond
    cfg = cfg_fig;
    cfg.figname = sprintf('CORR_sencat_best_fitting_siginf_%s',condnam{icond});
    cfg.xavg = X(:,icond);
    cfg.yavg = Y(:,icond);
    cfg.xstd = X_std(:,icond);
    cfg.ystd = Y_std(:,icond);

    cfg.xlabel = 'sensory siginf';
    cfg.ylabel = 'category siginf';

    cfg.xlim = [minY maxY]; cfg.ylim = cfg.xlim;
    cfg.xtick = xt; cfg.xticklabel = xtl;
    cfg.ytick = cfg.xtick; cfg.yticklabel = cfg.xticklabel;
    cfg.rgb = rgb(icond,:);

    cfg.verbose = false;
    [R,P,RLO,RUP] = corrcoef(cfg.xavg,cfg.yavg);
    fprintf('\t%s: siginf %-4s versus %-4s : RHO = %.5f ; PVAL = %.5f ; RHO^2 = %.5f ; CI = [%+.5f %+.5f]\n',condnam{icond},'sensory','category',R(1,2),P(1,2),R(1,2)^2,RLO(1,2),RUP(1,2));

    helper_plot_corr(cfg)
end

% correlation EPSI (one parameter for all conditions)
fprintf('<strong>Best-fitting response lapse epsi - correlation between sensory and category models</strong>\n')
X = epsi_sen;
Y = epsi_cat;
X_std = epsi_sen_std;
Y_std = epsi_cat_std;
acc = 0.02;
maxY = ceil(max([X(:);Y(:)])/acc)*acc;
minY = floor(min([X(:);Y(:)])/acc)*acc;
xt = minY:acc:maxY;
xtl = cell(1,numel(xt));
for it = 1:numel(xt)
    xtl{it} = sprintf('%3.2f',xt(it));
end
cfg = cfg_fig;
cfg.figname = 'CORR_sencat_best_fitting_epsi';
cfg.xavg = X;
cfg.yavg = Y;
cfg.xstd = X_std;
cfg.ystd = Y_std;

cfg.xlabel = 'sensory epsi';
cfg.ylabel = 'category epsi';

cfg.xlim = [minY maxY]; cfg.ylim = cfg.xlim;
cfg.xtick = xt; cfg.xticklabel = xtl;
cfg.ytick = cfg.xtick; cfg.yticklabel = cfg.xticklabel;
cfg.rgb = rgb(4,:);

cfg.verbose = false;
[R,P,RLO,RUP] = corrcoef(cfg.xavg,cfg.yavg);
fprintf('\tepsi %-4s versus %-4s : RHO = %.5f ; PVAL = %.5f ; RHO^2 = %.5f ; CI = [%+.5f %+.5f]\n','sensory','category',R(1,2),P(1,2),R(1,2)^2,RLO(1,2),RUP(1,2));

helper_plot_corr(cfg)

%% Supplementary Figure 4a: model and parameter recovery
HGT_rec = 2.5; % height of the recovery figures

% 8 models with following FREE parameters in condition 1 (same as run_fit_cond):
modnam = {'all','gamma+siginf','gamma+sigsel','siginf+sigsel','gamma','siginf','sigsel','epsi'};
nmodel = numel(modnam);

parnam = {'gamma','siginf','sigsel'};
npar = numel(parnam);
isfree = logical([ ...
    1 1 1;   % 1. all
    1 1 0;   % 2. gamma + siginf
    1 0 1;   % 3. gamma + sigsel
    0 1 1;   % 4. siginf + sigsel
    1 0 0;   % 5. gamma
    0 1 0;   % 6. siginf
    0 0 1;   % 7. sigsel
    0 0 0]); % 8. epsi

% load recovery
% elbo: participant x simulated model x fitted model; par_sim: parameters of the simulated full model;
% xavg, xstd: posterior mean and s.d. of the full model simulated and fitted (parameters xnam)
load('./results/results_rec_cond.mat','elbo','par_sim','parnam_sim','xnam','xavg','xstd','irep')
xnam_fit = xnam; 
xnam = {'gamma1' 'siginf1' 'sigsel1' 'epsi'};
nrep = 1;

% recovery: parameter recovery, model with all parameters, simulated and fitted
RHO = nan(numel(xnam),nrep); % correlation between simulated and recovered parameters
for ipar = 1:numel(xnam)
    R = corrcoef(par_sim(:,strcmp(parnam_sim,xnam{ipar})),xavg(:,strcmp(xnam_fit,xnam{ipar})));
    RHO(ipar,1) = R(1,2);
end
fprintf('<strong>Parameter recovery (full model): correlation between simulated and recovered parameters</strong>\n')
fprintf('mean [min max] over %d repetition(s):\n',nrep)
for ipar = 1:numel(xnam)
    fprintf('\t%-10s : RHO = %.5f [%.5f %.5f]\n',xnam{ipar},mean(RHO(ipar,:)),min(RHO(ipar,:)),max(RHO(ipar,:)));
end
fprintf('repetition %d:\n',irep)
for ipar = 1:numel(xnam)
    cfg = cfg_fig;
    cfg.HGT = HGT_rec;
    cfg.figname = sprintf('parameter_recovery_%s',xnam{ipar});
    if strcmp(xnam{ipar},'epsi')
        cfg.rgb = [.5 .5 .5];
    else
        cfg.rgb = rgb(1,:);
    end
    cfg.xavg = par_sim(:,strcmp(parnam_sim,xnam{ipar}));
    cfg.yavg = xavg(:,strcmp(xnam_fit,xnam{ipar}));
    cfg.ystd = xstd(:,strcmp(xnam_fit,xnam{ipar}));

    if contains(xnam{ipar},'gamma')
        acc = .6;
    elseif contains(xnam{ipar},'siginf')
        acc = .5;
    elseif contains(xnam{ipar},'sigsel')
        acc = 1;
    elseif contains(xnam{ipar},'epsi')
        acc = .02;
    end
    maxXY = ceil(max([cfg.xavg;cfg.yavg])/acc)*acc;
    minXY = floor(min([cfg.xavg;cfg.yavg])/acc)*acc;
    xt = minXY:acc:maxXY;

    cfg.ylim = [minXY maxXY];cfg.xlim = cfg.ylim;
    cfg.xtick = xt; cfg.ytick = xt;
    cfg.xlabel = sprintf('%s simulated',xnam{ipar});cfg.ylabel = sprintf('%s recovered',xnam{ipar});
    cfg.verbose = false;
    [R,P,RLO,RUP] = corrcoef(cfg.xavg,cfg.yavg);
    fprintf('\t%-10s : RHO = %.5f ; PVAL = %.5f ; RHO^2 = %.5f ; CI = [%+.5f %+.5f]\n',xnam{ipar},R(1,2),P(1,2),R(1,2)^2,RLO(1,2),RUP(1,2));
    helper_plot_corr(cfg)
end

% recovery: model recovery for each repetition
rng(1234,'twister'); % exceedance probabilities (Monte Carlo) reproducible
pexc     = nan(nmodel,nmodel,nrep);  % simulated model x fitted model x repetition
pexc_fam = nan(2,2,npar,nrep);       % simulated family x recovered family x parameter x repetition
for j = 1:nrep
    % random-effects BMS for each simulated model
    alphac = nan(nmodel,nmodel); % simulated model x fitted model
    for imodsim = 1:nmodel
        alphac(imodsim,:) = spm_BMS(squeeze(elbo(:,imodsim,:,j)));
    end

    % model recovery, all models: exceedance probability of each simulated model, for each recovered model
    for imodfit = 1:nmodel
        x = drchrnd(alphac(:,imodfit),1e6);
        [~,J] = max(x,[],2);
        for imodsim = 1:nmodel
            pexc(imodsim,imodfit,j) = mean(J == imodsim);
        end
    end

    % model recovery, factorized: parameter free vs. fixed (KO)
    for ipar = 1:npar
        on = isfree(:,ipar);
        % pooled model probabilities: simulated family (free,KO) x recovered family (free,KO)
        alpha_fam = [sum(alphac(on,on),'all')  sum(alphac(on,~on),'all'); ...
                     sum(alphac(~on,on),'all') sum(alphac(~on,~on),'all')];
        for ifit = 1:2
            x = drchrnd(alpha_fam(:,ifit),1e6);
            [~,J] = max(x,[],2);
            for isim = 1:2
                pexc_fam(isim,ifit,ipar,j) = mean(J == isim);
            end
        end
    end
end

% recovery: model recovery, all models (average over repetitions)
pexc_avg = mean(pexc,3);
fprintf('<strong>Model recovery (exceedance probability of the simulated model, diagonal; mean [min max] over %d repetition(s)):</strong>\n',nrep)
for imod = 1:nmodel
    fprintf('\t%-14s: %5.3f [%5.3f %5.3f]\n',modnam{imod},floor(1000*pexc_avg(imod,imod))/1000, ...
        floor(1000*min(pexc(imod,imod,:)))/1000,floor(1000*max(pexc(imod,imod,:)))/1000);
end

pbar = 1;
figname = 'model_recovery_all_models';
figure('Color','white','Name',figname);
imagesc(pexc_avg',[0 1]) % rows: recovered models, columns: simulated models
colormap(flipud(gray)); colorbar('off');
for imod = 1:nmodel
    text(imod,imod,sprintf('%5.3f',floor(1000*pexc_avg(imod,imod))/1000),'Color','w','HorizontalAlignment','center')
end
set(gca,'TickDir','out','XTick',1:nmodel,'YTick',1:nmodel,'XTickLabel',modnam,'YTickLabel',modnam)
xlabel('simulated models');
ylabel('recovered models');
helper_ILL(gcf,pbar,HGT_rec);
set(gca,'Layer','top','Box','on')
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,cfg_fig.figpath)
end

% recovery: model recovery, factorized, parameter free vs. fixed (KO) (average over repetitions)
fprintf('<strong>Factorized model recovery (exceedance probability of the simulated family, diagonal; mean [min max] over %d repetition(s)):</strong>\n',nrep)
for ipar = 1:npar
    pexc_fam_avg = mean(pexc_fam(:,:,ipar,:),4);
    fprintf('\t%-7s: free = %5.3f [%5.3f %5.3f] ; KO = %5.3f [%5.3f %5.3f]\n',parnam{ipar}, ...
        floor(1000*pexc_fam_avg(1,1))/1000,floor(1000*min(pexc_fam(1,1,ipar,:)))/1000,floor(1000*max(pexc_fam(1,1,ipar,:)))/1000, ...
        floor(1000*pexc_fam_avg(2,2))/1000,floor(1000*min(pexc_fam(2,2,ipar,:)))/1000,floor(1000*max(pexc_fam(2,2,ipar,:)))/1000);

    figname = sprintf('model_recovery_%s_or_not',parnam{ipar});
    figure('Color','white','Name',figname);
    imagesc(pexc_fam_avg',[0 1]) % rows: recovered family, columns: simulated family
    colormap(flipud(gray)); colorbar('off');
    for i = 1:2
        text(i,i,sprintf('%5.3f',floor(1000*pexc_fam_avg(i,i))/1000),'Color','w','HorizontalAlignment','center')
    end
    set(gca,'TickDir','out','XTick',1:2,'YTick',1:2,'XTickLabel',{'free','KO'},'YTickLabel',{'free','KO'})
    xlabel('simulated models');
    ylabel('recovered models');
    title(sprintf('%s?',parnam{ipar}))
    helper_ILL(gcf,pbar,HGT_rec);
    set(gca,'Layer','top','Box','on')
    if cfg_fig.saveornot
        print(fullfile(cfg_fig.figpath,figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',figname,cfg_fig.figpath)
    end
end

%% Supplementary Figure 4b-c: best-fitting response lapse, correlation between gamma and siginf
load('./results/results_fit.mat','out_fit')
ndata = size(out_fit,1);

xavg_cat = cell2mat(cellfun(@(x)x.xavg,out_fit,'UniformOutput',false));
xstd_cat = cell2mat(cellfun(@(x)x.xstd,out_fit,'UniformOutput',false));

gamma  = xavg_cat(:,contains(out_fit{1}.xnam,'gamma'));
siginf = xavg_cat(:,contains(out_fit{1}.xnam,'siginf'));
epsi   = xavg_cat(:,contains(out_fit{1}.xnam,'epsi'));
gamma_std  = xstd_cat(:,contains(out_fit{1}.xnam,'gamma'));
siginf_std = xstd_cat(:,contains(out_fit{1}.xnam,'siginf'));
epsi_std   = xstd_cat(:,contains(out_fit{1}.xnam,'epsi'));

% b. best-fitting response lapse (one parameter for all conditions)
fprintf('<strong>Best-fitting response lapse epsi:</strong> %5.4f±%5.4f\n',mean(epsi),std(epsi)/sqrt(ndata))
cfg = cfg_fig;
cfg.figname = 'best_fitting_epsi';
cfg.x = epsi;
cfg.rgb = [0.5 0.5 0.5];
cfg.plotx = true;
cfg.xticklabel = {'all conditions'};
cfg.xlabel = 'condition';
cfg.ylabel = 'response lapse';
cfg.ylim = [0 0.1];cfg.ytick = 0:.05:.1;
helper_plot_overall_bars(cfg);

% c. correlation between the best-fitting parameters inside each condition: inference imbalance gamma
% and inference noise siginf
fprintf('<strong>Correlation between best-fitting gamma and siginf in each condition:</strong>\n')
for icond = 1:ncond
    X = gamma(:,icond);
    X_std = gamma_std(:,icond);
    Y = siginf(:,icond);
    Y_std = siginf_std(:,icond);
    accX = 0.2; accY = 0.5;
    maxX = ceil(max(X(:))/accX)*accX;
    maxY = ceil(max(Y(:))/accY)*accY;
    minX = floor(min(X(:))/accX)*accX;
    minY = floor(min(Y(:))/accY)*accY;
    xt = minX:accX:maxX;
    yt = minY:accY:maxY;
    xtl = cell(1,numel(xt));
    for it = 1:numel(xt)
        xtl{it} = sprintf('%3.1f',xt(it));
    end
    ytl = cell(1,numel(yt));
    for it = 1:numel(yt)
        ytl{it} = sprintf('%3.1f',yt(it));
    end
    cfg = cfg_fig;
    cfg.figname = sprintf('CORR_gamma_siginf_%s',condnam{icond});
    cfg.rgb = rgb(icond,:);
    cfg.xavg = X;
    cfg.yavg = Y;
    cfg.xstd = X_std;
    cfg.ystd = Y_std;
    cfg.xlabel = sprintf('%s gamma',condnam{icond});
    cfg.ylabel = sprintf('%s siginf',condnam{icond});
    cfg.xlim = [minX maxX]; cfg.ylim = [minY maxY];
    cfg.xtick = xt; cfg.ytick = yt;
    cfg.xticklabel = xtl; cfg.yticklabel = ytl;
    cfg.verbose = true; cfg.identity = false;
    helper_plot_corr(cfg);
end

%% Supplementary Figure 5: model simulations versus human data
HGT = 4.5;
% data
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');
participants = unique(dat.participant);

% simulation of best-fitting model
nsim = 1;
B_l_sim = nan(ndata,5,ncond,nsim);
B_r_sim = nan(ndata,5,ncond,nsim);
out_sim = cell(ndata,1);
for idata = 1:ndata
    fprintf('Processing participant %02d/%02d ...\n',idata,ndata)
    ifilt = dat.participant == participants(idata);
    cfg_sim = [];
    cfg_sim.seqcat = dat.seqcat(ifilt);
    cfg_sim.seqtlt = dat.seqtlt(ifilt);
    cfg_sim.seqlen = dat.seqlen(ifilt);
    cfg_sim.tltllr = dat.LLR(ifilt);
    cfg_sim.condtn = dat.condtn(ifilt);
    cfg_sim.nsmp = nsim;
    cfg_sim.seed = 101*idata;
    
    cfg_sim.gamma1  = gamma(idata,1);
    cfg_sim.gamma2  = gamma(idata,2);
    cfg_sim.gamma3  = gamma(idata,3);
    cfg_sim.siginf1 = siginf(idata,1);
    cfg_sim.siginf2 = siginf(idata,2);
    cfg_sim.siginf3 = siginf(idata,3);
    cfg_sim.epsi    = epsi(idata);
    cfg_sim.sigsel1 = 0;
    cfg_sim.sigsel2 = 0;
    cfg_sim.sigsel3 = 0;

    out_sim{idata,1} = sim_model(cfg_sim);

    nseq = size(out_sim{idata}.tltllr,1);
    X_l = zeros(nseq,12);
    X_r = zeros(nseq,12);

    for iseq = 1:nseq
        X_l(iseq,1:out_sim{idata}.seqlen(iseq)) = out_sim{idata}.tltllr{iseq};
        X_r(iseq,(12-out_sim{idata}.seqlen(iseq)+1):12) = out_sim{idata}.tltllr{iseq};
    end
    for icond = 1:3
        ifilt = out_sim{idata}.condtn == icond;
        Y = out_sim{idata}.r(ifilt,:) == 1;
        % right aligned
        X = X_r(ifilt,:);
        X = [sum(X(:,1:8),2) X(:,9:end)];
        for isim = 1:nsim
            B_tmp  = glmfit(X,Y(:,isim),'binomial','link','probit');
            B_r_sim(idata,:,icond,isim) = B_tmp(2:end)';
        end
        % left aligned
        X = X_l(ifilt,:);
        X = [X(:,1:4) sum(X(:,5:end),2)];
        for isim = 1:nsim
            B_tmp  = glmfit(X,Y(:,isim),'binomial','link','probit');
            B_l_sim(idata,:,icond,isim) = B_tmp(2:end)';
        end
    end
end

% a. accuracy and fraction optimal: human data versus simulations of the best-fitting model
perf = zeros(ndata,ncond);
perf_opt = zeros(ndata,ncond);
for icond = 1:ncond
    for idata = 1:ndata
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        resp = dat.response(ifilt);
        seqcat = dat.seqcat(ifilt);
        seqllr = dat.seqllr(ifilt);
        perf(idata,icond) = mean(resp == seqcat);
        perf_opt(idata,icond) = mean(sign(1.5-resp) == sign(seqllr));
    end
end

out_sim_perf = cell(ndata,1);
for idata = 1:ndata
    fprintf('Processing participant %02d/%02d ...\n',idata,ndata)
    ifilt = dat.participant == participants(idata);
    cfg_sim = [];
    cfg_sim.seqcat = dat.seqcat(ifilt);
    cfg_sim.seqtlt = dat.seqtlt(ifilt);
    cfg_sim.seqlen = dat.seqlen(ifilt);
    cfg_sim.tltllr = dat.LLR(ifilt);
    cfg_sim.condtn = dat.condtn(ifilt);
    cfg_sim.nsmp = 1000;
    cfg_sim.seed = 100*idata;

    cfg_sim.gamma1  = gamma(idata,1);
    cfg_sim.gamma2  = gamma(idata,2);
    cfg_sim.gamma3  = gamma(idata,3);
    cfg_sim.siginf1 = siginf(idata,1);
    cfg_sim.siginf2 = siginf(idata,2);
    cfg_sim.siginf3 = siginf(idata,3);
    cfg_sim.epsi    = epsi(idata);
    cfg_sim.sigsel1 = 0;
    cfg_sim.sigsel2 = 0;
    cfg_sim.sigsel3 = 0;

    out_sim_perf{idata,1} = sim_model(cfg_sim);
end

perf_sim = cell2mat(cellfun(@(s)s.pcor_avg,out_sim_perf,'UniformOutput',false));
perf_sim_opt = cell2mat(cellfun(@(s)s.pcor_opt,out_sim_perf,'UniformOutput',false));

fprintf('<strong>Accuracy: human data versus simulations of the best-fitting model</strong>\n')
for icond = 1:ncond
    cfg = cfg_fig;
    cfg.figname = sprintf('perf_sub_sim_%s',condnam{icond});
    cfg.xavg = perf(:,icond);
    cfg.yavg = perf_sim(:,icond);
    cfg.rgb = rgb(icond,:);
    cfg.xlim = [0.5,1];
    cfg.ylim = [0.5,1];
    cfg.xlabel = sprintf('%s accuracy',condnam{icond});
    cfg.ylabel = sprintf('%s simulated accuracy',condnam{icond});
    helper_plot_corr(cfg)
end
fprintf('<strong>Fraction optimal: human data versus simulations of the best-fitting model</strong>\n')
for icond = 1:ncond
    cfg = cfg_fig;
    cfg.figname = sprintf('perf_opt_sub_sim_%s',condnam{icond});
    cfg.xavg = perf_opt(:,icond);
    cfg.yavg = perf_sim_opt(:,icond);
    cfg.rgb = rgb(icond,:);
    cfg.xlim = [0.5,1];
    cfg.ylim = [0.5,1];
    cfg.xlabel = sprintf('%s fraction optimal',condnam{icond});
    cfg.ylabel = sprintf('%s simulated fraction optimal',condnam{icond});
    helper_plot_corr(cfg)
end

% b. evidence weights (first and last 4 positions) after a median split on gamma or on siginf, for the
% human data and for the model simulations
X_l = zeros(size(dat,1),12);
X_r = zeros(size(dat,1),12);
for iseq = 1:size(dat,1)
    X_l(iseq,1:dat.seqlen(iseq)) = dat.LLR{iseq};
    X_r(iseq,(12-dat.seqlen(iseq)+1):12) = dat.LLR{iseq};
end
B_l_sub = nan(ndata,5,ncond);
B_r_sub = nan(ndata,5,ncond);
for idata = 1:ndata
    for icond = 1:ncond
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        Y = dat.response(ifilt,:)==1;
        % right aligned
        X = X_r(ifilt,:);
        X = [sum(X(:,1:8),2) X(:,9:end)];
        B_tmp  = glmfit(X,Y,'binomial','link','probit');
        B_r_sub(idata,:,icond) = B_tmp(2:end)';
        % left aligned
        X = X_l(ifilt,:);
        X = [X(:,1:4) sum(X(:,5:end),2)];
        B_tmp  = glmfit(X,Y,'binomial','link','probit');
        B_l_sub(idata,:,icond) = B_tmp(2:end)';
    end
end

% human evidence weights without split
pbar = 0.8;
figure('Color','white','Name','glmfit_seqpos_both_bis');
for icond = 1:3
    hold on
    % create concatenated vector
    yvec = nan(size(B_l_sub,1),12);
    yvec(:,1:4) = B_l_sub(:,1:4,icond);
    yvec(:,9:12) = B_r_sub(:,2:end,icond);
    for ihalf = 1:2
        if ihalf == 1 % first half
            xvec = 1:4; xpos = 1:4;
        else % second half
            xvec = 9:12; xpos = 5:8;
        end
        yavg = mean(yvec(:,xvec));
        % s.e.m.
        clearvars yerr
        yerr(1,:) = yavg+std(yvec(:,xvec),[],1)/sqrt(size(yvec,1));
        yerr(2,:) = yavg-std(yvec(:,xvec),[],1)/sqrt(size(yvec,1));
        patch([xpos,fliplr(xpos)],[yerr(1,:),fliplr(yerr(2,:))], ...
            0.5*(rgb(icond,:)+1),'EdgeColor','none');
        % avg
        plot(xpos,yavg,'LineWidth',2,'Color',rgb(icond,:));
    end
end
ylim([0 1.5])
xlim([.5 8.5])
set(gca,'XTick',1:8,'XTickLabel',{'1' '2' '3' '4' '-4' '-3' '-2' '-1'})
xlabel('sequence position')
ylabel('sensitivity to evidence')
fig = helper_ILL(gcf,pbar,HGT);
set(fig,'PaperPositionMode','manual', ...
    'PaperPosition',[2.5,13,16,HGT],'PaperUnits','centimeters', ...
    'PaperType','A4','PaperOrientation','portrait');
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,fig.Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n\n',fig.Name,cfg_fig.figpath)
end

parsplit = {'gamma','siginf'};
datsplit = {'','_model_simulations'}; % human data, model simulations
pbar = 0.8;
for ipar = 1:numel(parsplit)
    if ipar == 1
        par = gamma;
    else
        par = siginf;
    end
    for isplit = 1:numel(datsplit)
        if isplit == 1
            B_l_split = B_l_sub;
            B_r_split = B_r_sub;
        else
            B_l_split = mean(B_l_sim,4,'omitnan');
            B_r_split = mean(B_r_sim,4,'omitnan');
        end
        figure('Color','white','Name',sprintf('glmfit_seqpos_mediansplit_%s%s',parsplit{ipar},datsplit{isplit}));
        for icond = 1:3
            par_low = par(:,icond)<median(par(:,icond));
            hold on
            B_l_low = B_l_split(par_low,:,icond);
            B_r_low = B_r_split(par_low,:,icond);
            B_l_hgh = B_l_split(~par_low,:,icond);
            B_r_hgh = B_r_split(~par_low,:,icond);

            % create concatenated vector
            yvec_low            = nan(size(B_l_low,1),12);
            yvec_low(:,1:4)     = B_l_low(:,1:4);
            yvec_low(:,9:12)    = B_r_low(:,2:end);
            yvec_hgh            = nan(size(B_l_hgh,1),12);
            yvec_hgh(:,1:4)     = B_l_hgh(:,1:4);
            yvec_hgh(:,9:12)    = B_r_hgh(:,2:end);

            for ihalf = 1:2
                if ihalf == 1 % first half
                    xvec = 1:4; xpos = 1:4;
                else % second half
                    xvec = 9:12; xpos = 5:8;
                end
                yavg_low = mean(yvec_low(:,xvec));
                yavg_hgh = mean(yvec_hgh(:,xvec));
                % s.e.m. low
                clearvars yerr
                yerr(1,:) = yavg_low+std(yvec_low(:,xvec),[],1)/sqrt(size(yvec_low,1));
                yerr(2,:) = yavg_low-std(yvec_low(:,xvec),[],1)/sqrt(size(yvec_low,1));
                patch([xpos,fliplr(xpos)],[yerr(1,:),fliplr(yerr(2,:))], ...
                    0.5*(rgb(icond,:)+1),'EdgeColor','none');
                % s.e.m. high
                yerr(1,:) = yavg_hgh+std(yvec_hgh(:,xvec),[],1)/sqrt(size(yvec_hgh,1));
                yerr(2,:) = yavg_hgh-std(yvec_hgh(:,xvec),[],1)/sqrt(size(yvec_hgh,1));
                patch([xpos,fliplr(xpos)],[yerr(1,:),fliplr(yerr(2,:))], ...
                    0.5*(rgb(icond,:)+1),'EdgeColor','none');
                % avg
                plot(xpos,yavg_low,'LineWidth',2,'Color',rgb(icond,:));
                plot(xpos,yavg_hgh,'LineWidth',2,'LineStyle',':','Color',rgb(icond,:));
            end
        end
        ylim([0 1.5])
        xlim([.5 8.5])
        set(gca,'XTick',1:8,'XTickLabel',{'1' '2' '3' '4' '-4' '-3' '-2' '-1'})
        xlabel('sequence position')
        ylabel('sensitivity to evidence')
        fig = helper_ILL(gcf,pbar,HGT);
        set(fig,'PaperPositionMode','manual', ...
            'PaperPosition',[2.5,13,16,HGT],'PaperUnits','centimeters', ...
            'PaperType','A4','PaperOrientation','portrait');
        if cfg_fig.saveornot
            print(fullfile(cfg_fig.figpath,fig.Name),'-vector','-dpdf');
            fprintf('\nfigure %s saved in %s !\n\n',fig.Name,cfg_fig.figpath)
        end
    end
end
