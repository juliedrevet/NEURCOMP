%% analyses_figure_3_model
% model analyses of Figure 3:
%   * factorized Bayesian model selection over the 8 models of run_fit_cond.m (results_fit_cond.mat)
%   * best-fitting parameters of the winning model (run_fit.m, results_fit.mat): inference imbalance gamma
%     and inference noise siginf in each condition
%   * simulations of the best-fitting model (sim_model.m): sensitivity to evidence and integration kernel

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

HGT = 4.5;

ncond = 3;
ndata = 37;

%% load model fit variants
% lightornot = true: the VBMC posterior (vp) and the copies of the participant's data (cfg) are removed from
% each fit before saving (only the ELBO is used below); the full fits remain in ./results/tmp
lightornot = true;
if ~exist('./results/results_fit_cond.mat','file')
    out_fit = cell(ndata,1);
    parfor idata = 1:ndata
        tmp = run_fit_cond(idata);
        tmp(2:3,1) = tmp(1,1); % full model is shared across conditions
        out_fit{idata} = tmp;
    end
    if lightornot
        for idata = 1:ndata
            for k = 1:numel(out_fit{idata})
                out_fit{idata}{k} = rmfield(out_fit{idata}{k},'vp');
                out_fit{idata}{k}.cfg = rmfield(out_fit{idata}{k}.cfg,{'seqcat','tltllr','resp','condtn'});
            end
        end
    end
    save('./results/results_fit_cond.mat','out_fit')
else
    load('./results/results_fit_cond.mat','out_fit')
end

%% Factorized Bayesian Model Selection
% for each condition, 8 models with following FREE parameters (run_fit_cond):
%   1. gamma + siginf + sigsel + epsi (fitted once, shared across conditions)
%   2. gamma + siginf + epsi
%   3. gamma + sigsel + epsi
%   4. siginf + sigsel + epsi
%   5. gamma + epsi
%   6. siginf + epsi
%   7. sigsel + epsi
%   8. epsi

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

%
elbo = nan(ndata,nmodel,3);
for idata = 1:ndata
    for icond = 1:3
        elbo(idata,:,icond) = cellfun(@(s)s.elbo,out_fit{idata}(icond,:));
    end
end

% random-effects BMS per condition
alphac = nan(ncond,nmodel);
for icond = 1:ncond
    alphac(icond,:) = spm_BMS(elbo(:,:,icond));
end

% family-level: parameter free vs. fixed
% (use ~isfree(:,ipar) & any(isfree,2) to exclude the epsi-only model from the fixed family)
alpha_par = nan(ncond,2,npar); % condition x (free,fixed) x parameter
for ipar = 1:npar
    alpha_par(:,1,ipar) = sum(alphac(:, isfree(:,ipar)),2);
    alpha_par(:,2,ipar) = sum(alphac(:,~isfree(:,ipar)),2);
end

pavg = nan(ncond,2,npar); % condition x (free,fixed) x parameter
pstd = nan(ncond,npar);
pexc = nan(ncond,npar);
for ipar = 1:npar
    for icond = 1:ncond
        fprintf('COND%d - %s FREE?\n',icond,parnam{ipar});
        [pavg(icond,:,ipar),pvar] = drchstat(alpha_par(icond,:,ipar));
        x = drchrnd(alpha_par(icond,:,ipar),1e6);
        pexc(icond,ipar) = mean(x(:,1) > x(:,2));
        pstd(icond,ipar) = sqrt(pvar(1));
        fprintf('\tpavg = %4.4f; pvar = %4.4f; pexc = %4.4f; \n',pavg(icond,1,ipar),pvar(1),pexc(icond,ipar));
    end
end

%% plot family frequencies (POST and PRE only)
ordercond = [2 1 3];
nplot = 2;
wdth = 0.8;
nfam = 2; % free vs. fixed
pbar = .7;

for ipar = 1:npar
    figname = sprintf('model_BMS_%s_or_not',parnam{ipar});
    figure('Color','white','Name',figname);
    hold on
    xlim([0.4,(nfam*nplot)+1.6]);
    xt = [1:nplot (nplot+2):(2*nplot+1)];
    ylim([0 1])
    for i = 1:nplot
        icond = ordercond(i);
        bar(i,pavg(icond,1,ipar),wdth,'EdgeColor',rgb(icond,:),'FaceColor',0.5*(rgb(icond,:)+1),'LineWidth',1);
        plot(i*[1,1],pavg(icond,1,ipar)+pstd(icond,ipar)*[-1,+1],'k-');
        bar(i+nplot+1,pavg(icond,2,ipar),wdth,'EdgeColor',rgb(icond,:),'FaceColor',0.5*(0.5*(rgb(icond,:)+1)+1),'LineWidth',1);
        plot((i+nplot+1)*[1,1],pavg(icond,2,ipar)+pstd(icond,ipar)*[-1,+1],'k-');
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],'TickDir','out');
    set(gca,'XTick',xt,'YTick',0:0.2:1);
    xlabel('model');
    ylabel('model frequency');
    helper_ILL(gcf,pbar,HGT,false);
    if cfg_fig.saveornot
        print(fullfile(cfg_fig.figpath,figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',figname,cfg_fig.figpath)
    end
end


%% Run fit for winning model if not existing
if ~exist('./results/results_fit.mat','file')
    run_fit;
end

%% Best-fitting model parameters
load('./results/results_fit.mat','out_fit')
ndata = size(out_fit,1);
condnam = {'PRE','POST','REF'};
ncond = numel(condnam);

xavg = cell2mat(cellfun(@(x)x.xavg,out_fit,'UniformOutput',false));
xstd = cell2mat(cellfun(@(x)x.xstd,out_fit,'UniformOutput',false));

gamma  = xavg(:,contains(out_fit{1}.xnam,'gamma'));
siginf = xavg(:,contains(out_fit{1}.xnam,'siginf'));
gamma_std  = xstd(:,contains(out_fit{1}.xnam,'gamma'));
siginf_std = xstd(:,contains(out_fit{1}.xnam,'siginf'));
epsi_std   = xstd(:,contains(out_fit{1}.xnam,'epsi'));

fprintf('<strong>Best-fitting inference imbalance gamma:</strong>\n')
for icond = 1:ncond
    fprintf('\t%s: %4.3f±%4.3f\n',condnam{icond},mean(gamma(:,icond)),std(gamma(:,icond))/sqrt(ndata))
end
fprintf('ttest against 0:\n')
for icond = 1:ncond
    X = gamma(:,icond);
    [H,p,ci,stats] = ttest(X);
    effd = stats.tstat/sqrt(size(X,1));
    bf10 = bf.ttest(X);
    if H==0, BFval = 1/bf10; BFstr = 'BF01';
    elseif H == 1, BFval = bf10; BFstr = 'BF10';
    end
    fprintf('\t%s:\tH = %d ; t%02d = %+04.3f ; p = %4.5f ; Cohen''s d = %+04.3f ; CI = [%+5.3f %+5.3f] ; %s = %5.3f\n',condnam{icond},H,stats.df,stats.tstat,p,effd,ci(1),ci(2),BFstr,BFval);
end

X = gamma;
[~,rm] = anova_1way_ranova(X, 'condition'); 
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest between conditions:\n')
for icond = 1:ncond
    jcond = mod(icond,ncond)+1;
    [~,~,~,STATS] = ttest(X(:,icond),X(:,jcond));
    effd = STATS.tstat/sqrt(size(X,1));
    % identify row corresponding post-hoc ttest
    irow = double(c.condition_1) == icond & double(c.condition_2) == jcond;
    P  = c.pValue(irow);
    CI = [c.Lower(irow) c.Upper(irow)];
    H  = P < 0.05;
    bf10 = bf.ttest(X(:,icond),X(:,jcond));
    if H==0, BFval=1/bf10; BFstr='BF01';else,BFval=bf10;BFstr='BF10';end
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+09.5f ; p(corrected) = %5.5f ; Cohen''s d = %+09.5f ; CI(corrected) = [%+5.5f %+5.5f] ; %s = %5.5f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end

fprintf('<strong>Best-fitting inference noise siginf:</strong>\n')
for icond = 1:ncond
    fprintf('\t%s: %4.3f±%4.3f\n',condnam{icond},mean(siginf(:,icond)),std(siginf(:,icond))/sqrt(ndata))
end
X = siginf;
[~,rm] = anova_1way_ranova(X, 'condition'); 
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest between conditions:\n')
for icond = 1:ncond
    jcond = mod(icond,ncond)+1;
    [~,~,~,STATS] = ttest(X(:,icond),X(:,jcond));
    effd = STATS.tstat/sqrt(size(X,1));
    % identify row corresponding post-hoc ttest
    irow = double(c.condition_1) == icond & double(c.condition_2) == jcond;
    P  = c.pValue(irow);
    CI = [c.Lower(irow) c.Upper(irow)];
    H  = P < 0.05;
    bf10 = bf.ttest(X(:,icond),X(:,jcond));
    if H==0, BFval=1/bf10; BFstr='BF01';else,BFval=bf10;BFstr='BF10';end
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+09.5f ; p(corrected) = %5.5f ; Cohen''s d = %+09.5f ; CI(corrected) = [%+5.5f %+5.5f] ; %s = %5.5f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end

% plot inference imbalance
cfg = cfg_fig;
cfg.pbar = 0.5;
cfg.HGT = HGT;
cfg.figname = 'best_fitting_inference_imbalance';
cfg.x = gamma(:,[2 1]);
cfg.rgb = rgb([2 1],:);
cfg.ylim = [-.6 +.2]; 
cfg.ttest = true;
cfg.xoff = .6; 
cfg.width = 0.2;
cfg.mrk_siz = 2;
cfg.ytick = -.6:.2:.2;
cfg.ViolinFaceColor = 'None';
cfg.ViolinEdgeColor = 'None';
cfg.MarkerEdgeColor = 'None';
cfg.xlabel = 'condition';
cfg.ylabel = 'inference imbalance';
cfg.xticklabel = {'POST','PRE'};
cfg.plotsem = true;
cfg.plotmed = false;

helper_plot_violin(cfg)
% ADD REF
plot(xlim,mean(gamma(:,3))*[1 1],'Color',rgb(3,:),'Linewidth',1.5);
% ADD X AXIS
plot(xlim,[0 0],'k:','Linewidth',1);
helper_ILL(gcf,cfg.pbar,cfg.HGT);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end

% plot inference noise
cfg = cfg_fig;
cfg.pbar = 0.5;
cfg.HGT = HGT;
cfg.figname = 'best_fitting_inference_noise';
cfg.x = siginf(:,[2 1]);
cfg.rgb = rgb([2 1],:);
cfg.ylim = [0 2]; 
cfg.ttest = true;
cfg.xoff = .6; 
cfg.width = 0.2;
cfg.mrk_siz = 2;
cfg.ytick = 0:1:2;
cfg.yticklabel = {'0','1.0','2.0'};
cfg.ViolinFaceColor = 'None';
cfg.ViolinEdgeColor = 'None';
cfg.MarkerEdgeColor = 'None';
cfg.xlabel = 'condition';
cfg.ylabel = 'inference noise';
cfg.xticklabel = {'POST','PRE'};
cfg.plotsem = true;
cfg.plotmed = false;

helper_plot_violin(cfg)
% ADD REF
plot(xlim,mean(siginf(:,3))*[1 1],'Color',rgb(3,:),'Linewidth',1.5);
% ADD X AXIS
plot(xlim,[0 0],'k:','Linewidth',1);
helper_ILL(gcf,cfg.pbar,cfg.HGT);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end

%% Best-fitting model simulations
load('./results/results_fit.mat','out_fit')
ndata = size(out_fit,1);
condnam = {'PRE','POST','REF'};
ncond = numel(condnam);

xavg = cell2mat(cellfun(@(x)x.xavg,out_fit,'UniformOutput',false));
xstd = cell2mat(cellfun(@(x)x.xstd,out_fit,'UniformOutput',false));

gamma  = xavg(:,contains(out_fit{1}.xnam,'gamma'));
siginf = xavg(:,contains(out_fit{1}.xnam,'siginf'));
epsi   = xavg(:,contains(out_fit{1}.xnam,'epsi'));

load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = unique(dat.participant);
out_sim = cell(ndata,1);
nsim = 1;

% initialize variables for sigmoid curve sensitivity to evidence
maxbin = max(abs(dat.seqllr));
binedges = [-Inf -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 Inf];
xticks = -4:1:+4;
nbin = numel(binedges)-1;
% perf per condition + STATS
perf_bin_sim = nan(ndata,nbin,ncond);
beta_sim = nan(ndata,ncond,2);
perf_bin_sub = nan(ndata,nbin,ncond);
beta_sub = nan(ndata,ncond,2);

xvec_sim = -4:.01:+4;
y_sim = nan(ndata,numel(xvec_sim),3);
y_sub = nan(ndata,numel(xvec_sim),3);

fun = @(b,x) 1./(1+exp(-x*b));

% initialize variable for integration kernel
B_l_sim = nan(ndata,5,ncond,nsim);
B_r_sim = nan(ndata,5,ncond,nsim);
B_l = nan(ndata,5,ncond); 
B_r = nan(ndata,5,ncond);

for idata = 1:ndata
    fprintf('Processing participant %02d/%02d ...\n',idata,ndata)
    ifilt = dat.participant == participants(idata);
    resp = dat.response(ifilt);
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

    for icond = 1:3
        ifilt = out_sim{idata,1}.condtn == icond ;
        resp_sub = resp(ifilt);
        resp_sim = out_sim{idata,1}.r(ifilt);
        seqcat = out_sim{idata,1}.seqcat(ifilt);
        seqllr = cellfun(@sum,out_sim{idata,1}.tltllr(ifilt));
        beta_sim(idata,icond,:)= glmfit(seqllr, resp_sim==1, 'binomial', 'link', 'logit');

        y_sim(idata,:,icond) = glmval(squeeze(beta_sim(idata,icond,:)), xvec_sim', 'logit');

        binidx = discretize(seqllr,binedges);
        bins = unique(binidx);
        for ibin = 1:numel(bins)
            perf_bin_sim(idata,bins(ibin),icond) = mean(resp_sim(binidx==bins(ibin)) == 1);
        end

        beta_sub(idata,icond,:)= glmfit(seqllr, resp_sub==1, 'binomial', 'link', 'logit');
        y_sub(idata,:,icond) = glmval(squeeze(beta_sub(idata,icond,:)), xvec_sim', 'logit');
        binidx = discretize(seqllr,binedges);
        bins = unique(binidx);
        for ibin = 1:numel(bins)
            perf_bin_sub(idata,bins(ibin),icond) = mean(resp_sub(binidx==bins(ibin)) == 1);
        end

    end

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

% plot sigmoid curve
figure('Color','white','Name','sensitivity_evidence_curve_sim'); hold on
pbar = 1.3;
nbin = size(perf_bin_sim,2);
xvec = -4:1:+4;

Yavg = squeeze(mean(perf_bin_sim));
Yerr = squeeze(std(perf_bin_sim,[],1))/sqrt(ndata);

yavg = squeeze(mean(y_sim));
yerr = squeeze(std(y_sim,[],1)/sqrt(ndata));
yerr_vec = nan(2,numel(xvec_sim),ncond);

for icond = 1:ncond
    yerr_vec(1,:,icond) = reshape(yavg(:,icond),1,numel(xvec_sim)) + reshape(yerr(:,icond),1,numel(xvec_sim));
    yerr_vec(2,:,icond) = reshape(yavg(:,icond),1,numel(xvec_sim)) - reshape(yerr(:,icond),1,numel(xvec_sim));
    patch([xvec_sim,fliplr(xvec_sim)],[yerr_vec(1,:,icond),fliplr(yerr_vec(2,:,icond))], ...
        0.5*(rgb(icond,:)+1),'EdgeColor','none');
end
plot([0 0],ylim,'k')
plot(xlim,[0.5 0.5],'k')
for icond = 1:ncond
    plot(xvec_sim,yavg(:,icond),'Color',rgb(icond,:),'LineWidth',1.5);
end
for icond = 1:ncond
    for i = 1:nbin
        plot(xvec(i)*[1,1],Yavg(i,icond)+Yerr(i,icond)*[+1,-1],'k-');
    end
end
for icond = 1:ncond
    plot(xvec,Yavg(:,icond),'ko','MarkerSize',10,'MarkerFaceColor',rgb(icond,:),'LineWidth',.75);
end

xlim([-4 4])
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'XTick',-4:2:+4,'YTick',0:0.5:1);
xlabel('sequence evidence (logLR)');
ylabel('fraction category A');

fig = helper_ILL(gcf,pbar,HGT);

if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,fig.Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',fig.Name,cfg_fig.figpath)
end

% plot integration kernel
pbar = 0.8;
figure('Color','white','Name','integration_kernel_sim');
for icond = 1:3
    hold on
    % create concatenated vector
    yvec = nan(size(B_l_sim,1),12);
    yvec(:,1:4) = B_l_sim(:,1:4,icond);
    yvec(:,9:12) = B_r_sim(:,2:end,icond);
    xvec = 1:4; % first half
    yavg = mean(yvec(:,xvec));
    xpos = 1:4;
    % s.e.m.
    clearvars yerr
    yerr(1,:) = yavg+std(yvec(:,xvec),[],1)/sqrt(size(yvec,1)); 
    yerr(2,:) = yavg-std(yvec(:,xvec),[],1)/sqrt(size(yvec,1));
    patch([xpos,fliplr(xpos)],[yerr(1,:),fliplr(yerr(2,:))], ...
        0.5*(rgb(icond,:)+1),'EdgeColor','none');
    % avg
    plot(xpos,yavg,'LineWidth',2,'Color',rgb(icond,:));
    
    xvec = 9:12; % second half
    yavg = mean(yvec(:,xvec));
    xpos = 5:8;
    
    % s.e.m.
    yerr(1,:) = yavg+std(yvec(:,xvec),[],1)/sqrt(size(yvec,1)); 
    yerr(2,:) = yavg-std(yvec(:,xvec),[],1)/sqrt(size(yvec,1));
    patch([xpos,fliplr(xpos)],[yerr(1,:),fliplr(yerr(2,:))], ...
        0.5*(rgb(icond,:)+1),'EdgeColor','none');
    % avg
    plot(xpos,yavg,'LineWidth',2,'Color',rgb(icond,:));
end
ylim([0 1.5])
xlim([.5 8.5])
set(gca,'XTick',1:8,'XTickLabel',{'1' '2' '3' '4' '-4' '-3' '-2' '-1'})
xlabel('position in sequence')
ylabel('decision weight')
fig = helper_ILL(gcf,pbar,HGT);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath ,fig.Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n\n',fig.Name,cfg_fig.figpath )
end

anova_2way_ranova(cat(2,B_l_sim(:,1:4,:),B_r_sim(:,2:end,:)), {'seqpos','condition'}); 

% plot sensitivity coefficient bars
cfg = cfg_fig;
cfg.pbar = .5;
cfg.wdth = 0.5;
cfg.HGT  = HGT;
cfg.figname = 'sensitivity_evidence_bars_sim';
cfg.rgb = rgb([2 1],:);
cfg.x = beta_sim(:,[2 1],2);
cfg.plotx = true;
cfg.ttest = true;
cfg.xlabel = 'condition';
cfg.labels = condnam([2 1]);
cfg.ylabel = 'sensitivity coefficient';
cfg.ylim = [0 3.1];cfg.ytick = 0:1:3;cfg.yticklabel={'0','1.0','2.0','3.0'};
helper_plot_overall_bars(cfg);
plot(xlim,mean(beta_sim(:,3,2))*[1,1],'Color',rgb(3,:),'LineWidth',1);
if cfg.saveornot
    print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
end
X = beta_sim(:,:,2);
[stat,rm] = anova_1way_ranova(X, 'condition');
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest between conditions:\n')
for icond = 1:ncond
    jcond = mod(icond,ncond)+1;
    [~,~,~,STATS] = ttest(X(:,icond),X(:,jcond));
    effd = STATS.tstat/sqrt(size(X,1));
    % identify row corresponding post-hoc ttest
    irow = double(c.condition_1) == icond & double(c.condition_2) == jcond;
    P  = c.pValue(irow);
    CI = [c.Lower(irow) c.Upper(irow)];
    H  = P < 0.05;
    bf10 = bf.ttest(X(:,icond),X(:,jcond));
    if H==0, BFval=1/bf10; BFstr='BF01';else,BFval=bf10;BFstr='BF10';end
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+09.5f ; p(corrected) = %5.5f ; Cohen''s d = %+09.5f ; CI(corrected) = [%+5.5f %+5.5f] ; %s = %5.5f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end

figure('Color','white','Name','correlation_human_sim_sensitivity')
hold on
xlim([0 3.1])
ylim([0 3.1])
plot(xlim,ylim,'k:')
for idata = 1:ndata
    for icond = 1:ncond
        plot(beta_sub(idata,icond,2),beta_sim(idata,icond,2),'wo','MarkerFaceColor',rgb(icond,:),'MarkerSize', 4,'LineWidth',0.75);
    end
end
xlabel('human data')
ylabel('model simulations')
% compute correlation across conditions
[rho,pval] = corr(reshape(beta_sub(:,:,2),numel(beta_sub(:,:,2)),1),reshape(beta_sim(:,:,2),numel(beta_sim(:,:,2)),1));
text(.05*(max(xlim)-min(xlim))+min(xlim),.9*(max(ylim)-min(ylim))+min(ylim),sprintf('r^2 = %4.3f\np < %4.3g\n',rho^2,ceil(1000*pval)/1000))
set(gca,'XTick',0:1:3,'YTick',0:1:3,'XTickLabels',{'0','1.0','2.0','3.0'},'YTickLabels',{'0','1.0','2.0','3.0'});
fig = helper_ILL(gcf,1,HGT);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath ,fig.Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n\n',fig.Name,cfg_fig.figpath )
end

fprintf('Sensitivity to evidence: Human participants versus model simulations:\n')
for icond = 1:3
    [rho,pval] = corr(beta_sim(:,icond,2),beta_sub(:,icond,2));
    fprintf('\tcond%d: r2 = %4.3f; p = %4.3f\n',icond,rho^2,pval);
end
