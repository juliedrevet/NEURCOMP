%% analyses_figure_6
% this script plots:
%   * the main figure predicting compdif (conllr-nullr) from siginf and gamma
%   * the supplementary figure predicting siginf from  conllr and nulllr
%   * the supplementary figure predicting siginf from seqang, absllr, compdif (conllr-nullr)

clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

% only for one condition
icond = 1; % 1:PRE | 2:POST | 3:REF
condnam = {'PRE','POST','REF'};

figpath = './figs';
[~,~] = mkdir(figpath); 

saveornot = false;

if saveornot
    logfile = fullfile(figpath,'analyses_figure_6.txt');
    if exist(logfile,'file'), delete(logfile), end
    diary(logfile)
end

cfg_fig = [];
cfg_fig.figpath = figpath;
cfg_fig.saveornot = saveornot;

rgb = lines;
rgb(3,:) = [.7 .7 .7];

nres = 1e4; % number of bootstrap resample to compute s.e.m.

% best-fitting parameters
load(fullfile(cfg_path.results,'results_fit.mat'),'out_fit')
idx_siginf = contains(out_fit{1}.xnam,'siginf') & contains(out_fit{1}.xnam,num2str(icond));
idx_gamma = contains(out_fit{1}.xnam,'gamma') & contains(out_fit{1}.xnam,num2str(icond));
xavg = cell2mat(cellfun(@(s)getfield(s,'xavg'),out_fit,'UniformOutput',false));
siginf = xavg(:,idx_siginf);
gamma = xavg(:,idx_gamma);

ndata = size(siginf,1);
clearvars out_fit

alphathresh = 0.05;

% neural coding: CONLLR / NULLR / COMPDIF
load('./results/zdec_sensor/zdec_conllr_particle_08Hz.mat','zdec_conllr');
load('./results/zdec_sensor/zdec_nulllr_particle_08Hz.mat','zdec_nulllr');
zdec_conllr = cell2mat(zdec_conllr(:,icond));
zdec_nulllr = cell2mat(zdec_nulllr(:,icond));
zdec_compdif = zdec_conllr-zdec_nulllr;

% neural coding: SEQANG
load(fullfile('./results/zdec_sensor/','zdec_seqang_16Hz.mat'),'zdec_seqang');
zdec_seqang = mean(cell2mat(zdec_seqang(:,icond,:)),3);

% neural coding: ABSLLR
load(fullfile('./results/zdec_sensor/','zdec_absllr_08Hz.mat'),'zdec_absllr','tim_evt');
zdec_absllr = cell2mat(zdec_absllr(:,icond));

xvec = tim_evt;
nsmp = length(xvec);

% compute time window of interest
zdec = zdec_conllr;
yavg = mean(zdec);
ymax = max(yavg);
clust  = bwconncomp(yavg>=0.9*ymax);
ifirst = clust.PixelIdxList{1}(1);
ilast  = clust.PixelIdxList{1}(end);
idxwin_conllr = false(1,length(xvec));
idxwin_conllr(1,ifirst:ilast) = true;
winmin_conllr = xvec(ifirst);
winmax_conllr = xvec(ilast);

fprintf('%s - time window average for CONLLR: %03.1f - %03.1f ms\n',condnam{icond},winmin_conllr*1000,winmax_conllr*1000)

zdec = zdec_seqang;
yavg = mean(zdec);
ymax = max(yavg);
clust  = bwconncomp(yavg>=0.9*ymax);
ifirst = clust.PixelIdxList{1}(1);
ilast  = clust.PixelIdxList{1}(end);
idxwin_seqang = false(1,length(xvec));
idxwin_seqang(1,ifirst:ilast) = true;
winmin_seqang = xvec(ifirst);
winmax_seqang = xvec(ilast);

fprintf('%s - time window average for SEQANG: %03.1f - %03.1f ms\n',condnam{icond},winmin_seqang*1000,winmax_seqang*1000)

zdec = zdec_absllr;
yavg = mean(zdec);
ymax = max(yavg);
clust  = bwconncomp(yavg>=0.9*ymax);
ifirst = clust.PixelIdxList{1}(1);
ilast  = clust.PixelIdxList{1}(end);
idxwin_absllr = false(1,length(xvec));
idxwin_absllr(1,ifirst:ilast) = true;
winmin_absllr = xvec(ifirst);
winmax_absllr = xvec(ilast);

fprintf('%s - time window average for ABSLLR: %03.1f - %03.1f ms\n',condnam{icond},winmin_absllr*1000,winmax_absllr*1000)

%% PEAK LATENCY (jackknife): vertical lines of the median-split plots (SUPP)
xpeak_seqang = mean(jackknife_peak(zdec_seqang,xvec));
xpeak_absllr = mean(jackknife_peak(zdec_absllr,xvec));
xpeak_conllr = mean(jackknife_peak(zdec_conllr,xvec)); % for compression (time window of the true belief update)
fprintf('%s - jackknifed peak latency for SEQANG: %.1f ms\n',condnam{icond},1000*xpeak_seqang)
fprintf('%s - jackknifed peak latency for ABSLLR: %.1f ms\n',condnam{icond},1000*xpeak_absllr)
fprintf('%s - jackknifed peak latency for CONLLR: %.1f ms\n',condnam{icond},1000*xpeak_conllr)

%% COMPDIF: mediansplit on sigma (main)
zdec = zdec_compdif;
varnam = 'compdif';
winmin = winmin_conllr;
winmax = winmax_conllr;
medvar = siginf;
mednam = sprintf('siginf%d',icond);
idxlow = medvar <= median(medvar);
idxhgh = medvar >  median(medvar);
pbar = 1.3; HGT = 4;

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_medsplit_%s',varnam,condnam{icond},mednam);
cfg.xvec = xvec;
cfg.yavg = [mean(zdec(idxlow,:),'omitnan')' mean(zdec(idxhgh,:),'omitnan')'];
cfg.yerr = [std(zdec(idxlow,:),'omitnan')'/sqrt(sum(idxlow)) std(zdec(idxhgh,:),'omitnan')'/sqrt(sum(idxhgh))]; 
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ylabel = 'coding difference (z)';
cfg.pbar = pbar;
cfg.HGT = HGT;
helper_plot_decoding(cfg)

% SUPP: thiner
cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_medsplit_%s_thiner',varnam,condnam{icond},mednam);
cfg.xvec = xvec;
cfg.yavg = [mean(zdec(idxlow,:),'omitnan')' mean(zdec(idxhgh,:),'omitnan')'];
cfg.yerr = [std(zdec(idxlow,:),'omitnan')'/sqrt(sum(idxlow)) std(zdec(idxhgh,:),'omitnan')'/sqrt(sum(idxhgh))];
cfg.win  = [winmin winmax];
cfg.rgb  = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.pbar  = 1;
cfg.HGT   = 4;
cfg.peak = xpeak_conllr;
helper_plot_decoding(cfg)

%% COMPDIF: mediansplit on gamma (main)
zdec = zdec_compdif;
varnam = 'compdif';
winmin = winmin_conllr;
winmax = winmax_conllr;
medvar = gamma;
mednam = sprintf('gamma%d',icond);
idxlow = medvar <= median(medvar);
idxhgh = medvar >  median(medvar);
pbar = 1.3; HGT = 4;

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_medsplit_%s',varnam,condnam{icond},mednam);
cfg.xvec = xvec;
cfg.yavg = [mean(zdec(idxlow,:),'omitnan')' mean(zdec(idxhgh,:),'omitnan')'];
cfg.yerr = [std(zdec(idxlow,:),'omitnan')'/sqrt(sum(idxlow)) std(zdec(idxhgh,:),'omitnan')'/sqrt(sum(idxhgh))];
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ylabel = 'coding difference (z)';
cfg.pbar = pbar;
cfg.HGT = HGT;
helper_plot_decoding(cfg)

%% COMPDIF: time-resolved correlation with siginf (main)
X = log(siginf);
zdec = zdec_compdif;
varnam = 'compdif';
winmin = winmin_conllr;
winmax = winmax_conllr;
yl     = [-1 +2];
yt     = -1:.5:.5;
pbar   = 1.3; 
HGT    = 4;

fprintf('<strong>%s correlation siginf/%s</strong>\n',condnam{icond},varnam)

RHO_time = corr(X,zdec);
fprintf('%s: bootstrap estimate of the s.e.m.\n',varnam)
RHO_time_btstrp = nan(nres,nsmp);
rng(1234,'twister');
for ires = 1:nres % shuffle to estimate the s.e.m.
    if mod(ires,nres/10)==0
        fprintf('...%d/%d\n',ires,nres)
    end
    jres = randsample(ndata,ndata,true); % with replacement
    RHO_time_btstrp(ires,:) = corr(X(jres),zdec(jres,:));
end
% permutation: compute significance cluster
cfg_stat = [];
cfg_stat.alpha_form = alphathresh; % cluster-forming threshold
cfg_stat.nperm = nres;
cfg_stat.xvec = xvec;

stat = cluster_permutation(zdec,X,cfg_stat);

% plot time-resolved correlation
h = false(1,nsmp); 
for ic = 1:numel(stat.clusters)
    if stat.clusters(ic).p < alphathresh 
        h(stat.clusters(ic).idx) = true;
    end
end 
cfg = cfg_fig;
cfg.figname = sprintf('CORR_%s_siginf_%s_pcc',condnam{icond},varnam);
cfg.pbar = pbar;
cfg.HGT = HGT;
cfg.xvec = xvec;
cfg.yavg = RHO_time';
cfg.yerr = std(RHO_time_btstrp)';
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ytick = yt;
cfg.ylim  = yl;
cfg.xt = -0.2:.2:0.8;
cfg.ylabel = 'correlation (r)';
cfg.xlabel = 'time from stimulus onset (s)';
cfg.stats.h = h;
helper_plot_decoding(cfg)

% SUPP: thinner
cfg = cfg_fig;
cfg.figname = sprintf('CORR_%s_siginf_%s_pcc_thinner',condnam{icond},varnam);
cfg.pbar = 1;
cfg.HGT = 4;
cfg.xvec = xvec;
cfg.yavg = RHO_time';
cfg.yerr = std(RHO_time_btstrp)';
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ytick = yt;
cfg.ylim  = yl;
cfg.ylabel = 'correlation (r)';
cfg.xlabel = 'time from stimulus onset (s)';
cfg.stats.h = h;
helper_plot_decoding(cfg)

%% COMPDIF: mean correlation with siginf (SUPP)
zdec = zdec_compdif;
varnam = 'compdif';
idxwin = idxwin_conllr;

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_mean_CORR_siginf',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = log(siginf);
cfg.yavg = mean(zdec(:,idxwin),2);
cfg.rgb = rgb(icond,:);
cfg.xlabel = sprintf('siginf - %s',condnam{icond});cfg.ylabel = 'coding difference (z)';
cfg.ylim = [-.025 .4]; cfg.xlim = log([0.1,2]);
cfg.xtick = log([0.1:0.1:1,2]);
cfg.xticklabel = {'0.1','0.2','','','0.5','','','','','1','2'};
cfg.verbose = false;
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = true;
helper_plot_corr(cfg)

%% COMPDIF: time-resolved correlation with gamma (main)
X = gamma;
zdec = zdec_compdif;
varnam = 'compdif';
winmin = winmin_conllr;
winmax = winmax_conllr;
yl     = [-1 +2];
yt     = -1:.5:.5;
pbar   = 1.3; 
HGT    = 4;

RHO_time = corr(X,zdec);
fprintf('%s: bootstrap estimate of the s.e.m.\n',varnam)
RHO_time_btstrp = nan(nres,nsmp);
rng(1234,'twister');
for ires = 1:nres % shuffle to estimate the s.e.m.
    if mod(ires,nres/10)==0
        fprintf('...%d/%d\n',ires,nres)
    end
    jres = randsample(ndata,ndata,true); % with replacement
    RHO_time_btstrp(ires,:) = corr(X(jres),zdec(jres,:));
end

fprintf('<strong>%s correlation gamma/%s</strong>\n',condnam{icond},varnam)
% permutation: compute significance cluster
cfg_stat = [];
cfg_stat.alpha_form = alphathresh; % cluster-forming threshold
cfg_stat.nperm = nres;
cfg_stat.xvec = xvec;
stat = cluster_permutation(zdec,X,cfg_stat);

% plot time-resolved correlation
h = false(1,nsmp);
for ic = 1:numel(stat.clusters)
    if stat.clusters(ic).p < alphathresh % 
        h(stat.clusters(ic).idx) = true;
    end 
end 
cfg = cfg_fig;
cfg.figname = sprintf('CORR_%s_gamma_%s_pcc',condnam{icond},varnam);
cfg.pbar = pbar;
cfg.HGT = HGT;
cfg.xvec = xvec;
cfg.yavg = RHO_time';
cfg.yerr = std(RHO_time_btstrp)';
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ytick = yt;
cfg.ylim  = yl;
cfg.ylabel = 'correlation (r)';
cfg.xlabel = 'time from stimulus onset (s)';
cfg.stats.h = h;
helper_plot_decoding(cfg)

%% COMPDIF: multiple linear regression: zscore(DECODEDVAR) ~ zscore(GAMMA) + zscore(SIGINF (log)) on mean time window (main)
zdec = zdec_compdif;
X1   = zscore(log(siginf));
X2   = zscore(gamma);
Y    = zscore(mean(zdec(:,idxwin_conllr),2));

mdl = fitlm([X1 X2],Y,'varnames',{'log(siginf)','gamma','compdif'});
disp(mdl);
bavg = mdl.Coefficients.Estimate(2:end); 
berr = mdl.Coefficients.SE(2:end);
pval = mdl.Coefficients.pValue(2:end);
bcov = mdl.CoefficientCovariance(2:end,2:end);

% plot regression coefficients
figure('Color','white','Name',sprintf('compare_gamma_siginf_%s_bars',condnam{icond})); hold on
barh(1:2,bavg,0.5,'EdgeColor',rgb(icond,:),'FaceColor',0.5*(rgb(icond,:)+1),'LineWidth',1);
plot([bavg(1)-berr(1) bavg(1)+berr(1)],[1 1],'k-');
plot([bavg(2)-berr(2) bavg(2)+berr(2)],[2 2],'k-');
xlim([-1,0.5]); ylim([-1 4]);
xlabel('regression coefficient')
ylabel('model parameter')
helper_ILL(gcf,1.5,4.5);
if saveornot
    print(fullfile(figpath,gcf().Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',gcf().Name,figpath)
end
fprintf('STATS against 0:\n')
fprintf('\tlog(siginf): beta = %4.4f, \tp = %4.4f; \tse = %4.4f\n',bavg(1),pval(1),berr(1))
fprintf('\tgamma:       beta = %4.4f, \tp = %4.4f; \tse = %4.4f\n',bavg(2),pval(2),berr(2))

% stats coeff siginf and coeff gamma
beta_dif = bavg(1) - bavg(2);
se_dif = sqrt(bcov(1,1) + bcov(2,2) - 2*bcov(1,2));
t_value = beta_dif / se_dif;
df = mdl.DFE;
p_value = 2 * (1 - tcdf(abs(t_value), df));
fprintf('STATS betwen reg %d and reg %d:\n',1,2)
fprintf('\tt_value = %4.3f, p = %4.4f\n',t_value,p_value)

% partial r2 explained by each variable
R2_full = mdl.Rsquared.Ordinary;
mdl_reduced = fitlm(X2, Y);
R2_reduced = mdl_reduced.Rsquared.Ordinary;
partial_R2_X1 = R2_full - R2_reduced;
mdl_reduced = fitlm(X1, Y);
R2_reduced = mdl_reduced.Rsquared.Ordinary;
partial_R2_X2 = R2_full - R2_reduced;

fprintf('UNIQUE CONTRIBUTION:\n\tsiginf \tpartial R^2 = %.5f \n',partial_R2_X1);
fprintf('\tgamma  \tpartial R^2 = %.5f \n',partial_R2_X2);

% compute partial regression with siginf (no zscore)
zdec = zdec_compdif;
varnam = 'compdif';
X = log(siginf);
Y = mean(zdec(:,idxwin_conllr),2);

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_CORR_siginf_multireg',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = X;
cfg.yavg = Y;
cfg.rgb = rgb(icond,:);
cfg.xlabel = sprintf('siginf - %s',condnam{icond});cfg.ylabel = 'coding difference (z)';
cfg.ylim = [-0.025 .4];
cfg.xlim = log([0.1,2]);
cfg.xtick = log([0.1:0.1:1,2]);
cfg.xticklabel = {'0.1','0.2','','','0.5','','','','','1','2'};
cfg.verbose = false;
cfg.r2 = partial_R2_X1;
cfg.pval = pval(1);
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = false;
helper_plot_corr(cfg)
hold on 

% Prepare grid for prediction
xgrid = linspace(min(cfg.xlim), max(cfg.xlim), 200);
xgrid_z = (xgrid - mean(X)) ./ std(X); % standardize
Xz = [xgrid_z(:) zeros(numel(xgrid_z),1)];

% Predict and get confidence intervals
[Yz_pred, Yz_ci] = predict(mdl, Xz, 'Prediction', 'curve', 'Alpha', 0.05);

% Convert to raw units
Y_pred = Yz_pred * std(Y) + mean(Y);
Y_lo   = Yz_ci(:,1) * std(Y) + mean(Y);
Y_hi   = Yz_ci(:,2) * std(Y) + mean(Y);

plot(xgrid, Y_pred, 'Color', rgb(icond,:), 'LineWidth', 2);
fill([xgrid fliplr(xgrid)], ...
     [Y_lo' fliplr(Y_hi')], ...
     (1+rgb(icond,:))/2, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
if cfg.saveornot
    print(fullfile(figpath,gcf().Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',gcf().Name,figpath)
end

% plot partial regression with gamma (no zscore)
zdec = zdec_compdif;
varnam = 'compdif';
X = gamma;
Y = mean(zdec(:,idxwin_conllr),2);

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_CORR_gamma_multireg',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = X;
cfg.yavg = Y;
cfg.rgb = rgb(icond,:);
cfg.xlabel = sprintf('gamma - %s',condnam{icond});cfg.ylabel = 'coding difference (z)';
cfg.ylim = [-0.025 .4];
cfg.xlim = [-0.6 .2];
cfg.xtick = -0.6:.2:.2;
cfg.verbose = false;
cfg.r2 = partial_R2_X2;
cfg.pval = pval(2);
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = false;
helper_plot_corr(cfg)


%% CONLLR: mediansplit on sigma (SUPP)
zdec = zdec_conllr;
varnam = 'conllr';
winmin = winmin_conllr;
winmax = winmax_conllr;
medvar = siginf;
mednam = sprintf('siginf%d',icond);
idxlow = medvar <= median(medvar);
idxhgh = medvar >  median(medvar);
pbar = 1.3; HGT = 4;

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_medsplit_%s',varnam,condnam{icond},mednam);
cfg.xvec = xvec;
cfg.yavg = [mean(zdec(idxlow,:),'omitnan')' mean(zdec(idxhgh,:),'omitnan')'];
cfg.yerr = [std(zdec(idxlow,:),'omitnan')'/sqrt(sum(idxlow)) std(zdec(idxhgh,:),'omitnan')'/sqrt(sum(idxhgh))];
cfg.win  = [winmin winmax];
cfg.rgb  = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.pbar = pbar;
cfg.HGT  = HGT;
helper_plot_decoding(cfg)

%% CONLLR: time-resolved correlation with siginf (SUPP)
X = log(siginf);
zdec = zdec_conllr;
varnam = 'conllr';
winmin = winmin_conllr;
winmax = winmax_conllr;
yl     = [-1 +2];
yt     = -1:.5:.5;
pbar   = 1.3; 
HGT    = 4;

RHO_time = corr(X,zdec);
fprintf('%s: bootstrap estimate of the s.e.m.\n',varnam)
RHO_time_btstrp = nan(nres,nsmp);
rng(1234,'twister'); 
for ires = 1:nres % shuffle to estimate the s.e.m.
    if mod(ires,nres/10)==0
        fprintf('...%d/%d\n',ires,nres)
    end
    jres = randsample(ndata,ndata,true); % with replacement
    RHO_time_btstrp(ires,:) = corr(X(jres),zdec(jres,:));
end

fprintf('<strong>%s correlation siginf/%s</strong>\n',condnam{icond},varnam)
% permutation: compute significance cluster
cfg_stat = [];
cfg_stat.alpha_form = alphathresh; % cluster-forming threshold
cfg_stat.nperm = nres;
cfg_stat.xvec = xvec;
stat = cluster_permutation(zdec,X,cfg_stat);

% plot time-resolved correlation
h = false(1,nsmp);
for ic = 1:numel(stat.clusters)
    if stat.clusters(ic).p < alphathresh
        h(stat.clusters(ic).idx) = true;
    end
end 
cfg = cfg_fig;
cfg.figname = sprintf('CORR_%s_siginf_%s_pcc',condnam{icond},varnam);
cfg.pbar = pbar;
cfg.HGT = HGT;
cfg.xvec = xvec;
cfg.yavg = RHO_time';
cfg.yerr = std(RHO_time_btstrp)';
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ytick = yt;
cfg.ylim  = yl;
cfg.ylabel = 'correlation (r)';
cfg.xlabel = 'time from stimulus onset (s)';
cfg.stats.h = h;
helper_plot_decoding(cfg)

%% NULLLR: mediansplit on siginf (SUPP)
zdec = zdec_nulllr;
varnam = 'nulllr';
winmin = winmin_conllr;
winmax = winmax_conllr;
medvar = siginf;
mednam = sprintf('siginf%d',icond);
idxlow = medvar <= median(medvar);
idxhgh = medvar >  median(medvar);
pbar = 1.3; HGT = 4;

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_medsplit_%s',varnam,condnam{icond},mednam);
cfg.xvec = xvec;
cfg.yavg = [mean(zdec(idxlow,:),'omitnan')' mean(zdec(idxhgh,:),'omitnan')'];
cfg.yerr = [std(zdec(idxlow,:),'omitnan')'/sqrt(sum(idxlow)) std(zdec(idxhgh,:),'omitnan')'/sqrt(sum(idxhgh))]; 
cfg.win  = [winmin winmax];
cfg.rgb  = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.pbar = pbar;
cfg.HGT  = HGT;
helper_plot_decoding(cfg)

%% NULLLR: time-resolved correlation with siginf (SUPP)
X = log(siginf);
zdec = zdec_nulllr;
varnam = 'nulllr';
winmin = winmin_conllr;
winmax = winmax_conllr;
yl     = [-1 +2];
yt     = -1:.5:.5;
pbar   = 1.3; 
HGT    = 4;

RHO_time = corr(X,zdec);
fprintf('%s: bootstrap estimate of the s.e.m.\n',varnam)
RHO_time_btstrp = nan(nres,nsmp);
rng(1234,'twister');
for ires = 1:nres % shuffle to estimate the s.e.m.
    if mod(ires,nres/10)==0
        fprintf('...%d/%d\n',ires,nres)
    end
    jres = randsample(ndata,ndata,true); % with replacement
    RHO_time_btstrp(ires,:) = corr(X(jres),zdec(jres,:));
end

% permutation: compute significance cluster
cfg_stat = [];
cfg_stat.alpha_form = alphathresh; % cluster-forming threshold
cfg_stat.nperm = nres;
cfg_stat.xvec = xvec;
fprintf('<strong>%s correlation siginf/%s</strong>\n',condnam{icond},varnam)
stat = cluster_permutation(zdec,X,cfg_stat);

% plot time-resolved correlation
h = false(1,nsmp);
for ic = 1:numel(stat.clusters)
    if stat.clusters(ic).p < alphathresh 
        h(stat.clusters(ic).idx) = true; 
    end 
end 
cfg = cfg_fig;
cfg.figname = sprintf('CORR_%s_siginf_%s_pcc',condnam{icond},varnam);
cfg.pbar = pbar;
cfg.HGT = HGT;
cfg.xvec = xvec;
cfg.yavg = RHO_time';
cfg.yerr = std(RHO_time_btstrp)';
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ytick = yt;
cfg.ylim  = yl;
cfg.ylabel = 'correlation (r)';
cfg.xlabel = 'time from stimulus onset (s)';
cfg.stats.h = h;
helper_plot_decoding(cfg)


%% SIGINF: multiple linear regression: zscore(SIGINF (log)) ~ zscore(conllr) + zscore(nulllr) on mean time window (SUPP)
X1 = zscore(mean(zdec_conllr(:,idxwin_conllr),2));
X2 = zscore(mean(zdec_nulllr(:,idxwin_conllr),2));
X  = [X1,X2];
Y  = zscore(log(siginf));
nreg = size(X,2);

mdl = fitlm(X,Y,'varnames',{'conllr','nulllr','log(siginf)'});
bavg = mdl.Coefficients.Estimate(2:end); 
berr = mdl.Coefficients.SE(2:end);
pval = mdl.Coefficients.pValue(2:end);
bcov = mdl.CoefficientCovariance(2:end,2:end);

% plot regression coefficients
figure('Color','white','Name',sprintf('compare_conllr_nulllr_%s_bars',condnam{icond})); hold on
barh(1:2,bavg,0.6,'EdgeColor',rgb(icond,:),'FaceColor',0.5*(rgb(icond,:)+1),'LineWidth',1);
plot([bavg(1)-berr(1) bavg(1)+berr(1)],[1 1],'k-');
plot([bavg(2)-berr(2) bavg(2)+berr(2)],[2 2],'k-');
xlim([-1,0.5]); ylim([-1 4]);
xlabel('regression coefficient')
ylabel('belief update')
helper_ILL(gcf,1.5,4.5);
if saveornot
    print(fullfile(figpath,gcf().Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',gcf().Name,figpath)
end

fprintf('STATS against 0:\n')
for ireg = 1:nreg
    fprintf('\t%s:  beta = %4.4f, \tp = %4.4f; \tse = %4.4f\n',mdl.VariableNames{ireg},bavg(ireg),pval(ireg),berr(ireg))
end

% stats coeff absllr and compdif
beta_dif = bavg(2) - bavg(1);
se_dif = sqrt(bcov(2,2) + bcov(1,1) - 2*bcov(2,1));
t_value = beta_dif / se_dif;
df = mdl.DFE;
p_value = 2 * (1 - tcdf(abs(t_value), df));
fprintf('STATS betwen conllr and nulllr:\n')
fprintf('\tt_value = %4.3f, p = %4.4f\n',t_value,p_value)

% partial r2 explained by each variable
partial_R2 = nan(size(X,2),1);
R2_full = mdl.Rsquared.Ordinary;
for ireg = 1:nreg
    mdl_reduced = fitlm(X(:,setdiff(1:nreg,ireg)), Y);
    R2_reduced = mdl_reduced.Rsquared.Ordinary;
    partial_R2(ireg) = R2_full - R2_reduced;
end

fprintf('UNIQUE CONTRIBUTION to %s:\n',mdl.VariableNames{end})
for ireg = 1:numel(mdl.VariableNames)-1
    fprintf('\t%s:\tpartial R^2 = %.5f \n',mdl.VariableNames{ireg},partial_R2(ireg));
end

% compute partial regression with conllr (no zscore)
zdec = zdec_conllr;
varnam = 'conllr';
X = mean(zdec(:,idxwin_conllr),2);
Y = log(siginf);

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_CORR_siginf_multireg',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = X;
cfg.yavg = Y;
cfg.rgb = rgb(icond,:);
cfg.ylabel = sprintf('siginf - %s',condnam{icond});cfg.xlabel = 'coding precision (z)';
cfg.xlim = [0 .4];
cfg.ylim = log([0.1,2]);
cfg.ytick = log([0.1:0.1:1,2]);
cfg.yticklabel = {'0.1','0.2','','','0.5','','','','','1','2'};
cfg.verbose = false;
cfg.r2 = partial_R2(1);
cfg.pval = pval(1);
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = false;
helper_plot_corr(cfg)
hold on 

% Prepare grid for prediction
xgrid = linspace(min(cfg.xlim), max(cfg.xlim), 200);
xgrid_z = (xgrid - mean(X)) ./ std(X); % standardize
Xz = [xgrid_z(:) zeros(numel(xgrid_z),1)];

% Predict and get confidence intervals
[Yz_pred, Yz_ci] = predict(mdl, Xz, 'Prediction', 'curve', 'Alpha', 0.05);

% Convert to raw units
Y_pred = Yz_pred * std(Y) + mean(Y);
Y_lo   = Yz_ci(:,1) * std(Y) + mean(Y);
Y_hi   = Yz_ci(:,2) * std(Y) + mean(Y);

plot(xgrid, Y_pred, 'Color', rgb(icond,:), 'LineWidth', 2);
fill([xgrid fliplr(xgrid)], ...
     [Y_lo' fliplr(Y_hi')], ...
     (1+rgb(icond,:))/2, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

if cfg.saveornot
    print(fullfile(figpath,gcf().Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',gcf().Name,figpath)
end

% compute partial regression with nulllr (no zscore)
zdec = zdec_nulllr;
varnam = 'nulllr';
X = mean(zdec(:,idxwin_conllr),2);
Y = log(siginf);

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_CORR_siginf_multireg',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = X;
cfg.yavg = Y;
cfg.rgb = rgb(icond,:);
cfg.ylabel = sprintf('siginf - %s',condnam{icond});cfg.xlabel = 'coding precision (z)';
cfg.xlim = [0 .4];
cfg.ylim = log([0.1,2]);
cfg.ytick = log([0.1:0.1:1,2]);
cfg.yticklabel = {'0.1','0.2','','','0.5','','','','','1','2'};
cfg.verbose = false;
cfg.r2 = partial_R2(2);
cfg.pval = pval(2);
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = false;
helper_plot_corr(cfg)

%% SEQANG: mediansplit on siginf (SUPP)
zdec = zdec_seqang;
varnam = 'seqang';
winmin = winmin_seqang;
winmax = winmax_seqang;
yl     = [-0.05 0.8];
medvar = siginf;
mednam = sprintf('siginf%d',icond);
idxlow = medvar <= median(medvar);
idxhgh = medvar >  median(medvar);
pbar = 1; HGT = 4; % thinner for SUPP

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_medsplit_%s',varnam,condnam{icond},mednam);
cfg.xvec = xvec;
cfg.yavg = [mean(zdec(idxlow,:),'omitnan')' mean(zdec(idxhgh,:),'omitnan')'];
cfg.yerr = [std(zdec(idxlow,:),'omitnan')'/sqrt(sum(idxlow)) std(zdec(idxhgh,:),'omitnan')'/sqrt(sum(idxhgh))];
cfg.win  = [winmin winmax];
cfg.rgb  = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ylim  = yl;
cfg.pbar  = pbar;
cfg.HGT   = HGT;
cfg.peak = xpeak_seqang;
helper_plot_decoding(cfg)

%% SEQANG: time-resolved correlation with siginf (SUPP)
X = log(siginf);
zdec = zdec_seqang;
varnam = 'seqang';
winmin = winmin_seqang;
winmax = winmax_seqang;
yl     = [-1 +2];
yt     = -1:.5:.5;
pbar   = 1; 
HGT    = 4;

RHO_time = corr(X,zdec);
fprintf('%s: bootstrap estimate of the s.e.m.\n',varnam)
RHO_time_btstrp = nan(nres,nsmp);
rng(1234,'twister');
for ires = 1:nres % shuffle to estimate the s.e.m.
    if mod(ires,nres/10)==0
        fprintf('...%d/%d\n',ires,nres)
    end
    jres = randsample(ndata,ndata,true); % with replacement
    RHO_time_btstrp(ires,:) = corr(X(jres),zdec(jres,:));
end

fprintf('<strong>%s correlation siginf/%s</strong>\n',condnam{icond},varnam) 
% permutation: compute significance cluster
cfg_stat = [];
cfg_stat.alpha_form = alphathresh; % cluster-forming threshold
cfg_stat.nperm = nres;
cfg_stat.xvec = xvec;
stat = cluster_permutation(zdec,X,cfg_stat);

% plot time-resolved correlation
h = false(1,nsmp); 
for ic = 1:numel(stat.clusters) 
    if stat.clusters(ic).p < alphathresh 
        h(stat.clusters(ic).idx) = true;
    end 
end
cfg = cfg_fig;
cfg.figname = sprintf('CORR_%s_siginf_%s_pcc',condnam{icond},varnam);
cfg.pbar = pbar;
cfg.HGT = HGT;
cfg.xvec = xvec;
cfg.yavg = RHO_time';
cfg.yerr = std(RHO_time_btstrp)';
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ytick = yt;
cfg.ylim  = yl;
cfg.ylabel = 'correlation (r)';
cfg.xlabel = 'time from stimulus onset (s)';
cfg.stats.h = h;
helper_plot_decoding(cfg)

%% SEQANG: mean correlation with siginf (SUPP)
zdec = zdec_seqang;
varnam = 'seqang';
idxwin = idxwin_seqang;

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_mean_CORR_siginf',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = log(siginf);
cfg.yavg = mean(zdec(:,idxwin),2);
cfg.rgb = rgb(icond,:);
cfg.xlabel = sprintf('siginf - %s',condnam{icond});cfg.ylabel = 'coding precision (z)';
cfg.ylim = [0.2 1.2]; cfg.xlim = log([0.1,2]);
cfg.xtick = log([0.1:0.1:1,2]);
cfg.xticklabel = {'0.1','0.2','','','0.5','','','','','1','2'};
cfg.verbose = false;
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = true;
helper_plot_corr(cfg)

%% ABSLLR: mediansplit on siginf (SUPP)
zdec = zdec_absllr;
varnam = 'absllr';
winmin = winmin_absllr;
winmax = winmax_absllr;
medvar = siginf;
mednam = sprintf('siginf%d',icond);
idxlow = medvar <= median(medvar);
idxhgh = medvar >  median(medvar);
pbar = 1; HGT = 4; % thinner for SUPP

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_medsplit_%s',varnam,condnam{icond},mednam);
cfg.xvec = xvec;
cfg.yavg = [mean(zdec(idxlow,:),'omitnan')' mean(zdec(idxhgh,:),'omitnan')'];
cfg.yerr = [std(zdec(idxlow,:),'omitnan')'/sqrt(sum(idxlow)) std(zdec(idxhgh,:),'omitnan')'/sqrt(sum(idxhgh))];
cfg.win  = [winmin winmax];
cfg.rgb  = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.pbar  = pbar;
cfg.HGT   = HGT;
cfg.peak = xpeak_absllr;
helper_plot_decoding(cfg)

%% ABSLLR: time-resolved correlation with siginf (SUPP)
X = log(siginf);
zdec = zdec_absllr;
varnam = 'absllr';
winmin = winmin_absllr;
winmax = winmax_absllr;
yl     = [-1 +2];
yt     = -1:.5:.5;
pbar   = 1; 
HGT    = 4;

RHO_time = corr(X,zdec);
fprintf('%s : bootstrap estimate of the s.e.m.\n',varnam)
RHO_time_btstrp = nan(nres,nsmp);
rng(1234,'twister');
for ires = 1:nres % shuffle to estimate the s.e.m.
    if mod(ires,nres/10)==0
        fprintf('...%d/%d\n',ires,nres)
    end
    jres = randsample(ndata,ndata,true); % with replacement
    RHO_time_btstrp(ires,:) = corr(X(jres),zdec(jres,:));
end

% permutation: compute significance cluster
cfg_stat = [];
cfg_stat.alpha_form = alphathresh; % cluster-forming threshold
cfg_stat.nperm = nres;
cfg_stat.xvec = xvec;
fprintf('<strong>%s correlation siginf/%s</strong>\n',condnam{icond},varnam) 
stat = cluster_permutation(zdec,X,cfg_stat);

% plot time-resolved correlation
h = false(1,nsmp);
for ic = 1:numel(stat.clusters)
    if stat.clusters(ic).p < alphathresh
        h(stat.clusters(ic).idx) = true;
    end 
end 
cfg = cfg_fig;
cfg.figname = sprintf('CORR_%s_siginf_%s_pcc',condnam{icond},varnam);
cfg.pbar = pbar;
cfg.HGT = HGT;
cfg.xvec = xvec;
cfg.yavg = RHO_time';
cfg.yerr = std(RHO_time_btstrp)';
cfg.win = [winmin winmax];
cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
cfg.ytick = yt;
cfg.ylim  = yl;
cfg.ylabel = 'correlation (r)';
cfg.xlabel = 'time from stimulus onset (s)';
cfg.stats.h = h;
helper_plot_decoding(cfg);

%% ABSLLR: mean correlation with siginf (SUPP)
zdec = zdec_absllr;
varnam = 'absllr';
idxwin = idxwin_absllr;

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_mean_CORR_siginf',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = log(siginf);
cfg.yavg = mean(zdec(:,idxwin),2);
cfg.rgb = rgb(icond,:);
cfg.xlabel = sprintf('siginf - %s',condnam{icond});cfg.ylabel = 'coding precision (z)';
cfg.ylim = [0 .4]; cfg.xlim = log([0.1,2]);
cfg.xtick = log([0.1:0.1:1,2]);
cfg.xticklabel = {'0.1','0.2','','','0.5','','','','','1','2'};
cfg.verbose = false;
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = true;
helper_plot_corr(cfg)

%% NEURAL CODINGS: correlation matrix between neural codings (SUPP)
X = zscore([mean(zdec_seqang(:,idxwin_seqang),2) mean(zdec_absllr(:,idxwin_absllr),2) mean(zdec_compdif(:,idxwin_conllr),2)]);
corrX = corr(X);
R_name = {'seqang','absllr','compression'};
figure('Color','white','Name','correlation_matrix_neural_coding');
data = corrX;
nColors = 256;
gamma_cmap = 2.5;

t = linspace(0,1,nColors)'.^gamma_cmap;

r = 1 - t * (1 - 0.7); 
g = 1 - t * (1 - 0.0);
b = 1 - t * (1 - 0.0);
redmap = [r g b];

h = imagesc(corrX,[0 1]);
axis equal;
cmap = colormap(redmap);
colorbar;

[nRows, nCols] = size(corrX);

% Normalize data to [0,1] within the colormap scale
cmin = min(corrX(:));
cmax = max(corrX(:));
normData = (corrX - cmin) / (cmax - cmin);

% Draw numbers
for r = 1:nRows
    for c = 1:nCols
        % Text position (x,y)
        x = c;
        y = r;
        % Decide text color based on brightness of the colormap
        idx = max(1, round(normData(r,c) * (size(colormap,1)-1)) + 1);
        rgb_txt = cmap(idx,:);
        brightness = 0.299*rgb_txt(1) + 0.587*rgb_txt(2) + 0.114*rgb_txt(3);
        
        if brightness < 0.5
            txtColor = [1 1 1];   % background dark -> white text
        else
            txtColor = [0 0 0];   % background light -> black text
        end
        if r==c
            text(x, y, sprintf('%.0f', corrX(r,c)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'Color', txtColor);
        else
            text(x, y, sprintf('%.3f', corrX(r,c)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color', txtColor);
        end
    end
end

set(gca,'XTick',1:nCols,'XTickLabel',R_name,'YTick',1:nRows,'YTickLabel',R_name);
helper_ILL(gcf,1,4);
if saveornot
    print(fullfile(figpath,gcf().Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',gcf().Name,figpath)
end

%% SIGINF: multiple linear regression: zscore(SIGINF (log)) ~ zscore(seqang) + zscore(absllr) + zscore(compdif) on mean time window (SUPP)
X1 = zscore(mean(zdec_seqang(:,idxwin_seqang),2));
X2 = zscore(mean(zdec_absllr(:,idxwin_absllr),2));
X3 = zscore(mean(zdec_compdif(:,idxwin_conllr),2));
X  = [X1,X2,X3];
Y  = zscore(log(siginf));
nreg = size(X,2);

mdl = fitlm(X,Y,'varnames',{'seqang','absllr','compdif','log(siginf)'});
bavg = mdl.Coefficients.Estimate(2:end); 
berr = mdl.Coefficients.SE(2:end);
pval = mdl.Coefficients.pValue(2:end);
bcov = mdl.CoefficientCovariance(2:end,2:end);

figure('Color','white','Name',sprintf('siginf_predicted_from_seqang_absllr_compdif_%s',condnam{icond}));
bar(1:numel(bavg),bavg,0.5,'EdgeColor',rgb(icond,:),'FaceColor',0.5*(rgb(icond,:)+1),'LineWidth',1);
hold on
for ibar = 1:numel(bavg)
    plot([ibar ibar],[bavg(ibar)-berr(ibar) bavg(ibar)+berr(ibar)],'k-');
end
ylim([-.8,0.16]);
xlim([.5 numel(bavg)+0.5]);
ylabel('regression coefficient')
xticklabels({'orientation','evidence strength','compression'})
helper_ILL(gcf,1,4)
if saveornot
    print(fullfile(figpath,gcf().Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',gcf().Name,figpath)
end

fprintf('STATS against 0:\n')
for ireg = 1:nreg
    fprintf('\t%s:  beta = %4.4f, \tp = %4.4f; \tse = %4.4f\n',mdl.VariableNames{ireg},bavg(ireg),pval(ireg),berr(ireg))
end

% stats between coeff absllr and compdif
beta_dif = bavg(3) - bavg(2);
se_dif = sqrt(bcov(3,3) + bcov(2,2) - 2*bcov(3,2));
t_value = beta_dif / se_dif;
df = mdl.DFE;
p_value = 2 * (1 - tcdf(abs(t_value), df));
fprintf('STATS betwen absllr and compdif:\n')
fprintf('\tt_value = %4.3f, p = %4.4f\n',t_value,p_value)

% partial r2 explained by each variable
partial_R2 = nan(size(X,2),1);
R2_full = mdl.Rsquared.Ordinary;
for ireg = 1:nreg
    mdl_reduced = fitlm(X(:,setdiff(1:nreg,ireg)), Y);
    R2_reduced = mdl_reduced.Rsquared.Ordinary;
    partial_R2(ireg) = R2_full - R2_reduced;
end

fprintf('UNIQUE CONTRIBUTION to %s:\n',mdl.VariableNames{end})
for ireg = 1:numel(mdl.VariableNames)-1
    fprintf('\t%s\tpartial R^2 = %.5f \n',mdl.VariableNames{ireg},partial_R2(ireg));
end


%% correlation between participants' compression (average) and fraction optimal
zdec = zdec_compdif;
varnam = 'compdif';
idxwin = idxwin_conllr;
pbar   = 1; 
HGT    = 4;

% load data table organized per sequence
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = unique(dat.participant);
ndata = length(participants);
ncond = 3;

popt = zeros(ndata,1);
popt_std = zeros(ndata,1);
for idata = 1:ndata
    ifilt = dat.participant == participants(idata) & dat.condtn == icond;
    resp = dat.response(ifilt);
    optcat = (3-sign(dat.seqllr(ifilt)))/2;
    popt(idata) = mean(resp == optcat);
    popt_std(idata) = sqrt(popt(idata)*(1-popt(idata))/sum(ifilt));
end

cfg = cfg_fig;
cfg.figname = sprintf('%s_%s_mean_CORR_fraction_optimal',varnam,condnam{icond});
cfg.HGT = 4;
cfg.xavg = mean(zdec(:,idxwin),2);
cfg.yavg = popt;
cfg.rgb = rgb(icond,:);
cfg.ylabel = sprintf('fraction optimal - %s',condnam{icond});cfg.xlabel = 'mean compression (z)';
cfg.ylim = [0.6 1];
cfg.xlim = [-0.05 .3];
cfg.verbose = false;
cfg.pbar = 1;
cfg.identity = false;
cfg.ci = true;
helper_plot_corr(cfg)

if saveornot
    diary off
end
