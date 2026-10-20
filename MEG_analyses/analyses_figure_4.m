%% analyses_figure_4
% this script runs the analyses to reproduce the main Figure 4
%   seqang: time course
%   seqang: source for each condition
%   seqang: OPT/SUB: no time window problem for seqang
%   absllr: time course, decoding against 0 in each condition (cluster-based permutation test)
%   absllr: source for each condition
%   absllr: OPT/SUB: time window condition-specific
%   supp: seqang (POST, PRE) and absllr (PRE) OPT/SUB difference at the source level

clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

figpath = './figs';
[~,~] = mkdir(figpath);
saveornot = false;
plotsource = false;

if saveornot
    logfile = fullfile(figpath,'analyses_figure_4.txt');
    if exist(logfile,'file'), delete(logfile), end
    diary(logfile)
end


cfg_fig = [];
cfg_fig.figpath = figpath;
cfg_fig.saveornot = saveornot;

rgb = lines;
rgb(3,:) = [.7 .7 .7];

alphathresh = 0.05;
xl = [-0.2 .8];
xt = -0.2:.2:.8;
xtl = {'-0.2' '0' '0.2' '0.4' '0.6' '0.8'};

condnam = {'PRE','POST','REF'};

% neural coding: SEQANG (+ compute mean cos|sin)
load(fullfile('./results/zdec_sensor/','zdec_seqang_16Hz.mat'),'zdec_seqang','tim_evt')
ndata = size(zdec_seqang,1);
xvec  = tim_evt;
nsmp  = length(xvec);
zdec = cell(ndata,3);
for idata = 1:ndata
    for icond = 1:3
        zdec{idata,icond} = mean([zdec_seqang{idata,icond,1};zdec_seqang{idata,icond,2}]);
    end
end
zdec_seqang = zdec; clearvars zdec

% neural coding: ABSLLR
load(fullfile('./results/zdec_sensor/','zdec_absllr_08Hz.mat'),'zdec_absllr')


%% SEQANG: compute time windows for averaging
% either percond or global
ifirst_seqang_percond = nan(3,1);
ilast_seqang_percond = nan(3,1);
winmin_seqang_percond = nan(3,1);
winmax_seqang_percond = nan(3,1);
idxwin_seqang_percond = false(3,nsmp);

for icond = 1:3
    yavg_seqang = mean(cell2mat(zdec_seqang(:,icond)));
    ymax_seqang = max(yavg_seqang);
    idxwin = yavg_seqang>=0.9*ymax_seqang;
    clust = bwconncomp(idxwin);
    ifirst_seqang_percond(icond) = clust.PixelIdxList{1}(1); % first cluster
    ilast_seqang_percond(icond) = clust.PixelIdxList{1}(end);
    clearvars idxwin
    winmin_seqang_percond(icond) = xvec(ifirst_seqang_percond(icond));
    winmax_seqang_percond(icond) = xvec(ilast_seqang_percond(icond));
    idxwin_seqang_percond(icond,ifirst_seqang_percond(icond):ilast_seqang_percond(icond)) = true;
    fprintf('SEQANG averaging time window %s: %d-%d ms\n',condnam{icond},round(1000*winmin_seqang_percond(icond)),round(1000*winmax_seqang_percond(icond)))
end

% global
ifirst_seqang_global = repmat(min(ifirst_seqang_percond),3,1);
ilast_seqang_global  = repmat(max(ilast_seqang_percond),3,1);
winmin_seqang_global = repmat(xvec(min(ifirst_seqang_percond)),3,1);
winmax_seqang_global = repmat(xvec(max(ilast_seqang_percond)),3,1);
idxwin_seqang_global = false(3,nsmp);
idxwin_seqang_global(:,min(ifirst_seqang_percond):max(ilast_seqang_percond)) = true;
fprintf('SEQANG averaging time window global: %d-%d ms\n',round(1000*winmin_seqang_global(1)),round(1000*winmax_seqang_global(1)))

%% SEQANG peak latency (jackknife)
% decoding latency PEAK (1 per condition)
fprintf('SEQANG: Jackknifed mean peak for each condition:\n')
ncond = 3;
xpeak_jk = nan(ndata,ncond);
ypeak_jk = nan(ndata,ncond);
for icond = 1:ncond
    x = cell2mat(zdec_seqang(:,icond));
    [xpeak_jk(:,icond),ypeak_jk(:,icond)] = jackknife_peak(x,xvec);
    fprintf('\t%s : xpeak = %5.5f ms\n',condnam{icond},1000*mean(xpeak_jk(:,icond)))
end

xpeak_seqang = mean(mean(xpeak_jk));
fprintf('\tMEAN LATENCY: xpeak = %5.5f ms\n',1000*mean(xpeak_seqang))

% repeated-measures ANOVA across conditions on jackknifed estimates:
% F corrected by (n-1)^2, p Greenhouse-Geisser corrected
jknam = {'LATENCY','AMPLITUDE'};
for ijk = 1:2
    if ijk == 1
        X = xpeak_jk;
    else
        X = ypeak_jk;
    end
    evalc('stat = anova_1way_ranova(X,''condition'');'); % do not print uncorrected F
    Fc = stat.F/(ndata-1)^2;
    pc = 1-fcdf(Fc,stat.eps*stat.df(1),stat.eps*stat.df(2));
    eta2p = Fc*stat.df(1)/(Fc*stat.df(1)+stat.df(2));
    fprintf('\tANOVA (jackknife-corrected) %s: F(%d,%d) = %.4f, p = %.4f (eps = %.2f, partial eta2 = %.3f)\n', ...
        jknam{ijk},stat.df(1),stat.df(2),Fc,pc,stat.eps,eta2p);
end

% POST-HOC PRE vs POST: latency
[~,~,~,STATS] = ttest(xpeak_jk(:,1),xpeak_jk(:,2));
t_jk = STATS.tstat/(ndata-1); % jackknife-corrected t
PVAL = 2*tcdf(-abs(t_jk),ndata-1);
bf10 = bf.ttest('T',t_jk,'N',ndata); % BF from the jackknife-corrected t
if PVAL < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
fprintf('\tPOST-HOC ttest %s vs %s LATENCY: t_jk(%d) = %+.5f ; p = %.5f ; %s = %.5f\n',condnam{1},condnam{2},STATS.df,t_jk,PVAL,BFstr,BFval)

% POST-HOC PRE vs POST: amplitude
[~,~,~,STATS] = ttest(ypeak_jk(:,1),ypeak_jk(:,2));
t_jk = STATS.tstat/(ndata-1); % jackknife-corrected t
PVAL = 2*tcdf(-abs(t_jk),ndata-1);
bf10 = bf.ttest('T',t_jk,'N',ndata); % BF from the jackknife-corrected t
if PVAL < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
fprintf('\tPOST-HOC ttest %s vs %s AMPLITUDE: t_jk(%d) = %+.5f ; p = %.5f ; %s = %.5f\n',condnam{1},condnam{2},STATS.df,t_jk,PVAL,BFstr,BFval)

%% SEQANG: plot time course
varnam = 'seqang';
win = [winmin_seqang_global(1) winmax_seqang_global(1)];
zdec_avg = zeros(nsmp,3);
zdec_err = zeros(nsmp,3);
for icond = 1:3
    zdec_avg(:,icond) = mean(cell2mat(zdec_seqang(:,icond)),'omitnan');
    zdec_err(:,icond) = std(cell2mat(zdec_seqang(:,icond)),'omitnan')/sqrt(ndata);
end

pbar = 1.3; HGT = 4;
cfg = cfg_fig;
cfg.figname = sprintf('%s_time_course',varnam);
cfg.xvec = xvec;
cfg.yavg = zdec_avg;
cfg.yerr = zdec_err;
cfg.win = win;
cfg.ylim = [-0.05,0.8];
cfg.ytick = 0:.2:.8;
cfg.peak = xpeak_seqang;
cfg.rgb = rgb;
cfg.pbar = pbar;
cfg.HGT = HGT;
helper_plot_decoding(cfg)

% plot inset with mean over the time window (bars: PRE and POST, line: REF)
idxwin = idxwin_seqang_global(1,:);
X = nan(ndata,3);
for icond = 1:3
    x = cell2mat(zdec_seqang(:,icond));
    X(:,icond) = mean(x(:,idxwin),2);
end

% ttest PRE vs POST
[H,P,CI,STATS] = ttest(X(:,1),X(:,2));
effd = STATS.tstat/sqrt(ndata);
bf10 = bf.ttest(X(:,1),X(:,2));
if H == 1, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
fprintf('TTEST ON MEAN TIME WINDOW %d-%d ms PRE-POST: t%02d = %+.5f ; p = %.5f ; Cohen''s d = %+.5f ; CI = [%+.5f %+.5f] ; %s = %.5f\n', ...
    round(win*1000),STATS.df,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);

figname = sprintf('%s_mean_%03d_%03d',varnam,round(win*1000));
figure('Color','white','Name',figname)
pbar = 1; HGT = 1.5;
hold on
bavg = mean(X);
berr = std(X)/sqrt(ndata);
xlim([0.4,2.6]);
ylim([0 .8]);
for icond = 1:2
    bar(icond,bavg(icond),0.5,'EdgeColor',rgb(icond,:),'FaceColor',0.5*(rgb(icond,:)+1),'LineWidth',1);
    plot([icond,icond],bavg(icond)+berr(icond)*[-1,+1],'k-');
end
plot(xlim,bavg(3)*[1,1],'Color',rgb(3,:))
set(gca,'XTick',1:2,'YTick',[0 .4 .8]);
helper_ILL(gcf,pbar,HGT,false);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

%% SEQANG: plot OPT/SUB bars
% neural coding: seqang OPT/SUB
load(fullfile('./results/zdec_sensor/zdec_seqang_16Hz_optsub.mat'),'zdec_seqang','tim_evt')
xvec = tim_evt;
% size(seqang) = 37 x 3 x 2 = ndata x ncond x cos/sin x opt/sub
zdec = cell(ndata,3,2); % compute mean cos/sin
for idata = 1:ndata
    for icond = 1:3
        for icor = 1:2
            zdec{idata,icond,icor} = mean([zdec_seqang{idata,icond,1,icor};zdec_seqang{idata,icond,2,icor}]);
        end
    end
end
zdec_seqang = zdec; clearvars zdec

idxwin = idxwin_seqang_global;
winmin = winmin_seqang_global;
winmax = winmax_seqang_global;

yl = [-0.05 0.8];

zdec_cond3 = squeeze(cell2mat(zdec_seqang(:,3,:)));
X3 = squeeze(mean(zdec_cond3(:,idxwin(3,:),:),2));
for icond = 1:2
    zdec_cond = squeeze(cell2mat(zdec_seqang(:,icond,:)));
    X = squeeze(mean(zdec_cond(:,idxwin(icond,:),:),2));
    cfg = cfg_fig;
    cfg.pbar = .6;
    cfg.HGT = 4;
    cfg.figname = sprintf('seqang_optsub_mean_%s_%03d_%03d',condnam{icond},round(winmin(icond)*1000),round(winmax(icond)*1000)); 
    cfg.rgb = [rgb(icond,:);rgb(icond,:);];
    cfg.x = X;
    cfg.plotx = true;
    cfg.xlabel = condnam{icond};
    cfg.labels = {'correct' 'error'};
    cfg.ttest = true; cfg.tail = 'right'; % one-tailed paired t-test OPT > SUB, shown on the plot
    cfg.ylabel = 'coding precision (z)';
    cfg.ylim = yl; cfg.ytick = 0:.4:.8;
    helper_plot_overall_bars(cfg);
    plot(xlim,mean(X3(:,1)).*[1 1],'Color',rgb(3,:),'LineWidth',1)
    plot(xlim,mean(X3(:,2)).*[1 1],'Color',rgb(3,:),'LineWidth',1,'LineStyle',':')
    helper_ILL(gcf,cfg.pbar,cfg.HGT,3);
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end

% ONE-TAILED TTEST: COR>ERR
for icond = 1:3
    zdec_cond = squeeze(cell2mat(zdec_seqang(:,icond,:)));
    X = squeeze(mean(zdec_cond(:,idxwin(icond,:),:),2));
    diffX = X(:,1)-X(:,2);
    [H,P,CI,STATS]=ttest(diffX,0,0.05,'right');
    fprintf('\t%s one-tailed ttest OPT>SUB ? : H = %d ; t = %+09.5f ; p = %5.5f ; CI = [%+5.5f %+5.5f]\n',condnam{icond},H,STATS.tstat,P,CI(1),CI(2));
end

%% SEQANG: OPT/SUB plot time course (supp)
load(fullfile('./results/zdec_sensor/zdec_seqang_16Hz_optsub.mat'),'zdec_seqang','tim_evt')
xvec = tim_evt;
% size(seqang) = 37 x 3 x 2 = ndata x ncond x cos/sin x opt/sub
zdec = cell(ndata,3,2); % compute mean cos/sin
for idata = 1:ndata
    for icond = 1:3
        for icor = 1:2
            zdec{idata,icond,icor} = mean([zdec_seqang{idata,icond,1,icor};zdec_seqang{idata,icond,2,icor}]);
        end
    end
end
zdec_seqang = zdec; clearvars zdec

varnam = 'seqang';
win = [winmin_seqang_percond(1) winmax_seqang_percond(1)]; % PRE only
zdec_avg_opt = zeros(nsmp,3);
zdec_err_opt = zeros(nsmp,3);
zdec_avg_sub = zeros(nsmp,3);
zdec_err_sub = zeros(nsmp,3);
for icond = 1:3
    zdec_avg_opt(:,icond) = mean(cell2mat(zdec_seqang(:,icond,1)),'omitnan');
    zdec_err_opt(:,icond) = std(cell2mat(zdec_seqang(:,icond,1)),'omitnan')/sqrt(ndata);
    zdec_avg_sub(:,icond) = mean(cell2mat(zdec_seqang(:,icond,2)),'omitnan');
    zdec_err_sub(:,icond) = std(cell2mat(zdec_seqang(:,icond,2)),'omitnan')/sqrt(ndata);
end

pbar = 1.3; HGT = 4;
cfg = cfg_fig;
cfg.figname = sprintf('%s_time_course_OPTSUB',varnam);
cfg.xvec = xvec;
cfg.yavg = zdec_avg_opt;
cfg.yerr = zdec_err_opt;
cfg.yavg_dashed = zdec_avg_sub;
cfg.yerr_dashed = zdec_err_sub;
cfg.win = win;
cfg.ylim = [-0.05,0.8];
cfg.ytick = 0:.2:.8;
cfg.rgb = rgb;
cfg.pbar = pbar;
cfg.HGT = HGT;
helper_plot_decoding(cfg)

%% ABSLLR: compute time windows for averaging
% either percond or global
ifirst_absllr_percond = nan(3,1);
ilast_absllr_percond = nan(3,1);
winmin_absllr_percond = nan(3,1);
winmax_absllr_percond = nan(3,1);
idxwin_absllr_percond = false(3,nsmp);

for icond = [1 3] % cond POST not decoded
    yavg_absllr = mean(cell2mat(zdec_absllr(:,icond)));
    ymax_absllr = max(yavg_absllr);
    idxwin = yavg_absllr>=0.9*ymax_absllr;
    clust = bwconncomp(idxwin);
    ifirst_absllr_percond(icond) = clust.PixelIdxList{1}(1); % first cluster
    ilast_absllr_percond(icond) = clust.PixelIdxList{1}(end);
    clearvars idxwin
    winmin_absllr_percond(icond) = xvec(ifirst_absllr_percond(icond));
    winmax_absllr_percond(icond) = xvec(ilast_absllr_percond(icond));
    idxwin_absllr_percond(icond,ifirst_absllr_percond(icond):ilast_absllr_percond(icond)) = true;
    fprintf('ABSLLR averaging time window %s: %d-%d ms\n',condnam{icond},round(1000*winmin_absllr_percond(icond)),round(1000*winmax_absllr_percond(icond))) 
end

% global
ifirst_absllr_global = repmat(min(ifirst_absllr_percond),3,1);
ilast_absllr_global  = repmat(max(ilast_absllr_percond),3,1);
winmin_absllr_global = repmat(xvec(min(ifirst_absllr_percond)),3,1);
winmax_absllr_global = repmat(xvec(max(ilast_absllr_percond)),3,1);
idxwin_absllr_global = false(3,nsmp);
idxwin_absllr_global(:,min(ifirst_absllr_percond):max(ilast_absllr_percond)) = true;
fprintf('ABSLLR averaging time window global (PRE and REF): %d-%d ms\n',round(1000*winmin_absllr_global(1)),round(1000*winmax_absllr_global(1)))

%% ABSLLR: peak latency (jackknife)
% decoding latency PEAK (1 per condition)
fprintf('ABSLLR: Jackknifed mean peak for each condition:\n')
ncond = 3;
xpeak_jk = nan(ndata,ncond);
xpeak = nan(1,ncond);
ypeak_jk = nan(ndata,ncond);
for icond = [1 3]
    x = cell2mat(zdec_absllr(:,icond));
    [xpeak_jk(:,icond),ypeak_jk(:,icond)] = jackknife_peak(x,xvec);
    fprintf('\t%s : xpeak = %5.5f ms\n',condnam{icond},1000*mean(xpeak_jk(:,icond)))
end

xpeak_absllr = mean(mean(xpeak_jk(:,[1 3])));
fprintf('\tMEAN LATENCY: xpeak = %5.5f ms\n',1000*mean(xpeak_absllr))

% POST-HOC PRE vs REF: latency
[~,~,~,STATS] = ttest(xpeak_jk(:,1),xpeak_jk(:,3));
t_jk = STATS.tstat/(ndata-1); % jackknife-corrected t
PVAL = 2*tcdf(-abs(t_jk),ndata-1);
bf10 = bf.ttest('T',t_jk,'N',ndata); % BF from the jackknife-corrected t
if PVAL < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
fprintf('\tPOST-HOC ttest %s vs %s LATENCY: t_jk(%d) = %+.5f ; p = %.5f ; %s = %.5f\n',condnam{1},condnam{3},STATS.df,t_jk,PVAL,BFstr,BFval)

% POST-HOC PRE vs REF: amplitude
[~,~,~,STATS] = ttest(ypeak_jk(:,1),ypeak_jk(:,3));
t_jk = STATS.tstat/(ndata-1); % jackknife-corrected t
PVAL = 2*tcdf(-abs(t_jk),ndata-1);
bf10 = bf.ttest('T',t_jk,'N',ndata); % BF from the jackknife-corrected t
if PVAL < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
fprintf('\tPOST-HOC ttest %s vs %s AMPLITUDE: t_jk(%d) = %+.5f ; p = %.5f ; %s = %.5f\n',condnam{1},condnam{3},STATS.df,t_jk,PVAL,BFstr,BFval)


%% ABSLLR: plot time course
varnam = 'absllr';
win = [winmin_absllr_global(1) winmax_absllr_global(1)];
zdec_avg = zeros(nsmp,3);
zdec_err = zeros(nsmp,3);
for icond = 1:3
    zdec_avg(:,icond) = mean(cell2mat(zdec_absllr(:,icond)),'omitnan');
    zdec_err(:,icond) = std(cell2mat(zdec_absllr(:,icond)),'omitnan')/sqrt(ndata);
end

pbar = 1.3; HGT = 4;
cfg = cfg_fig;
cfg.figname = sprintf('%s_time_course',varnam);
cfg.xvec = xvec;
cfg.yavg = zdec_avg;
cfg.yerr = zdec_err;
cfg.win = win;
cfg.peak = xpeak_absllr;
cfg.rgb = rgb;
cfg.pbar = pbar;
cfg.HGT = HGT;
helper_plot_decoding(cfg)

%% ABSLLR: decoding against 0 in each condition (cluster-based permutation test across time)
cfg_stat = [];
cfg_stat.xvec = xvec;
for icond = 1:3
    fprintf('<strong>ABSLLR %s: decoding against 0, cluster-based permutation test</strong>\n',condnam{icond})
    stat_absllr(icond) = cluster_permutation(cell2mat(zdec_absllr(:,icond)),[],cfg_stat);
end

%% ABSLLR: plot OPT/SUB bars
% neural coding: ABSLLR OPT/SUB
load(fullfile('./results/zdec_sensor/zdec_absllr_08Hz_optsub.mat'),'zdec_absllr','tim_evt')
% size(absllr) = 37 x 3 x 2 = ndata x ncond x opt/sub

idxwin = idxwin_absllr_global;
winmin = winmin_absllr_global;
winmax = winmax_absllr_global;

yl = [-0.025 0.4];

zdec_cond3 = squeeze(cell2mat(zdec_absllr(:,3,:)));
X3 = squeeze(mean(zdec_cond3(:,idxwin(3,:),:),2));
for icond = 1:2
    zdec_cond = squeeze(cell2mat(zdec_absllr(:,icond,:)));
    X = squeeze(mean(zdec_cond(:,idxwin(icond,:),:),2));

    cfg = cfg_fig;
    cfg.pbar = .6;
    cfg.HGT = 4;
    cfg.figname = sprintf('absllr_optsub_mean_%s_%03d_%03d',condnam{icond},round(winmin(icond)*1000),round(winmax(icond)*1000));
    cfg.rgb = [rgb(icond,:);rgb(icond,:);];
    cfg.x = X;
    cfg.plotx = true;
    cfg.xlabel = condnam{icond};
    cfg.labels = {'correct' 'error'};
    cfg.ttest = true; cfg.tail = 'right'; % one-tailed paired t-test OPT > SUB, shown on the plot
    cfg.ylabel = 'coding precision (z)';
    cfg.ylim = yl; cfg.ytick = 0:.2:.4;
    helper_plot_overall_bars(cfg);
    plot(xlim,mean(X3(:,1)).*[1 1],'Color',rgb(3,:),'LineWidth',1)
    plot(xlim,mean(X3(:,2)).*[1 1],'Color',rgb(3,:),'LineWidth',1,'LineStyle',':')
    helper_ILL(gcf,cfg.pbar,cfg.HGT,3);
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end

% ONE-TAILED TTEST: COR>ERR
for icond = [1 2 3]
    zdec_cond = squeeze(cell2mat(zdec_absllr(:,icond,:)));
    X = squeeze(mean(zdec_cond(:,idxwin(icond,:),:),2));
    diffX = X(:,1)-X(:,2);
    [H,P,CI,STATS]=ttest(diffX,0,0.05,'right');
    fprintf('\t%s one-tailed ttest OPT>SUB ? : H = %d ; t = %+09.5f ; p = %5.5f ; CI = [%+5.5f %+5.5f]\n',condnam{icond},H,STATS.tstat,P,CI(1),CI(2));
end

%% ABSLLR: OPT/SUB plot time course (supp)
% neural coding: ABSLLR OPT/SUB
load(fullfile('./results/zdec_sensor/zdec_absllr_08Hz_optsub.mat'),'zdec_absllr','tim_evt')
% size(absllr) = 37 x 3 x 2 = ndata x ncond x opt/sub
xvec = tim_evt;

varnam = 'absllr';
win = [winmin_absllr_percond(1) winmax_absllr_percond(1)];
zdec_avg_opt = zeros(nsmp,3);
zdec_err_opt = zeros(nsmp,3);
zdec_avg_sub = zeros(nsmp,3);
zdec_err_sub = zeros(nsmp,3);
for icond = 1:3
    zdec_avg_opt(:,icond) = mean(cell2mat(zdec_absllr(:,icond,1)),'omitnan');
    zdec_err_opt(:,icond) = std(cell2mat(zdec_absllr(:,icond,1)),'omitnan')/sqrt(ndata);
    zdec_avg_sub(:,icond) = mean(cell2mat(zdec_absllr(:,icond,2)),'omitnan');
    zdec_err_sub(:,icond) = std(cell2mat(zdec_absllr(:,icond,2)),'omitnan')/sqrt(ndata);
end

pbar = 1.3; HGT = 4;
cfg = cfg_fig;
cfg.figname = sprintf('%s_time_course_OPTSUB',varnam);
cfg.xvec = xvec;
cfg.yavg = zdec_avg_opt;
cfg.yerr = zdec_err_opt;
cfg.yavg_dashed = zdec_avg_sub;
cfg.yerr_dashed = zdec_err_sub;
cfg.win = win;
cfg.ylim = [-0.025,0.4];
cfg.ytick = 0:.1:.4;
cfg.rgb = rgb;
cfg.pbar = pbar;
cfg.HGT = HGT;
helper_plot_decoding(cfg)

%% ABSLLR: plot OPT/SUB bars (main: GLOBAL time window; supp: PRE time window)
% neural coding: ABSLLR OPT/SUB
load(fullfile('./results/zdec_sensor/zdec_absllr_08Hz_optsub.mat'),'zdec_absllr','tim_evt')
% size(absllr) = 37 x 3 x 2 = ndata x ncond x opt/sub
xvec = tim_evt;

idxwin = idxwin_absllr_global(1,:);
winmin = winmin_absllr_global(1);
winmax = winmax_absllr_global(1);

yl = [-0.025 0.4];

zdec_cond3 = squeeze(cell2mat(zdec_absllr(:,3,:)));
X3 = squeeze(mean(zdec_cond3(:,idxwin_absllr_percond(3,:),:),2));
for icond = 1:2
    zdec_cond = squeeze(cell2mat(zdec_absllr(:,icond,:)));
    X = squeeze(mean(zdec_cond(:,idxwin,:),2));

    cfg = cfg_fig;
    cfg.pbar = .6;
    cfg.HGT = 4;
    cfg.figname = sprintf('absllr_optsub_mean_%s_%03d_%03d',condnam{icond},round(winmin*1000),round(winmax*1000));
    cfg.rgb = [rgb(icond,:);rgb(icond,:);];
    cfg.x = X;
    cfg.plotx = true;
    cfg.xlabel = condnam{icond};
    cfg.labels = {'correct' 'error'};
    cfg.ttest = true; cfg.tail = 'right';
    cfg.ylabel = 'coding precision (z)';
    cfg.ylim = yl; cfg.ytick = 0:.2:.4;
    helper_plot_overall_bars(cfg);
    plot(xlim,mean(X3(:,1)).*[1 1],'Color',rgb(3,:),'LineWidth',1)
    plot(xlim,mean(X3(:,2)).*[1 1],'Color',rgb(3,:),'LineWidth',1,'LineStyle',':')
    helper_ILL(gcf,cfg.pbar,cfg.HGT,3);
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end

% ONE-TAILED TTEST: COR>ERR
idxwin = idxwin_absllr_percond;
idxwin(2,:) = idxwin(1,:);
for icond = [1 2 3]
    zdec_cond = squeeze(cell2mat(zdec_absllr(:,icond,:)));
    X = squeeze(mean(zdec_cond(:,idxwin(icond,:),:),2));
    diffX = X(:,1)-X(:,2);
    [H,P,CI,STATS]=ttest(diffX,0,0.05,'right');
    fprintf('\t%s one-tailed ttest OPT>SUB ? : H = %d ; t = %+09.5f ; p = %5.5f ; CI = [%+5.5f %+5.5f]\n',condnam{icond},H,STATS.tstat,P,CI(1),CI(2));
end

%% SEQANG and ABSLLR SOURCE: AVERAGE (main) and OPT/SUB difference (supp)
%   * AVERAGE (main): decoding precision on the global window
%   * OPT/SUB difference (supp): correct-error
if plotsource
    srcpath = './results/zdec_source';
    cortex = load('./default_cortex.mat'); % cortex of the default anatomy
    addpath(fullfile(cfg_path.brainstorm,'toolbox','misc')) % Brainstorm colormaps

    % AVERAGE (main): decoding precision
    clim_seqang = [0 0.5];
    clim_absllr = [0 0.1];
    % seqang: template vertex x time sample x participant x condition x cos/sin
    tmp = load(fullfile(srcpath,'zdec_seqang_16Hz.mat'),'zdec','tim_evt');
    if max(abs(tmp.tim_evt-xvec)) > 1e-6
        error('not the same time samples at the sensor and source levels!')
    end
    for icond = [2 1]
        iwin = round([winmin_seqang_global(icond)*1000,winmax_seqang_global(icond)*1000]);
        x = tmp.zdec(:,idxwin_seqang_global(icond,:),:,icond,:);
        cfg = cfg_fig;
        cfg.cortex = cortex;
        cfg.map  = mean(mean(mean(x,5),3),2); % mean over cos/sin, participants and time window
        cfg.clim = clim_seqang;
        cfg.cmap = cmap_royal_gramma(256);
        cfg.views = {'top','left'};
        cfg.figname = sprintf('seqang_time_avg_%s_%03d-%03d',condnam{icond},iwin(1),iwin(2));
        helper_plot_source(cfg)
    end
    clearvars tmp x
    % absllr: template vertex x time sample x participant x condition (POST: not decoded)
    tmp = load(fullfile(srcpath,'zdec_absllr_08Hz.mat'),'zdec','tim_evt');
    if max(abs(tmp.tim_evt-xvec)) > 1e-6
        error('not the same time samples at the sensor and source levels!')
    end
    for icond = [2 1]
        iwin = round([winmin_absllr_global(icond)*1000,winmax_absllr_global(icond)*1000]);
        x = tmp.zdec(:,idxwin_absllr_global(icond,:),:,icond);
        cfg = cfg_fig;
        cfg.cortex = cortex;
        cfg.map  = mean(mean(x,3),2); % mean over participants and time window
        cfg.clim = clim_absllr;
        cfg.cmap = cmap_royal_gramma(256);
        cfg.views = {'top','left'};
        cfg.figname = sprintf('absllr_time_avg_%s_%03d-%03d',condnam{icond},iwin(1),iwin(2));
        helper_plot_source(cfg)
    end
    clearvars tmp x

    % OPT/SUB difference (supp)
    newmin = -0.05;
    newmax = 0.05;
    % seqang: template vertex x time sample x participant x condition x cos/sin x opt/sub
    tmp = load(fullfile(srcpath,'zdec_seqang_16Hz_optsub.mat'),'zdec','tim_evt');
    if max(abs(tmp.tim_evt-xvec)) > 1e-6
        error('not the same time samples at the sensor and source levels!')
    end
    idxwin = idxwin_seqang_percond;
    winmin = winmin_seqang_percond;
    winmax = winmax_seqang_percond;
    for icond = [2 1]
        iwin = round([winmin(icond)*1000,winmax(icond)*1000]);
        x = tmp.zdec(:,idxwin(icond,:),:,icond,:,1)-tmp.zdec(:,idxwin(icond,:),:,icond,:,2);
        cfg = cfg_fig;
        cfg.cortex = cortex;
        cfg.map  = mean(mean(mean(x,5),3),2); % mean over cos/sin, participants and time window
        cfg.clim = [newmin newmax];
        cfg.cmap = cmap_ns_grey(256);
        cfg.views = {'top','left'};
        cfg.figname = sprintf('seqang_optsub_time_avg_%s_%03d-%03d',condnam{icond},iwin(1),iwin(2));
        helper_plot_source(cfg)
    end
    clearvars tmp x
    % absllr: template vertex x time sample x participant x condition x opt/sub
    tmp = load(fullfile(srcpath,'zdec_absllr_08Hz_optsub.mat'),'zdec','tim_evt');
    if max(abs(tmp.tim_evt-xvec)) > 1e-6
        error('not the same time samples at the sensor and source levels!')
    end
    idxwin = idxwin_absllr_percond;
    winmin = winmin_absllr_percond;
    winmax = winmax_absllr_percond;
    icond = 1; % PRE
    iwin = round([winmin(icond)*1000,winmax(icond)*1000]);
    x = tmp.zdec(:,idxwin(icond,:),:,icond,1)-tmp.zdec(:,idxwin(icond,:),:,icond,2);
    cfg = cfg_fig;
    cfg.cortex = cortex;
    cfg.map  = mean(mean(x,3),2); % mean over participants and time window
    cfg.clim = [newmin newmax];
    cfg.cmap = cmap_ns_grey(256);
    cfg.views = {'top','left'};
    cfg.figname = sprintf('absllr_optsub_time_avg_%s_%03d-%03d',condnam{icond},iwin(1),iwin(2));
    helper_plot_source(cfg)
    clearvars tmp x
end

if saveornot
    diary off
end
