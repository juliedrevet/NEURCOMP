%% analyses_figure_5
% this script runs the analyses to reproduce the main Figure 5 (and corresponding supp. figures)
%   (a) conllr: time course, peak latency (jackknife), PRE vs POST (cluster-based permutation test), inset PRE vs POST
%   (a) nulllr: time course, peak latency (jackknife), PRE vs POST (cluster-based permutation test), inset PRE vs POST
%   (b) compdif: time course, bars (conllr global time window)
%   (c) comprel: ROI-specific compression in PRE (source level)
%   (d) compdif: OPT/SUB time course, bars (conllr time window of each condition)
%   supp: conllr and nulllr OPT/SUB time courses and bars
%   supp: conllr, nulllr and compdif sources, conllr and nulllr OPT/SUB sources
%   supp: conllr and nulllr time courses in the ROIs (SPL, ITC, LPFC) in PRE, source level
% conllr and nulllr are decoded from the same particles (particle filter
% resampled according to the participant's choice)

clearvars
close all
clc

plotsource = false;
saveornot = false;

addpath('..'); cfg_path = config_paths;

figpath = './figs';
[~,~] = mkdir(figpath);

if saveornot
logfile = fullfile(figpath,'analyses_figure_5.txt');
if exist(logfile,'file'), delete(logfile), end
diary(logfile)
end


cfg_fig = [];
cfg_fig.figpath = figpath;
cfg_fig.saveornot = saveornot;

rgb = lines;
rgb(3,:) = [.7 .7 .7];

condnam = {'PRE','POST','REF'};
ordercond = [2 1]; % bars: POST then PRE

% neural coding: CONLLR and NULLLR (ndata x ncond)
load('./results/zdec_sensor/zdec_conllr_particle_08Hz.mat','zdec_conllr','tim_evt')
load('./results/zdec_sensor/zdec_nulllr_particle_08Hz.mat','zdec_nulllr')
ndata = size(zdec_conllr,1);
xvec  = tim_evt;
nsmp  = length(xvec);
% compression: conllr - nulllr
zdec_compdif = cellfun(@minus,zdec_conllr,zdec_nulllr,'UniformOutput',false);

% neural coding OPT/SUB (ndata x ncond x opt/sub)
tmp = load('./results/zdec_sensor/zdec_conllr_particle_08Hz_optsub.mat','zdec_conllr');
zdec_conllr_optsub = tmp.zdec_conllr;
tmp = load('./results/zdec_sensor/zdec_nulllr_particle_08Hz_optsub.mat','zdec_nulllr');
zdec_nulllr_optsub = tmp.zdec_nulllr; clearvars tmp
zdec_compdif_optsub = cellfun(@minus,zdec_conllr_optsub,zdec_nulllr_optsub,'UniformOutput',false);

% neural coding at the source level: CONLLR and NULLLR projected on the default
% anatomy (template vertex x time sample x participant x condition)
if plotsource
    srcpath = './results/zdec_source';
    tmp = load(fullfile(srcpath,'zdec_conllr_particle_08Hz.mat'),'zdec','tim_evt','subjlist');
    zsrc_conllr = tmp.zdec;
    if ~isequal(tmp.subjlist(:),(0:ndata-1)') || max(abs(tmp.tim_evt-tim_evt)) > 1e-6
        error('not the same participants or time samples at the sensor and source levels!')
    end
    tmp = load(fullfile(srcpath,'zdec_nulllr_particle_08Hz.mat'),'zdec');
    zsrc_nulllr = tmp.zdec; clearvars tmp
    % cortex of the default anatomy: Vertices, Faces, Atlas (regions of interest), SulciMap
    cortex = load('./default_cortex.mat');
    addpath(fullfile(cfg_path.brainstorm,'toolbox','misc')) % Brainstorm colormaps
end

%% CONLLR: compute time windows for averaging
% either percond or global
ifirst_conllr_percond = nan(3,1);
ilast_conllr_percond = nan(3,1);
winmin_conllr_percond = nan(3,1);
winmax_conllr_percond = nan(3,1);
idxwin_conllr_percond = false(3,nsmp);

for icond = 1:3
    yavg = mean(cell2mat(zdec_conllr(:,icond)));
    ymax = max(yavg);
    idxwin = yavg>=0.9*ymax;
    clust = bwconncomp(idxwin);
    ifirst_conllr_percond(icond) = clust.PixelIdxList{1}(1); % first cluster
    ilast_conllr_percond(icond) = clust.PixelIdxList{1}(end);
    clearvars idxwin
    winmin_conllr_percond(icond) = xvec(ifirst_conllr_percond(icond));
    winmax_conllr_percond(icond) = xvec(ilast_conllr_percond(icond));
    idxwin_conllr_percond(icond,ifirst_conllr_percond(icond):ilast_conllr_percond(icond)) = true;
    fprintf('CONLLR averaging time window %s: %d-%d ms\n',condnam{icond},round(1000*winmin_conllr_percond(icond)),round(1000*winmax_conllr_percond(icond))) 
end

% global
ifirst_conllr_global = repmat(min(ifirst_conllr_percond),3,1);
ilast_conllr_global  = repmat(max(ilast_conllr_percond),3,1);
winmin_conllr_global = repmat(xvec(min(ifirst_conllr_percond)),3,1);
winmax_conllr_global = repmat(xvec(max(ilast_conllr_percond)),3,1);
idxwin_conllr_global = false(3,nsmp);
idxwin_conllr_global(:,min(ifirst_conllr_percond):max(ilast_conllr_percond)) = true;
fprintf('CONLLR averaging time window global: %d-%d ms\n',round(1000*winmin_conllr_global(1)),round(1000*winmax_conllr_global(1)))

%% CONLLR and NULLLR: peak latency (jackknife)
varnam = {'conllr','nulllr'};
zdec_all = {zdec_conllr,zdec_nulllr};
xpeak_avg = nan(1,2); % mean jackknifed latency across conditions (plots)
jknam = {'LATENCY','AMPLITUDE'};
for ivar = 1:2
    fprintf('%s: Jackknifed mean peak for each condition:\n',upper(varnam{ivar}))
    xpeak_jk = nan(ndata,3);
    ypeak_jk = nan(ndata,3);
    for icond = 1:3
        [xpeak_jk(:,icond),ypeak_jk(:,icond)] = jackknife_peak(cell2mat(zdec_all{ivar}(:,icond)),xvec);
        fprintf('\t%s : xpeak = %5.5f ms\n',condnam{icond},1000*mean(xpeak_jk(:,icond)))
    end
    xpeak_avg(ivar) = mean(xpeak_jk(:));
    fprintf('\tMEAN LATENCY: xpeak = %5.5f ms\n',1000*xpeak_avg(ivar))

    % repeated-measures ANOVA across conditions on jackknifed estimates:
    % F corrected by (n-1)^2, p Greenhouse-Geisser corrected
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

        % POST-HOC PRE vs POST
        [~,~,~,STATS] = ttest(X(:,1),X(:,2));
        t_jk = STATS.tstat/(ndata-1); % jackknife-corrected t
        PVAL = 2*tcdf(-abs(t_jk),ndata-1);
        bf10 = bf.ttest('T',t_jk,'N',ndata); % BF from the jackknife-corrected t
        if PVAL < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
        fprintf('\tPOST-HOC ttest %s vs %s %s: t_jk(%d) = %+.5f ; p = %.5f ; %s = %.5f\n',condnam{1},condnam{2},jknam{ijk},STATS.df,t_jk,PVAL,BFstr,BFval)
    end
end

%% CONLLR, NULLLR, COMPDIF: plot time courses (conllr global time window)
varnam = {'conllr','nulllr','compdif'};
zdec_all = {zdec_conllr,zdec_nulllr,zdec_compdif};
ylab = {'coding precision (z)','coding precision (z)','coding difference (z)'};
win = [winmin_conllr_global(1) winmax_conllr_global(1)];
for ivar = 1:3
    zdec_avg = zeros(nsmp,3);
    zdec_err = zeros(nsmp,3);
    for icond = 1:3
        x = cell2mat(zdec_all{ivar}(:,icond));
        zdec_avg(:,icond) = mean(x,'omitnan');
        zdec_err(:,icond) = std(x,'omitnan')/sqrt(ndata);
    end
    cfg = cfg_fig;
    cfg.figname = sprintf('%s_time_course',varnam{ivar});
    cfg.xvec = xvec;
    cfg.yavg = zdec_avg;
    cfg.yerr = zdec_err;
    cfg.win = win;
    if ivar <= 2
        cfg.peak = xpeak_avg(ivar);
    end
    cfg.ylim = [-0.025,0.4];
    cfg.ytick = 0:.1:.4;
    cfg.ylabel = ylab{ivar};
    cfg.rgb = rgb;
    cfg.pbar = 1.3;
    cfg.HGT = 4;
    helper_plot_decoding(cfg)
end

%% CONLLR and NULLLR: PRE vs POST (cluster-based permutation test across time on PRE - POST)
varnam = {'conllr','nulllr'};
zdec_all = {zdec_conllr,zdec_nulllr};
cfg_stat = [];
cfg_stat.xvec = xvec;
for ivar = 1:2
    fprintf('<strong>%s: %s - %s, cluster-based permutation test</strong>\n',upper(varnam{ivar}),condnam{1},condnam{2})
    x = cell2mat(zdec_all{ivar}(:,1))-cell2mat(zdec_all{ivar}(:,2));
    stat_prepost.(varnam{ivar}) = cluster_permutation(x,[],cfg_stat);
end

%% CONLLR and NULLLR: insets with mean over the conllr global time window (bars: POST and PRE, line: REF)
varnam = {'conllr','nulllr'};
zdec_all = {zdec_conllr,zdec_nulllr};
idxwin = idxwin_conllr_global(1,:);
win = [winmin_conllr_global(1) winmax_conllr_global(1)];
for ivar = 1:2
    X = nan(ndata,3);
    for icond = 1:3
        x = cell2mat(zdec_all{ivar}(:,icond));
        X(:,icond) = mean(x(:,idxwin),2);
    end

    % ttest PRE vs POST
    [H,P,CI,STATS] = ttest(X(:,1),X(:,2));
    effd = STATS.tstat/sqrt(ndata);
    bf10 = bf.ttest(X(:,1),X(:,2));
    if H == 1, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
    fprintf('%s: TTEST ON MEAN TIME WINDOW %d-%d ms PRE-POST: t%02d = %+.5f ; p = %.5f ; Cohen''s d = %+.5f ; CI = [%+.5f %+.5f] ; %s = %.5f\n', ...
        upper(varnam{ivar}),round(win*1000),STATS.df,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);

    figname = sprintf('%s_mean_%03d_%03d',varnam{ivar},round(win*1000));
    figure('Color','white','Name',figname)
    pbar = 1; HGT = 1.5;
    hold on
    bavg = mean(X);
    berr = std(X)/sqrt(ndata);
    xlim([0.4,2.6]);
    ylim([0 .4]);
    for ibar = 1:2
        icond = ordercond(ibar);
        bar(ibar,bavg(icond),0.5,'EdgeColor',rgb(icond,:),'FaceColor',0.5*(rgb(icond,:)+1),'LineWidth',1);
        plot([ibar,ibar],bavg(icond)+berr(icond)*[-1,+1],'k-');
    end
    plot(xlim,bavg(3)*[1,1],'Color',rgb(3,:))
    set(gca,'XTick',1:2,'YTick',[0 .2 .4]);
    helper_ILL(gcf,pbar,HGT,false);
    if saveornot
        print(fullfile(figpath,figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',figname,figpath)
    end
end

%% COMPDIF: bars (conllr global time window)
idxwin = idxwin_conllr_global(1,:);
win = [winmin_conllr_global(1) winmax_conllr_global(1)];
X = nan(ndata,3);
for icond = 1:3
    x = cell2mat(zdec_compdif(:,icond));
    X(:,icond) = mean(x(:,idxwin),2);
end

% ttest compression vs 0 in each condition and PRE vs POST
fprintf('COMPDIF (conllr - nulllr) ON MEAN TIME WINDOW %d-%d ms:\n',round(win*1000))
for icond = 1:3
    [~,P,CI,STATS] = ttest(X(:,icond));
    bf10 = bf.ttest(X(:,icond));
    if P < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
    fprintf('\t%s vs 0: mean = %+.5f ; t%02d = %+.5f ; p = %.5f ; Cohen''s d = %+.5f ; CI = [%+.5f %+.5f] ; %s = %.5f\n', ...
        condnam{icond},mean(X(:,icond)),STATS.df,STATS.tstat,P,STATS.tstat/sqrt(ndata),CI(1),CI(2),BFstr,BFval);
end
[~,P,CI,STATS] = ttest(X(:,1),X(:,2));
bf10 = bf.ttest(X(:,1),X(:,2));
if P < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
fprintf('\t%s vs %s: t%02d = %+.5f ; p = %.5f ; Cohen''s d = %+.5f ; CI = [%+.5f %+.5f] ; %s = %.5f\n', ...
    condnam{1},condnam{2},STATS.df,STATS.tstat,P,STATS.tstat/sqrt(ndata),CI(1),CI(2),BFstr,BFval);

% bars: POST and PRE, line: REF
cfg = cfg_fig;
cfg.figname = sprintf('compdif_mean_%03d_%03d',round(win*1000));
cfg.x = X(:,ordercond);
cfg.plotx = true;
cfg.rgb = rgb(ordercond,:);
cfg.labels = condnam(ordercond);
cfg.xlabel = 'condition';
cfg.ylabel = 'coding difference (z)';
cfg.ylim = [-.2 .4]; cfg.ytick = -.2:.2:.4;
cfg.pbar = .5; cfg.HGT = 4;
helper_plot_overall_bars(cfg);
plot(xlim,mean(X(:,3))*[1 1],'Color',rgb(3,:),'LineWidth',1)
helper_ILL(gcf,cfg.pbar,cfg.HGT,3);
if cfg.saveornot
    print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
end

%% COMPREL: ROI-specific compression in PRE (source level)
% compression = 1 - nulllr/conllr, computed on the leave-one-out means of
% conllr and nulllr averaged over the ROI (mean of left and right
% hemispheres) and over the conllr global time window
if plotsource
    zdec_conllr_source = zsrc_conllr; % template vertex x time sample x participant x condition
    zdec_nulllr_source = zsrc_nulllr;
    itime = idxwin_conllr_global(1,:);

    % ROIs from the Desikan-Killiany atlas of the default anatomy
    Atlas = cortex.Atlas;
    iatlas = find(contains({Atlas.Name},'Desikan-Killiany') & ~contains({Atlas.Name},'volume'));
    Scout = Atlas(iatlas).Scouts;
    roinam = {'rostralmiddlefrontal','inferiortemporal','superiorparietal'};
    roilab = {'LPFC','ITC','SPL'};
    nroi = numel(roinam);

    icond = 1; % PRE
    conllr_roi = nan(ndata,nroi);
    nulllr_roi = nan(ndata,nroi);
    for iroi = 1:nroi
        ivL = Scout(strcmp({Scout.Label},sprintf('%s L',roinam{iroi}))).Vertices;
        ivR = Scout(strcmp({Scout.Label},sprintf('%s R',roinam{iroi}))).Vertices;
        % mean over vertices of each hemisphere and time window, then mean of left and right
        conllr_roi(:,iroi) = (squeeze(mean(mean(zdec_conllr_source(ivL,itime,:,icond),1),2)) + ...
            squeeze(mean(mean(zdec_conllr_source(ivR,itime,:,icond),1),2)))/2;
        nulllr_roi(:,iroi) = (squeeze(mean(mean(zdec_nulllr_source(ivL,itime,:,icond),1),2)) + ...
            squeeze(mean(mean(zdec_nulllr_source(ivR,itime,:,icond),1),2)))/2;
    end

    % jackknifed compression
    comprel_jk = nan(ndata,nroi);
    for idata = 1:ndata
        subset = setdiff(1:ndata,idata);
        comprel_jk(idata,:) = 1-mean(nulllr_roi(subset,:),1)./mean(conllr_roi(subset,:),1);
    end
    comprel_avg = mean(comprel_jk,1);
    comprel_err = std(comprel_jk,[],1)*(ndata-1)/sqrt(ndata); % jackknife s.e.m.
    fprintf('COMPREL %s on %d-%d ms:\n',condnam{icond},round(1000*[winmin_conllr_global(1) winmax_conllr_global(1)]))
    for iroi = 1:nroi
        fprintf('\t%s: compression = %.1f ± %.1f %%\n',roilab{iroi},100*comprel_avg(iroi),100*comprel_err(iroi))
    end
    % pairwise comparisons between ROIs (jackknife-corrected t)
    pairs = nchoosek(1:nroi,2);
    for ipair = 1:size(pairs,1)
        i1 = pairs(ipair,1); i2 = pairs(ipair,2);
        [~,~,~,STATS] = ttest(comprel_jk(:,i1),comprel_jk(:,i2));
        t_jk = STATS.tstat/(ndata-1); % jackknife-corrected t
        PVAL = 2*tcdf(-abs(t_jk),ndata-1);
        bf10 = bf.ttest('T',t_jk,'N',ndata);
        if PVAL < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
        fprintf('\t%s vs %s: t_jk(%d) = %+.5f ; p = %.5f ; %s = %.5f\n',roilab{i1},roilab{i2},STATS.df,t_jk,PVAL,BFstr,BFval)
    end

    % horizontal bars (LPFC on top), error bars = jackknife s.e.m.
    iorder = nroi:-1:1;
    cfg = cfg_fig;
    cfg.figname = sprintf('comprel_ROI_%s',condnam{icond});
    cfg.type = 'horizontal';
    cfg.x = 100*comprel_jk(:,iorder);
    cfg.xqt = 100*[comprel_avg(iorder)-comprel_err(iorder);comprel_avg(iorder)+comprel_err(iorder)]';
    cfg.rgb = repmat(rgb(icond,:),nroi,1);
    cfg.labels = roilab(iorder);
    cfg.xlabel = '';
    cfg.ylabel = 'neural compression (%)'; % x-axis for horizontal bars
    cfg.ylim = [0 100];
    cfg.pbar = 1.6; cfg.HGT = 3;
    helper_plot_overall_bars(cfg);
    helper_ILL(gcf,cfg.pbar,cfg.HGT,3);
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end

%% CONLLR and NULLLR SOURCE: time course per ROI in PRE (supp)
if plotsource
    Atlas = cortex.Atlas;
    iatlas = find(contains({Atlas.Name},'Desikan-Killiany') & ~contains({Atlas.Name},'volume'));
    Scout = Atlas(iatlas).Scouts;
    roinam = {'superiorparietal','inferiortemporal','rostralmiddlefrontal'};
    roilab = {'SPL','ITC','LPFC'};
    icond = 1; % PRE
    for iroi = 1:numel(roinam)
        ivL = Scout(strcmp({Scout.Label},sprintf('%s L',roinam{iroi}))).Vertices;
        ivR = Scout(strcmp({Scout.Label},sprintf('%s R',roinam{iroi}))).Vertices;
        % participant x time sample
        yc = double(squeeze(mean(zsrc_conllr(ivL,:,:,icond),1)+mean(zsrc_conllr(ivR,:,:,icond),1))'/2);
        yn = double(squeeze(mean(zsrc_nulllr(ivL,:,:,icond),1)+mean(zsrc_nulllr(ivR,:,:,icond),1))'/2);
        cfg_stat = [];
        cfg_stat.xvec = xvec;
        y = {yc,yn};
        ynam = {'conllr','nulllr'};
        h = false(2,nsmp); % significant clusters of conllr and nulllr against 0
        for ivar = 1:2
            fprintf('%s %s: %s against 0, cluster-based permutation test\n',roilab{iroi},condnam{icond},ynam{ivar})
            stat = cluster_permutation(y{ivar},[],cfg_stat);
            for ic = 1:numel(stat.clusters)
                if stat.clusters(ic).p < 0.05
                    h(ivar,stat.clusters(ic).idx) = true;
                end
            end
        end
        cfg = cfg_fig;
        cfg.figname = sprintf('conllr_nulllr_source_%s_%s',roilab{iroi},condnam{icond});
        cfg.xvec = xvec;
        cfg.yavg = [mean(yc);mean(yn)]';
        cfg.yerr = [std(yc);std(yn)]'/sqrt(ndata);
        cfg.rgb = [rgb(icond,:);(1+rgb(icond,:))/2];
        cfg.ylim = [-0.025 0.1];
        cfg.ytick = 0:.05:.1;
        cfg.ylabel = 'coding precision (z)';
        cfg.pbar = 1.3;
        cfg.HGT = 4;
        cfg.saveornot = false; % saved below, with the bars of the significant clusters
        helper_plot_decoding(cfg)
        for ivar = 1:2 % conllr below nulllr, at the bottom of the plot
            c = bwconncomp(h(ivar,:));
            for ic = 1:c.NumObjects
                plot(xvec(c.PixelIdxList{ic}([1 end])),(cfg.ylim(1)+0.005*ivar)*[1 1],'Color',cfg.rgb(ivar,:),'LineWidth',2)
            end
        end
        if cfg_fig.saveornot
            print(fullfile(cfg_fig.figpath,cfg.figname),'-vector','-dpdf');
            fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg_fig.figpath)
        end
    end
end

%% CONLLR, NULLLR, COMPDIF: OPT/SUB time courses (compdif: main, conllr/nulllr: supp)
varnam = {'conllr','nulllr','compdif'};
zdec_all = {zdec_conllr_optsub,zdec_nulllr_optsub,zdec_compdif_optsub};
ylab = {'coding precision (z)','coding precision (z)','coding difference (z)'};
yl = [-0.05 0.4; -0.05 0.4; -0.2 0.4];
for ivar = 1:3
    zdec_avg_opt = zeros(nsmp,3);
    zdec_err_opt = zeros(nsmp,3);
    zdec_avg_sub = zeros(nsmp,3);
    zdec_err_sub = zeros(nsmp,3);
    for icond = 1:3
        x = cell2mat(zdec_all{ivar}(:,icond,1));
        zdec_avg_opt(:,icond) = mean(x,'omitnan');
        zdec_err_opt(:,icond) = std(x,'omitnan')/sqrt(ndata);
        x = cell2mat(zdec_all{ivar}(:,icond,2));
        zdec_avg_sub(:,icond) = mean(x,'omitnan');
        zdec_err_sub(:,icond) = std(x,'omitnan')/sqrt(ndata);
    end
    cfg = cfg_fig;
    cfg.figname = sprintf('%s_time_course_OPTSUB',varnam{ivar});
    cfg.xvec = xvec;
    cfg.yavg = zdec_avg_opt;
    cfg.yerr = zdec_err_opt;
    cfg.yavg_dashed = zdec_avg_sub;
    cfg.yerr_dashed = zdec_err_sub;
    cfg.ylim = yl(ivar,:);
    cfg.ytick = -.2:.2:.4;
    cfg.ylabel = ylab{ivar};
    cfg.rgb = rgb;
    cfg.pbar = 1.3;
    cfg.HGT = 4;
    helper_plot_decoding(cfg)
end

%% CONLLR, NULLLR, COMPDIF: OPT/SUB bars (conllr time window of each condition)
varnam = {'conllr','nulllr','compdif'};
zdec_all = {zdec_conllr_optsub,zdec_nulllr_optsub,zdec_compdif_optsub};
ylab = {'coding precision (z)','coding precision (z)','coding difference (z)'};
yl = [-0.05 0.4; -0.05 0.4; -0.2 0.4];
idxwin = idxwin_conllr_percond;
winmin = winmin_conllr_percond;
winmax = winmax_conllr_percond;
for ivar = 1:3
    X = nan(ndata,3,2); % participant x condition x opt/sub
    for icond = 1:3
        for iopt = 1:2
            x = cell2mat(zdec_all{ivar}(:,icond,iopt));
            X(:,icond,iopt) = mean(x(:,idxwin(icond,:)),2);
        end
    end
    for icond = 1:2
        cfg = cfg_fig;
        cfg.pbar = .6;
        cfg.HGT = 4;
        cfg.figname = sprintf('%s_optsub_mean_%s_%03d_%03d',varnam{ivar},condnam{icond},round(winmin(icond)*1000),round(winmax(icond)*1000));
        cfg.rgb = [rgb(icond,:);rgb(icond,:)];
        cfg.x = squeeze(X(:,icond,:));
        cfg.plotx = true;
        cfg.xlabel = condnam{icond};
        cfg.labels = {'correct' 'error'};
        cfg.ttest = true; cfg.tail = 'right'; % one-tailed paired t-test OPT > SUB, shown on the plot
        cfg.ylabel = ylab{ivar};
        cfg.ylim = yl(ivar,:); cfg.ytick = -.2:.2:.4;
        helper_plot_overall_bars(cfg);
        plot(xlim,mean(X(:,3,1))*[1 1],'Color',rgb(3,:),'LineWidth',1)
        plot(xlim,mean(X(:,3,2))*[1 1],'Color',rgb(3,:),'LineWidth',1,'LineStyle',':')
        helper_ILL(gcf,cfg.pbar,cfg.HGT,3);
        if cfg.saveornot
            print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
            fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
        end
    end
    fprintf('%s OPT/SUB (conllr time window of each condition):\n',upper(varnam{ivar}))
    % ttest OPT and SUB vs 0 (compression only)
    if strcmp(varnam{ivar},'compdif')
        optnam = {'OPT','SUB'};
        for icond = 1:3
            for iopt = 1:2
                [~,P,CI,STATS] = ttest(X(:,icond,iopt));
                bf10 = bf.ttest(X(:,icond,iopt));
                if P < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
                fprintf('\t%s %s vs 0: mean = %+.5f ; t%02d = %+.5f ; p = %.5f ; Cohen''s d = %+.5f ; CI = [%+.5f %+.5f] ; %s = %.5f\n', ...
                    condnam{icond},optnam{iopt},mean(X(:,icond,iopt)),STATS.df,STATS.tstat,P,STATS.tstat/sqrt(ndata),CI(1),CI(2),BFstr,BFval);
            end
        end
    end
    % ONE-TAILED TTEST: OPT>SUB
    for icond = 1:3
        [H,P,CI,STATS] = ttest(X(:,icond,1),X(:,icond,2),'Tail','right');
        fprintf('\t%s one-tailed ttest OPT>SUB ? : H = %d ; t = %+09.5f ; p = %5.5f ; CI = [%+5.5f %+5.5f]\n',condnam{icond},H,STATS.tstat,P,CI(1),CI(2));
    end
end

%% CONLLR, NULLLR, COMPDIF SOURCE (supp)
if plotsource
    varnam = {'conllr','nulllr','compdif'};
    clim = [0 0.1; 0 0.05; -0.1 .1];
    cmap = {cmap_royal_gramma(256),cmap_royal_gramma(256),cmap_ns_grey(256)}; % Brainstorm colormaps
    idxwin = idxwin_conllr_percond;
    winmin = winmin_conllr_percond;
    winmax = winmax_conllr_percond;
    for ivar = 3%1:3
        for icond = 1:2
            iwin = round([winmin(icond)*1000,winmax(icond)*1000]);
            switch varnam{ivar}
                case 'conllr'
                    x = zsrc_conllr(:,idxwin(icond,:),:,icond);
                case 'nulllr'
                    x = zsrc_nulllr(:,idxwin(icond,:),:,icond);
                case 'compdif'
                    x = zsrc_conllr(:,idxwin(icond,:),:,icond)-zsrc_nulllr(:,idxwin(icond,:),:,icond);
            end
            fprintf('max decoding %s cond %s : %4.4f\n',varnam{ivar},condnam{icond},max(mean(mean(x,3),2)))
            cfg = cfg_fig;
            cfg.cortex = cortex;
            cfg.map  = mean(mean(x,3),2); % mean over participants and time window
            cfg.clim = clim(ivar,:);
            cfg.cmap = cmap{ivar};
            cfg.views = {'top','left'};
            cfg.figname = sprintf('%s_time_avg_%s_%03d-%03d',varnam{ivar},condnam{icond},iwin(1),iwin(2));
            helper_plot_source(cfg)
        end
    end
end

%% CONLLR, NULLLR SOURCE OPT/SUB difference (supp)
if plotsource
    varnam = {'conllr','nulllr'};
    idxwin = idxwin_conllr_percond;
    winmin = winmin_conllr_percond;
    winmax = winmax_conllr_percond;
    newmin = -0.05;
    newmax = 0.05;
    for ivar = 1:2
        % template vertex x time sample x participant x condition x opt/sub
        tmp = load(fullfile(srcpath,sprintf('zdec_%s_particle_08Hz_optsub.mat',varnam{ivar})),'zdec');
        for icond = 1:3
            iwin = round([winmin(icond)*1000,winmax(icond)*1000]);
            x = tmp.zdec(:,idxwin(icond,:),:,icond,1)-tmp.zdec(:,idxwin(icond,:),:,icond,2);
            cfg = cfg_fig;
            cfg.cortex = cortex;
            cfg.map  = mean(mean(x,3),2); % mean over participants and time window
            cfg.clim = [newmin newmax];
            cfg.cmap = cmap_ns_grey(256);
            cfg.views = {'top','left'};
            cfg.figname = sprintf('%s_optsub_time_avg_%s_%03d-%03d',varnam{ivar},condnam{icond},iwin(1),iwin(2));
            helper_plot_source(cfg)
        end
        clearvars tmp
    end
end

if saveornot
    diary off
end
