%% analyses_supp_RT
% response time analysis (Supplementary Information)
% Supplementary Figure 2
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

% load data table
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = unique(dat.participant);
ndata = length(participants);
ncond = 3;
condnam = {'PRE' 'POST' 'REF'};
seqlen0 = [4,8,12];
nseqlen = length(seqlen0);

%% check overall performance to determine if participant should remain in
% the analysis (not more than 3*std away from mean)
nstd = 3;
pcor = nan(ndata,1);
for idata = 1:ndata
    ifilt = dat.participant == participants(idata);% & dat.ncard == 1;
    r = dat.response(ifilt);
    c = dat.seqllr(ifilt);
    pcor(idata) = mean(sign(1.5-r) == sign(c),'omitnan');
end
pcor_loo_avg = nan(ndata,1);
pcor_loo_std = nan(ndata,1);
for idata = 1:ndata
    pcor_loo_avg(idata) = mean(pcor(setdiff(1:ndata,idata)));
    pcor_loo_std(idata) = std(pcor(setdiff(1:ndata,idata)));
end
I = find((pcor<(pcor_loo_avg+pcor_loo_std*(-nstd))));
if ~isempty(I)
    for i = 1:numel(I)
        fprintf('remove participant S%02d (idata %02d/%02d) from analysis.\n',participants(I(i)),I(i),ndata)
    end
else
    fprintf('keep all participant for analysis (N=%02d)\n',ndata)
end
participants_all = participants;
participants = setdiff(participants,participants(I));
ndata = numel(participants);
[~,idata2keep] = ismember(participants,participants_all);

% compute accuracy per condition per seqlen
pcor_seqlen = nan(ndata,ncond,nseqlen);
for idata = 1:ndata
    for icond = 1:ncond
        for iseqlen = 1:3
            ifilt = dat.participant == participants(idata) & dat.condtn == icond & dat.seqlen == seqlen0(iseqlen);
            resp = dat.response(ifilt);
            seqcat = dat.seqcat(ifilt);
            pcor_seqlen(idata,icond,iseqlen) = mean(resp == seqcat);
        end
    end
end

for icond = 1:ncond
    fprintf('COND %s: accuracy (4, 8, 12) = %s\n',condnam{icond},mat2str(squeeze(mean(pcor_seqlen(:,icond,:),1))',3));
    anova_1way_ranova(squeeze(pcor_seqlen(:,icond,:)), 'seqlen');
end

%% Compare conditions (ANOVA+ttest)
RT = nan(ndata,ncond);
RT_std = nan(ndata,ncond);
for idata = 1:ndata
    for icond = 1:ncond
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        RT(idata,icond) = mean(dat.RT(ifilt));
        RT_std(idata,icond) = std(dat.RT(ifilt));
    end
end

% plot RT
X = RT;
cfg = cfg_fig;
cfg.figname = 'RT_comparison_condition_bars';
cfg.rgb = rgb([2 1],:);
cfg.x = X(:,[2 1]);
cfg.plotx = true;
cfg.wdth = 0.5;
cfg.HGT = 4.5;
cfg.pbar = 0.5;
cfg.xlabel = 'condition';
cfg.labels = condnam([2 1]);
cfg.ylabel = 'Response time (s)';
cfg.ylim = [0 2.5];cfg.ytick = 0:0.5:2.5;
helper_plot_overall_bars(cfg);
plot(xlim,mean(X(:,3))*[1 1],'Color',rgb(3,:),'LineWidth',1) % add REF
helper_ILL(gcf,cfg.pbar,cfg.HGT,false);
% re-save
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end

% analyses on log(RT)
X = log(RT);
fprintf('<strong>LOG(Response time) - compare conditions - ANOVA+ post-hoc TTEST</strong>\n')
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
fprintf('\n')

% CORRELATION PRE-POST
fprintf('<strong>Response time - PRE versus POST - CORR (Pearson, on log(RT))</strong>\n')
cfg = cfg_fig;
cfg.figname = 'RT_CORR_PRE_POST';
cfg.models = {'POST' 'PRE'};
cfg.HGT = 4.5;
cfg.xavg = RT(:,2);
cfg.yavg = RT(:,1);
cfg.xstd = RT_std(:,2);cfg.ystd = RT_std(:,1);
cfg.xlabel = 'RT (s) POST';cfg.ylabel = 'RT (s) PRE';
cfg.ylim = [0 2.5];cfg.xlim = [0 2.5];
cfg.xtick = 0:.5:2.5;
cfg.verbose = false;
[R,P,RLO,RUP]=corrcoef(log(cfg.xavg),log(cfg.yavg));
fprintf('\tRHO = %.5f ; PVAL = %.5f ; RHO^2 = %.5f ; CI = [%+.5f %+.5f]\n',R(1,2),P(1,2),R(1,2)^2,RLO(1,2),RUP(1,2));
helper_plot_corr(cfg)
fprintf('\n')

%% RT per seqlen and condition
fprintf('<strong>LOG(Response time) - inside each condition per sequence length</strong>\n')
RT_seq = nan(ndata,ncond,nseqlen);
for icond = 1:ncond
    for idata = 1:ndata
        for iseqlen = 1:nseqlen
            ifilt = dat.participant == participants(idata) & dat.condtn == icond & dat.seqlen == seqlen0(iseqlen);
            RT_seq(idata,icond,iseqlen) = mean(dat.RT(ifilt));
        end
    end
end

pbar = .8; HGT = 4.5;
figname = 'RT_per_seqlen';
figure('Color','white','Name',figname);
for icond = 1:2%3
    X = squeeze(RT_seq(:,icond,:));
    nstim = size(X,2);
    
    xpos = [1 2 3];
    xoff = .4;
    
    rgb_mod = rgb(icond,:);
    xlim([min(xpos)-xoff,max(xpos)+xoff]);
    hold on
    line(xlim,[0 0],'Color','k')
    for istim =1:nstim
        x = X(:,istim);
        pos = xpos(istim);
        wid = 0.2;
        xk = linspace(min(x),max(x),100);
        [pk,~,u] = ksdensity(x,xk);
        while true
            z = findpeaks(pk,xk);
            if numel(z) <= 1
                break
            else
                u = u*1.1;
                pk = ksdensity(x,xk,'BandWidth',u);
            end
        end
        pk = pk/max(pk);
        str = interp1(xk,pk,x);
        xmed = interp1(xk,pk,median(x));
        jit = linspace(-.75,.75,numel(x));
        for i = 1:numel(x)
            plot(pos+jit(i)*str(i)*wid,x(i),'o','MarkerFaceColor',(1+rgb_mod)/2,'MarkerEdgeColor','no','MarkerSize',2);
        end
        plot(pos*[1,1],mean(x)+std(x)*[-1,+1],'k-','LineWidth',0.75);
        plot(pos,mean(x),'ko','MarkerSize',5,'MarkerFaceColor',rgb_mod,'LineWidth',0.75);

    end
end
line(xpos,mean(squeeze(RT_seq(:,3,:))),'Color',rgb(3,:))
set(gca,'XTick',xpos,'XTickLabel',{'4' '8' '12'},'YLim',[0 2.5]);%,'YTick',-1:0.5:1);
xlabel('sequence length')
ylabel('Response time (s)')
helper_ILL(gcf,pbar,HGT,false);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end

anova_2way_ranova(log(RT_seq), {'condition','seqlen'});
for icond = 1:3
    fprintf('COND %s\n',condnam{icond})
    [stat,rm] = anova_1way_ranova(squeeze(log(RT_seq(:,icond,:))), 'seqlen', '', true); 
end
fprintf('\n')

%% GLMFIT to predict RT
B_LLR_seqlen = nan(ndata,ncond,3);
B_LLR_TRUE_NULL = nan(ndata,ncond,3);
B_LLR_TRUE_NULL_seqlen = nan(ndata,ncond,4);
kappa = 0.5;

for idata = 1:ndata
    for icond = 1:3
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        RT_tmp = log(dat.RT(ifilt));
        resp = 1.*(dat.response(ifilt) == 1) -1.*(dat.response(ifilt) == 2);
        seqllr = dat.seqllr(ifilt);
        seqlen = dat.seqlen(ifilt);
        seqang = dat.seqang(ifilt);
        catang = dat.catang(ifilt);
        
        absllr_true = abs(cellfun(@(x,y) sum(2*kappa*sin(2*(x-y))),seqang,num2cell(catang)));
        absllr_null = abs(cellfun(@(x,y) sum(2*kappa*sin(2*(x-y))),seqang,num2cell(catang+pi/4)));
        
        % predict RT:
        Y = RT_tmp;
        
        X = [seqllr.*resp seqlen];
        B_LLR_seqlen(idata,icond,:) = glmfit(zscore(X),Y);% z-score to have comparable reg coefficients

        X = [absllr_true absllr_null];
        B_LLR_TRUE_NULL(idata,icond,:) = glmfit(X,Y);% NO NEED TO Z-SCORE

        X = [absllr_true absllr_null seqlen];
        B_LLR_TRUE_NULL_seqlen(idata,icond,:) = glmfit(zscore(X),Y);% z-score to have comparable reg coefficients
    end
end


%% Compare regressors LLR_SIGNED SEQLEN
pbar = .5; HGT = 4;
coeffnam = {'intercept','signed LLR' 'SEQLEN'};
X2plot = B_LLR_seqlen;
for icoeff = 2:size(X2plot,3) % skip intercept
    X = X2plot(:,:,icoeff);
    fprintf('<strong>Response time - predicted from signed LLR and SEQLEN - glmfit coeff %s</strong>\n',coeffnam{icoeff})
    [H,p,ci,stats]=ttest(X); effd = stats.tstat/sqrt(size(X,1));
    for icond = 1:ncond
        bf10 = bf.ttest(X(:,icond));
        if H(icond)==0, BFval = 1/bf10; BFstr = 'BF01';
        elseif H(icond) == 1, BFval = bf10; BFstr = 'BF10';
        end
        fprintf('\tttest against 0 (%s): t = %+09.5f ; p = %5.5f ; Cohen''s d = %+09.5f ; CI = [%+5.5f %+5.5f]; %s = %5.5f\n',condnam{icond},stats.tstat(icond),p(icond),effd(icond),ci(1,icond),ci(2,icond),BFstr,BFval)
    end
    anova_1way_ranova(X, 'condition', '', true);

    figname = sprintf('RT_glmfit_LLR_signed_seqlen_coeff%d',icoeff);
    figure('Color','white','Name',figname);
    xpos = [1 2];
    xoff = .5;
    xlim([min(xpos)-xoff,max(xpos)+xoff]);
    hold on
    line(xlim,[0 0],'Color','k','LineWidth',0.75)

    ordercond = [2 1];
    for icond=1:2
        x = X(:,ordercond(icond));
        rgb_mod = rgb(ordercond(icond),:);
        pos = xpos(icond);
        wid = 0.2;
        xk = linspace(min(x),max(x),100);
        [pk,~,u] = ksdensity(x,xk);
        while true
            z = findpeaks(pk,xk);
            if numel(z) <= 1
                break
            else
                u = u*1.1;
                pk = ksdensity(x,xk,'BandWidth',u);
            end
        end
        pk = pk/max(pk);
        str = interp1(xk,pk,x);
        xmed = interp1(xk,pk,median(x));
        jit = linspace(-.75,+.75,numel(x));
        for i = 1:numel(x)
            plot(pos+jit(i)*str(i)*wid,x(i),'o','MarkerFaceColor',(1+rgb_mod)/2,'MarkerEdgeColor','no','MarkerSize',2);
        end
        plot(pos*[1,1],mean(x)+std(x)*[-1,+1],'k-','LineWidth',0.75);
        plot(pos,mean(x),'ko','MarkerSize',5,'MarkerFaceColor',rgb_mod,'LineWidth',0.75);
    end
    if icoeff == 1
        yl = [-1.2 1.2];yt = -1:.5:1;
    elseif icoeff == 2
        yl = [-.24 .24];yt = -.2:0.2:.2;
    elseif icoeff == 3
        yl = [-.12 .12];yt = -.1:0.1:.1;
    end
    plot(xlim,mean(X2plot(:,3,icoeff))*[1 1],'Color',rgb(3,:),'LineWidth',1.5)
    set(gca,'XTick',xpos,'XTickLabel',condnam(ordercond),'YLim',yl,'YTick',yt);
    xlabel('condition')
    ylabel(sprintf('%s',coeffnam{icoeff}))
    helper_ILL(gcf,pbar,HGT,false);
    if cfg_fig.saveornot
        print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
    end
    fprintf('\n')
end


%% Compare regressors TRUE and NULL projection axis
pbar = .5; HGT = 4;
coeffnam = {'intercept','absLLR(TRUE)' 'absLLR(NULL)'};
X2plot = B_LLR_TRUE_NULL;

for icoeff = 2:size(X2plot,3) % skip intercept
    X = X2plot(:,:,icoeff);
    fprintf('<strong>Response time - predicted from ABSLLR TRUE AND NULL - glmfit coeff %s</strong>\n',coeffnam{icoeff})
    [H,p,ci,stats]=ttest(X); effd = stats.tstat/sqrt(size(X,1));
    for icond = 1:ncond
        bf10 = bf.ttest(X(:,icond));
        if H(icond)==0, BFval = 1/bf10; BFstr = 'BF01';
        elseif H(icond) == 1, BFval = bf10; BFstr = 'BF10';
        end
        fprintf('\tttest against 0 (%s): t = %+09.5f ; p = %5.5f ; Cohen''s d = %+09.5f ; CI = [%+5.5f %+5.5f]; %s = %5.5f\n',condnam{icond},stats.tstat(icond),p(icond),effd(icond),ci(1,icond),ci(2,icond),BFstr,BFval)
    end
    [stat,rm] = anova_1way_ranova(X, 'condition');
    fprintf('post-hoc ttest between conditions:\n')
    c = multcompare(rm,'condition'); % correct for multiple comparisons
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
    fprintf('\n')

    figname = sprintf('RT_glmfit_ABSLLR_TRUE_NULL_coeff%d',icoeff);
    figure('Color','white','Name',figname);
    xpos = [1 2];
    xoff = .5;
    xlim([min(xpos)-xoff,max(xpos)+xoff]);
    hold on
    line(xlim,[0 0],'Color','k','LineWidth',0.75)
    
    ordercond = [2 1];
    for icond=1:2
        x = X2plot(:,ordercond(icond),icoeff);
        rgb_mod = rgb(ordercond(icond),:);
        pos = xpos(icond);
        wid = 0.2;
        xk = linspace(min(x),max(x),100);
        [pk,~,u] = ksdensity(x,xk);
        while true
            z = findpeaks(pk,xk);
            if numel(z) <= 1
                break
            else
                u = u*1.1;
                pk = ksdensity(x,xk,'BandWidth',u);
            end
        end
        pk = pk/max(pk);
        str = interp1(xk,pk,x);
        xmed = interp1(xk,pk,median(x));
        jit = linspace(-.75,+.75,numel(x));
        for i = 1:numel(x)
            plot(pos+jit(i)*str(i)*wid,x(i),'o','MarkerFaceColor',(1+rgb_mod)/2,'MarkerEdgeColor','no','MarkerSize',2);
        end
        plot(pos*[1,1],mean(x)+std(x)*[-1,+1],'k-','LineWidth',0.75);
        plot(pos,mean(x),'ko','MarkerSize',5,'MarkerFaceColor',rgb_mod,'LineWidth',0.75);
    end
    if icoeff == 1
        yl = [-1.2 1.2];yt = -1:.5:1;
    else
        yl = [-.12 .12];yt = -.1:0.1:.1;
    end
    plot(xlim,mean(X2plot(:,3,icoeff))*[1 1],'Color',rgb(3,:),'LineWidth',1.5)
    set(gca,'XTick',xpos,'XTickLabel',condnam(ordercond),'YLim',yl,'YTick',yt);
    xlabel('condition')
    ylabel(sprintf('%s',coeffnam{icoeff}))
    helper_ILL(gcf,pbar,HGT,false);
    if cfg_fig.saveornot
        print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
    end
end



%% Compare regressors TRUE and NULL projection axis (cotrolling for seqlen)
pbar = .5; HGT = 4;
coeffnam = {'intercept','absLLR(TRUE)' 'absLLR(NULL)' 'seqlen'};
X2plot = B_LLR_TRUE_NULL_seqlen;

for icoeff = 2:size(X2plot,3) % skip intercept
    X = X2plot(:,:,icoeff);
    fprintf('<strong>Response time - predicted from ABSLLR TRUE AND NULL AND SEQLEN - glmfit coeff %s</strong>\n',coeffnam{icoeff})
    [H,p,ci,stats]=ttest(X); effd = stats.tstat/sqrt(size(X,1));
    for icond = 1:ncond
        bf10 = bf.ttest(X(:,icond));
        if H(icond)==0, BFval = 1/bf10; BFstr = 'BF01';
        elseif H(icond) == 1, BFval = bf10; BFstr = 'BF10';
        end
        fprintf('\tttest against 0 (%s): t = %+09.5f ; p = %5.5f ; Cohen''s d = %+09.5f ; CI = [%+5.5f %+5.5f]; %s = %5.5f\n',condnam{icond},stats.tstat(icond),p(icond),effd(icond),ci(1,icond),ci(2,icond),BFstr,BFval)
    end
    [stat,rm] = anova_1way_ranova(X, 'condition');
    fprintf('post-hoc ttest between conditions:\n')
    c = multcompare(rm,'condition'); % correct for multiple comparisons
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
    fprintf('\n')

    figname = sprintf('RT_glmfit_ABSLLR_TRUE_NULL_seqlen_coeff%d',icoeff);
    figure('Color','white','Name',figname);
    xpos = [1 2];
    xoff = .5;
    xlim([min(xpos)-xoff,max(xpos)+xoff]);
    hold on
    line(xlim,[0 0],'Color','k','LineWidth',0.75)
    
    ordercond = [2 1];
    for icond=1:2
        x = X2plot(:,ordercond(icond),icoeff);
        rgb_mod = rgb(ordercond(icond),:);
        pos = xpos(icond);
        wid = 0.2;
        xk = linspace(min(x),max(x),100);
        [pk,~,u] = ksdensity(x,xk);
        while true
            z = findpeaks(pk,xk);
            if numel(z) <= 1
                break
            else
                u = u*1.1;
                pk = ksdensity(x,xk,'BandWidth',u);
            end
        end
        pk = pk/max(pk);
        str = interp1(xk,pk,x);
        xmed = interp1(xk,pk,median(x));
        jit = linspace(-.75,+.75,numel(x));
        for i = 1:numel(x)
            plot(pos+jit(i)*str(i)*wid,x(i),'o','MarkerFaceColor',(1+rgb_mod)/2,'MarkerEdgeColor','no','MarkerSize',2);
        end
        plot(pos*[1,1],mean(x)+std(x)*[-1,+1],'k-','LineWidth',0.75);
        plot(pos,mean(x),'ko','MarkerSize',5,'MarkerFaceColor',rgb_mod,'LineWidth',0.75);
    end
    if icoeff == 1
        yl = [-1.2 1.2];yt = -1:.5:1;
    else
        yl = [-.12 .12];yt = -.1:0.1:.1;
    end
    plot(xlim,mean(X2plot(:,3,icoeff))*[1 1],'Color',rgb(3,:),'LineWidth',1.5)
    set(gca,'XTick',xpos,'XTickLabel',condnam(ordercond),'YLim',yl,'YTick',yt);
    xlabel('condition')
    ylabel(sprintf('%s',coeffnam{icoeff}))
    helper_ILL(gcf,pbar,HGT,false);
    if cfg_fig.saveornot
        print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
    end
end