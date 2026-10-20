%% analyses_figure_2_behavior
% this script runs behavioral analysis corresponding to Figure 2
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

ordercond = [2 1 3];

% load data table organized per sequence
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = unique(dat.participant);
ndata = length(participants);

condnam = {'PRE' 'POST' 'REF'};
ncond = numel(condnam);
seqlen0 = [4 8 12];
nseqlen = numel(seqlen0);

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

%% perf per condition + STATS
fprintf('<strong>Accuracy (true) - compare conditions - ANOVA + TTEST </strong>\n')
perf = zeros(ndata,ncond);
perf_std = zeros(ndata,ncond);
for icond = 1:ncond
    for idata = 1:ndata
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        resp = dat.response(ifilt);
        seqcat = dat.seqcat(ifilt);
        perf(idata,icond) = mean(resp == seqcat);
        perf_std(idata,icond) = sqrt(perf(idata,icond)*(1-perf(idata,icond))/sum(ifilt));
    end
end

fprintf('Human performance:\n')
for icond = 1:ncond
    fprintf('\t%s:\t%4.1f ±%4.1f %%\n',condnam{icond},mean(perf(:,icond))*100,100*std(perf(:,icond))/sqrt(ndata))
end

X = perf;
[~,rm] = anova_1way_ranova(X, 'condition'); 
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest:\n')
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
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+07.3f ; p(corrected) = %4.4f ; Cohen''s d = %+04.4f ; CI(corrected) = [%+4.4f %+4.4f] ; %s = %4.4f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end
fprintf('\n')


%% compute Bayes optimal performance
fprintf('<strong>Accuracy (true) - compare to Bayes-optimal - TTEST </strong>\n')
perf_bayes = nan(ndata,ncond);
for icond = 1:ncond
    for idata = 1:ndata
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        seqllr = dat.seqllr(ifilt);
        seqcat = dat.seqcat(ifilt);
        perf_bayes(idata,icond) = mean(sign(seqllr) == sign(3-2*seqcat));
    end
    [H,P,CI,STATS] = ttest(perf(:,icond)-perf_bayes(:,1)); effd = STATS.tstat/sqrt(ndata);
    fprintf('\tcond %s vs optimal \t: H = %d ; t%02d = %+09.5f ; p = %5.5f ; Cohen''s d = %+09.5f ; CI = [%+5.5f %+5.5f]\n',condnam{icond},H,STATS.df,STATS.tstat,P,effd,CI(1),CI(2));
end
% actually no need to differentiate the conditions
perf_bayes = mean(perf_bayes(:,1));

%% Plot overall accuracy
% TRUE
cfg = cfg_fig;
cfg.pbar = .5;
cfg.wdth = 0.5;
cfg.HGT = HGT;
cfg.figname = 'performance_bars';
cfg.rgb = rgb([2 1],:);
cfg.x = perf(:,[2 1]);
cfg.plotx = true;
cfg.ttest = true;
cfg.xlabel = 'condition';
cfg.labels = condnam([2 1]);
cfg.ylabel = 'accuracy';
cfg.ylim = [0.5 1];cfg.ytick = .5:.1:1;
helper_plot_overall_bars(cfg);
plot(xlim,mean(perf(:,3))*[1 1],'Color',rgb(3,:),'LineWidth',1) % add REF
plot(xlim,perf_bayes*[1 1],'k:') % add Bayes optimal
% re-save
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end

%% Sensitivity to Bayes-optimal evidence (binning of seqllr)
fprintf('<strong>Sensitivity to Bayes-optimal evidence</strong>\n')
binedges = [-Inf -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 Inf];
nbin = numel(binedges)-1;
% perf per condition + STATS
perf_bin = nan(ndata,nbin,ncond);
beta = nan(ndata,ncond,2);

xvec_sim = -4:.01:+4;
y_sim = nan(ndata,numel(xvec_sim),3);

for icond = 1:ncond
    for idata = 1:ndata
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        resp = dat.response(ifilt);
        seqcat = dat.seqcat(ifilt);
        seqllr = dat.seqllr(ifilt);
        beta(idata,icond,:)= glmfit(seqllr, resp==1, 'binomial', 'link', 'logit');

        y_sim(idata,:,icond) = glmval(squeeze(beta(idata,icond,:)), xvec_sim', 'logit');
        binidx = discretize(seqllr,binedges);
        bins = unique(binidx);
        for ibin = 1:numel(bins)
            perf_bin(idata,bins(ibin),icond) = mean(resp(binidx==bins(ibin)) == 1);
        end
    end
end

% plot sigmoid curve
figure('Color','white','Name','sensitivity_evidence_curve'); hold on
pbar = 1.3;
nbin = size(perf_bin,2);
xvec = -4:1:+4;

Yavg = squeeze(mean(perf_bin));
Yerr = squeeze(std(perf_bin,[],1))/sqrt(ndata);

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
    print(fullfile(cfg.figpath,fig.Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',fig.Name,cfg.figpath)
end

% plot sensitivity coefficient bars
cfg = cfg_fig;
cfg.pbar = .5;
cfg.wdth = 0.5;
cfg.HGT  = HGT;
cfg.figname = 'sensitivity_evidence_bars';
cfg.rgb = rgb([2 1],:);
cfg.x = beta(:,[2 1],2);
cfg.plotx = true;
cfg.ttest = true;
cfg.xlabel = 'condition';
cfg.labels = condnam([2 1]);
cfg.ylabel = 'sensitivity coefficient';
cfg.ylim = [0 2.5];cfg.ytick = 0:1:2;cfg.yticklabel={'0','1.0','2.0'};
helper_plot_overall_bars(cfg);
plot(xlim,mean(beta(:,3,2))*[1,1],'Color',rgb(3,:),'LineWidth',1);
if cfg.saveornot
    print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
end

X = beta(:,:,2);
[~,rm] = anova_1way_ranova(X, 'condition'); 
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest:\n')
for icond = 1:ncond
    jcond = mod(icond,ncond)+1;
    [~,~,~,STATS] = ttest(X(:,icond),X(:,jcond));
    effd = STATS.tstat/sqrt(size(X,1));
    % identify row corresponding to post-hoc ttest in the multcompare
    % output
    irow = double(c.condition_1) == icond & double(c.condition_2) == jcond;
    P  = c.pValue(irow);
    CI = [c.Lower(irow) c.Upper(irow)];
    H  = P < 0.05;
    bf10 = bf.ttest(X(:,icond),X(:,jcond));
    if H==0, BFval=1/bf10; BFstr='BF01';else,BFval=bf10;BFstr='BF10';end
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+06.2f ; p(corrected) = %4.4f ; Cohen''s d = %+04.4f ; CI(corrected) = [%+4.4f %+4.4f] ; %s = %4.4f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end
fprintf('\n')

%% Choice variability
fprintf('<strong>Choice variability</strong>\n')
B = nan(ndata,nseqlen,ncond);
for idata = 1:ndata
    fprintf('Processing participant %02d/%02d ...',idata,ndata)
    for icond = 1:ncond
        for iseqlen = 1:nseqlen
            ifilt = dat.participant == participants(idata) & dat.condtn == icond & dat.seqlen == seqlen0(iseqlen);
            X = dat.seqllr(ifilt);
            Y = dat.response(ifilt) == 1;
            B_tmp = glmfit(X,Y,'binomial','link','probit');
            B(idata,iseqlen,icond) = B_tmp(2);
        end
    end
    fprintf('done.\n')
end
sigma = 1./B;

% load model with siginf only
if ~exist('./results/results_fit_choice_variability.mat','file')
    % run git
    out_fit = cell(numel(participants_all),1);
    parfor idata = 1:numel(participants_all)
        out_fit{idata} = run_fit_choice_variability(idata);
    end
    save('./results/results_fit_choice_variability.mat','out_fit')
else
    load('./results/results_fit_choice_variability.mat','out_fit')
end
out_fit = out_fit(idata2keep);
siginf1 = cellfun(@(s) s.siginf1,out_fit);
siginf2 = cellfun(@(s) s.siginf2,out_fit);
siginf3 = cellfun(@(s) s.siginf3,out_fit);
siginf = cat(2,siginf1,siginf2,siginf3);
sigma_model = nan(ndata,3,3);
for iseqlen = 1:numel(seqlen0)
    for icond = 1:ncond
        sigma_model(:,iseqlen,icond) = sqrt(seqlen0(iseqlen)) * siginf(:,icond);
    end
end
%
figure('Color','white','Name','choice_variability_seqlen')
pbar = .8;
xpos = [1 2 3];
xoff = .4;

xlim([min(xpos)-xoff,max(xpos)+xoff]);
xvec = xpos;

for icond = 1:3
    rgb_cond = [(1+(1+rgb(icond,:))/2)/2;...
        (1+rgb(icond,:))/2;...
        rgb(icond,:)];
    xvec = xpos; 
    ntrl = numel(xvec);
    Yavg = mean(sigma(:,:,icond).^2); Yerr = std(sigma(:,:,icond).^2)/sqrt(ndata);
    yavg = mean(sigma_model(:,:,icond).^2);
    yerr = std(sigma_model(:,:,icond).^2)/sqrt(ndata);
    clearvars yerr_vec
    yerr_vec(1,:) = reshape(yavg,1,ntrl) + reshape(yerr,1,ntrl);
    yerr_vec(2,:) = reshape(yavg,1,ntrl) - reshape(yerr,1,ntrl);
    hold on
    patch([xvec,fliplr(xvec)],[yerr_vec(1,:),fliplr(yerr_vec(2,:))], ...
        (1+rgb(icond,:))/2,'EdgeColor','none','FaceAlpha',.5);
    plot(xvec,yavg,'Color',rgb(icond,:),'LineWidth',1.5);
    for i = 1:numel(seqlen0)
        plot(xvec(i)*[1,1],Yavg(i)+Yerr(i)*[+1,-1],'k-','Linewidth',0.75);
        plot(xvec(i),Yavg(i),'ko','MarkerSize',10,'MarkerFaceColor',rgb_cond(i,:),'Linewidth',0.75);
    end
end
ylim([-1 15])
set(gca,'Layer','top','Box','off','TickDir','out','PlotBoxAspectRatio',[pbar,1,1]);
set(gca,'XTick',xvec,'XTickLabel',{'4' '8' '12'},'YTick',0:5:15);
xlabel('sequence length');
ylabel('choice variability');
helper_ILL(gcf,pbar,HGT);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end
fprintf('\n')

anova_2way_ranova(sigma.^2, {'seqlen','condition'});


%% integration kernel
fprintf('<strong>Integration kernel</strong>\n')
X_l = zeros(size(dat,1),12);
X_r = zeros(size(dat,1),12);
for iseq = 1:size(dat,1)
    X_l(iseq,1:dat.seqlen(iseq)) = dat.LLR{iseq};
    X_r(iseq,(12-dat.seqlen(iseq)+1):12) = dat.LLR{iseq};
end

B_l = nan(ndata,5,ncond);
B_r = nan(ndata,5,ncond);
for idata = 1:ndata
    fprintf('Processing participant %02d/%02d',idata,ndata)

    for icond = 1:ncond
        fprintf('.')
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;

        Y = dat.response(ifilt,:)==1;
        % right aligned
        X = X_r(ifilt,:);
        X = [sum(X(:,1:8),2) X(:,9:end)];
        B_tmp  = glmfit(X,Y,'binomial','link','probit');
        B_r(idata,:,icond) = B_tmp(2:end)';
        % left aligned
        X = X_l(ifilt,:);
        X = [X(:,1:4) sum(X(:,5:end),2)];
        B_tmp  = glmfit(X,Y,'binomial','link','probit');
        B_l(idata,:,icond) = B_tmp(2:end)';
    end
    fprintf('done.\n')
end

pbar = 1.3;
figure('Color','white','Name','integration_kernel');
for icond = 1:3
    hold on
    % create concatenated vector
    yvec = nan(size(B_l,1),12);
    yvec(:,1:4) = B_l(:,1:4,icond);
    yvec(:,9:12) = B_r(:,2:end,icond);
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
    print(fullfile(cfg_fig.figpath,fig.Name),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n\n',fig.Name,cfg_fig.figpath)
end
anova_2way_ranova(cat(2,B_l(:,1:4,:),B_r(:,2:end,:)), {'seqpos','condition'}); 


%% Bias-variance decomposition
% load data table all conditions and repetitions aligned
% per sequence
load(fullfile(cfg_path.data,'data_table_cond.mat'),'dat');

% COMPUTE FRACTION CONSISTENT for repeated presentation of sequence for each condition
fprintf('<strong>Fraction consistent responses on repeated trials depending on condition</strong>\n')
prop_same = nan(ndata,ncond);
for idata = 1:ndata
    ifilt = dat.participant == participants(idata);
    resp = dat.response(ifilt,:,:); % iseq x icond x irepetition
    prop_same(idata,:) = mean(resp(:,:,1) == resp(:,:,2));
end
for icond = 1:3
    fprintf('cond%d: %4.3f ± %4.3f\n',icond,mean(prop_same(:,icond))*100,std(prop_same(:,icond))/sqrt(ndata)*100)
end

% plot fraction consistent
cfg = cfg_fig;
cfg.pbar = .5;
cfg.HGT = HGT;
cfg.wdth = 0.5;
cfg.figname = 'fraction_consistent_repeated_trials';
cfg.rgb = rgb(ordercond(1:2),:);
cfg.x = prop_same(:,ordercond(1:2));
cfg.plotx = true;
cfg.ttest = true;
cfg.xlabel = 'condition';
cfg.labels = condnam(ordercond(1:2));
cfg.ylabel = 'fraction consistent';
cfg.ylim = [0.5 1];cfg.ytick = 0.5:.1:1;
helper_plot_overall_bars(cfg);
plot(xlim,mean(prop_same(:,3))*[1 1],'Color',rgb(3,:),'Linewidth',1);
helper_ILL(gcf,cfg.pbar,cfg.HGT,false);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end

X = prop_same;
[~,rm] = anova_1way_ranova(X, 'condition'); 
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest:\n')
for icond = 1:ncond
    jcond = mod(icond,ncond)+1;
    [~,~,~,STATS] = ttest(X(:,icond),X(:,jcond));
    effd = STATS.tstat/sqrt(size(X,1));
    % identify row corresponding to post-hoc ttest in the multcompare
    % output
    irow = double(c.condition_1) == icond & double(c.condition_2) == jcond;
    P  = c.pValue(irow);
    CI = [c.Lower(irow) c.Upper(irow)];
    H  = P < 0.05;
    bf10 = bf.ttest(X(:,icond),X(:,jcond));
    if H==0, BFval=1/bf10; BFstr='BF01';else,BFval=bf10;BFstr='BF10';end
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+06.2f ; p(corrected) = %4.4f ; Cohen''s d = %+04.4f ; CI(corrected) = [%+4.4f %+4.4f] ; %s = %4.4f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end
fprintf('\n')

if ~exist('./results/results_sim_rsq.mat','file')
    out_sim = run_bias_variance;
else
    load('./results/results_sim_rsq.mat','out_sim')
end
out_sim = out_sim(idata2keep,:);

siginf = cell2mat(cellfun(@(x)x.siginf,out_sim(:,1),'UniformOutput',0));
rsq    = cellfun(@(x)x.rsq,out_sim(1,:));
nrsq   = length(rsq);
pcon_rsq = nan(ndata,nrsq,ncond);
for icond = 1:ncond
    pcon_rsq(:,:,icond) = cellfun(@(x) x(icond),cellfun(@(x)x.pcon_rsq,out_sim,'UniformOutput',0));
end

% find participants' r_squared depending from its fraction consistent for each condition
rsq_dat = zeros(ndata,ncond);
for idata = 1:ndata
    for icond = 1:ncond
        if prop_same(idata,icond) < pcon_rsq(idata,1,icond)
            continue % no interpolation possible
        else
            rsq_dat(idata,icond) = interp1(pcon_rsq(idata,:,icond),rsq,prop_same(idata,icond),'spline');
            if rsq_dat(idata,icond)<0
                rsq_dat(idata,icond) = 0;
            end
        end
    end
end


% bias variance decomposition
bv = cat(3,sqrt((siginf.^2).*(rsq_dat)),sqrt((siginf.^2).*(1-rsq_dat)));

% bias variance decomp - bias
fprintf('<strong>bias-variance decomposition - BIAS - ANOVA</strong>\n')
cfg = cfg_fig;
cfg.pbar = .5;
cfg.HGT = HGT;
cfg.wdth = 0.5;
cfg.figname = 'bias_variance_siginf_bias';
cfg.rgb = rgb(ordercond(1:2),:);
cfg.x = bv(:,ordercond(1:2),1);
cfg.plotx = true;
cfg.ttest = true;
cfg.xlabel = 'condition';
cfg.labels = condnam(ordercond(1:2));
cfg.ylabel = 'bias term (s.d.)';
cfg.ylim = [0 2];cfg.ytick = 0:.5:2;
helper_plot_overall_bars(cfg);
plot(xlim,mean(bv(:,3,1))*[1 1],'Color',rgb(3,:),'Linewidth',1);
helper_ILL(gcf,cfg.pbar,cfg.HGT,false);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end


X = bv(:,:,1);
[~,rm] = anova_1way_ranova(X, 'condition'); 
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest:\n')
for icond = 1:ncond
    jcond = mod(icond,ncond)+1;
    [~,~,~,STATS] = ttest(X(:,icond),X(:,jcond));
    effd = STATS.tstat/sqrt(size(X,1));
    % identify row corresponding to post-hoc ttest in the multcompare
    % output
    irow = double(c.condition_1) == icond & double(c.condition_2) == jcond;
    P  = c.pValue(irow);
    CI = [c.Lower(irow) c.Upper(irow)];
    H  = P < 0.05;
    bf10 = bf.ttest(X(:,icond),X(:,jcond));
    if H==0, BFval=1/bf10; BFstr='BF01';else,BFval=bf10;BFstr='BF10';end
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+06.2f ; p(corrected) = %4.4f ; Cohen''s d = %+04.4f ; CI(corrected) = [%+4.4f %+4.4f] ; %s = %4.4f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end
fprintf('\n')

% bias variance decomp - variance
fprintf('<strong>bias-variance decomposition - VARIANCE - ANOVA</strong>\n')
cfg = cfg_fig;
cfg.pbar = .5;
cfg.HGT = HGT;
cfg.wdth = 0.5;
cfg.figname = 'bias_variance_siginf_variance';
cfg.rgb = rgb(ordercond(1:2),:);
cfg.x = bv(:,ordercond(1:2),2);
cfg.plotx = true;
cfg.ttest = true;
cfg.xlabel = 'condition';
cfg.labels = condnam(ordercond(1:2));
cfg.ylabel = 'variance term (s.d.)';
cfg.ylim = [0 2];cfg.ytick = 0:.5:2;
helper_plot_overall_bars(cfg);
plot(xlim,mean(bv(:,3,2))*[1 1],'Color',rgb(3,:),'Linewidth',1);
helper_ILL(gcf,cfg.pbar,cfg.HGT,false);
if cfg_fig.saveornot
    print(fullfile(cfg_fig.figpath,get(gcf,'Name')),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',get(gcf,'Name'),cfg_fig.figpath)
end

X = bv(:,:,2);
[~,rm] = anova_1way_ranova(X, 'condition'); 
c = multcompare(rm,'condition'); % correct for multiple comparisons
fprintf('post-hoc ttest:\n')
for icond = 1:ncond
    jcond = mod(icond,ncond)+1;
    [~,~,~,STATS] = ttest(X(:,icond),X(:,jcond));
    effd = STATS.tstat/sqrt(size(X,1));
    % identify row corresponding to post-hoc ttest in the multcompare
    % output
    irow = double(c.condition_1) == icond & double(c.condition_2) == jcond;
    P  = c.pValue(irow);
    CI = [c.Lower(irow) c.Upper(irow)];
    H  = P < 0.05;
    bf10 = bf.ttest(X(:,icond),X(:,jcond));
    if H==0, BFval=1/bf10; BFstr='BF01';else,BFval=bf10;BFstr='BF10';end
    fprintf('\tcond%d vs cond%d : H = %d ; t = %+06.2f ; p(corrected) = %4.4f ; Cohen''s d = %+04.4f ; CI(corrected) = [%+4.4f %+4.4f] ; %s = %4.4f\n',icond,jcond,H,STATS.tstat,P,effd,CI(1),CI(2),BFstr,BFval);
end
fprintf('\n')

% plot horizontal bars
fprintf('<strong>R-squared estimated from fraction consistent for each condition</strong>\n')
cfg = cfg_fig;
cfg.pbar = 1.2;
cfg.HGT = HGT;
cfg.wdth = 0.3;
cfg.figname = 'bias_variance_rsq';
cfg.rgb = rgb(ordercond,:);
cfg.x = rsq_dat(:,ordercond);
cfg.plotx = true;
cfg.xlabel = 'condition';
cfg.labels = condnam(ordercond);
cfg.ylabel = 'r-squared';
cfg.xlim = [0 4];
cfg.ylim = [0 1];cfg.ytick = 0:.2:1;
cfg.type = 'horizontal';
helper_plot_overall_bars(cfg);
fprintf('\n')
for icond = 1:3
    fprintf('cond %s:\tvariance term (1-R2) = %4.3f ± %4.3f %%\n',condnam{icond},100-mean(rsq_dat(:,icond))*100,std(rsq_dat(:,icond))/sqrt(ndata)*100)
end
stat = anova_1way_ranova(rsq_dat, 'condition');