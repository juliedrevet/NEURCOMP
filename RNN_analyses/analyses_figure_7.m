%% analyses_figure_7
% this script runs the statistics and plots of Figure 7 (RNNs trained in PRE and POST) from the results of
% the Python analyses (python/scripts/analyses_PRE_POST.ipynb, exported in .mat files by its export cells):
%   (b) accuracy along training (PRE vs POST, cluster-based permutation test), compression, and partial
%       correlation between compression and accuracy controlling for the coding of the mean orientation
%       (cluster-based permutation test)
%   (c) resilience to computation noise: accuracy difference PRE - POST (one-tailed cluster-based permutation
%       test), mean accuracy, histograms of human accuracy and of the noise matching it, coding precision on
%       correct and error sequences
%   (d) causal perturbation: accuracy drop, change in mean orientation coding and in compression (paired t-tests)
% arrays are networks x (epochs or noise levels) x conditions (1:PRE | 2:POST)

clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

sigcat = 15;
condnam = {'PRE','POST'};

figpath = './figs';
[~,~] = mkdir(figpath);
saveornot = false;

cfg_fig = [];
cfg_fig.figpath = figpath;
cfg_fig.saveornot = saveornot;

rgb = lines;

alphathresh = 0.05;

nrnn = Inf; % number of networks used (networks 1:nrnn), Inf: all networks

% cluster-based permutation tests
cfg_stat = [];
cfg_stat.alpha_form = 0.05;
cfg_stat.nperm = 1e4;
cfg_stat.seed = 1234;
cfg_stat.verbose = false;

nres = 1e4; % number of bootstrap resamples to compute the s.e.m. of the partial correlation

% results of the Python analyses (networks trained without noise)
respath = fullfile('python','rnn_simulations','RNN_train_PRE_POST',sprintf('sigcat_%d',sigcat));

%% Figure 7b: along training
load(fullfile(respath,'test_epochs.mat'))
N = min(nrnn,size(acc,1)); % networks 1:N
fprintf('%d/%d networks (sigcat = %d), %d saved epochs (%d to %d)\n',N,size(acc,1),sigcat,numel(epochs),epochs(1),epochs(end))
acc = acc(1:N,:,:);
acc_bayes = acc_bayes(1:N);
compression = compression(1:N,:,:);
cumang = cumang(1:N,:,:);

pbar = 1.2;
HGT = 4.5;

% accuracy along training
figname = sprintf('PRE_POST_accuracy_per_epoch_sigcat%d',sigcat);
figure('Color','white','Name',figname);
hold on
% plot single RNN
for irnn = 1:N
    for icond = 1:2
        plot(epochs,acc(irnn,:,icond),'color',(1+rgb(icond,:))/2)
    end
end
plot([0 epochs(end)],[.5 .5],'Color',[.8 .8 .8])
for icond = 1:2
    plot(epochs,mean(acc(:,:,icond)),'color',rgb(icond,:),'linewidth',2)
end
plot([0 epochs(end)],mean(acc_bayes)*[1 1],'k:') % Bayes-optimal accuracy
% cluster-based permutation test on PRE - POST
cfg = cfg_stat;
cfg.tail = 'both';
cfg.xvec = epochs;
stat_acc = cluster_permutation(acc(:,:,1)-acc(:,:,2),[],cfg);
% significant clusters (PRE - POST): PRE > POST in blue, POST > PRE in orange
fprintf('<strong>accuracy along training: %s vs %s</strong>\n',condnam{1},condnam{2})
fprintf('\tBayes-optimal accuracy: %.3f\n',mean(acc_bayes))
% final accuracy (last saved epoch)
fprintf('\tfinal accuracy (epoch %d): %s %.2f ± %.2f %%, %s %.2f ± %.2f %% (mean ± s.e.m.)\n',epochs(end), ...
    condnam{1},100*mean(acc(:,end,1)),100*std(acc(:,end,1))/sqrt(N), ...
    condnam{2},100*mean(acc(:,end,2)),100*std(acc(:,end,2))/sqrt(N))
for ic = 1:numel(stat_acc.clusters)
    c = stat_acc.clusters(ic);
    jcond = 1+(c.tstat<0); % condition with the higher accuracy
    fprintf('\tepochs %d-%d: %s > %s, cluster tstat = %+.3f, dz = %+.2f, p = %.5f\n', ...
        c.tbeg,c.tend,condnam{jcond},condnam{3-jcond},c.tstat,mean(stat_acc.t(c.idx))/sqrt(N),c.p)
    if c.p < alphathresh
        plot([c.tbeg c.tend],[.97 .97],'Color',rgb(jcond,:),'LineWidth',3)
    end
end
xlim([0 epochs(end)])
ylim([0.4,1]);
set(gca,'Xtick',200:200:1000)
xlabel('training epoch');
ylabel('decision accuracy');
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% compression along training
figname = sprintf('PRE_POST_compression_per_epoch_sigcat%d',sigcat);
figure('Color','white','Name',figname)
hold on
% plot single RNN
for irnn = 1:N
    for icond = 1:2
        plot(epochs,compression(irnn,:,icond),'color',(1+rgb(icond,:))/2)
    end
end
plot([0 epochs(end)],[0 0],'Color',[.8 .8 .8])
for icond = 1:2
    plot(epochs,mean(compression(:,:,icond)),'Color',rgb(icond,:),'LineWidth',1.5)
end
xlim([0 epochs(end)])
set(gca,'Xtick',200:200:1000)
xlabel('training epoch')
ylabel('compression')
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% partial correlation between compression and accuracy (controlling for the coding of the mean orientation), PRE
icond = 1;
% partial correlation across networks and cluster-based permutation test (permutation of networks)
cfg = cfg_stat;
cfg.tail = 'both';
cfg.xvec = epochs;
cfg.Z = cumang(:,:,icond);
stat_par = cluster_permutation(acc(:,:,icond),compression(:,:,icond),cfg);
r_par = stat_par.r;
% s.e.m.: bootstrap over networks
rs = RandStream('mt19937ar','Seed',1234);
r_boot = nan(nres,numel(epochs));
for ires = 1:nres
    jres = randi(rs,N,N,1);
    r_boot(ires,:) = partial_corr(compression(jres,:,icond),acc(jres,:,icond),cumang(jres,:,icond));
end
r_sem = std(r_boot,[],1);
figname = sprintf('partial_corr_comp_cumang_sigcat%d',sigcat);
figure('Color','white','Name',figname);
hold on
% s.e.m. (bootstrap over networks)
xvec = epochs;
patch([xvec,fliplr(xvec)],[r_par+r_sem,fliplr(r_par-r_sem)], ...
    (1+rgb(icond,:))/2,'EdgeColor','none','FaceAlpha',.5);
plot([0 epochs(end)],[0 0],'Color',[.8 .8 .8]);
fprintf('<strong>%s partial correlation between compression and accuracy</strong>\n',condnam{icond})
for ic = 1:numel(stat_par.clusters)
    c = stat_par.clusters(ic);
    fprintf('\tepochs %d-%d: mean partial r = %.2f, cluster tstat = %+.3f, p = %.5f\n', ...
        c.tbeg,c.tend,mean(r_par(c.idx)),c.tstat,c.p)
    if c.p < alphathresh
        plot([c.tbeg c.tend],[0 0],'Color',rgb(icond,:),'LineWidth',3)
    end
end
plot(epochs,r_par,'Color',rgb(icond,:),'LineWidth',2);
xlim([0 epochs(end)])
set(gca,'Xtick',200:200:1000)
xlabel('training epoch')
ylabel('partial correlation')
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

%% Figure 7c: resilience to computation noise
load(fullfile(respath,'test_noise.mat'))
N = min(nrnn,size(acc,1)); % networks 1:N
acc = acc(1:N,:,:);
rdec_cor = structfun(@(x) x(1:N,:,:),rdec_cor,'UniformOutput',false);
rdec_err = structfun(@(x) x(1:N,:,:),rdec_err,'UniformOutput',false);

HGT = 4;
pbar = 1.5;

% accuracy difference PRE - POST
figname = sprintf('PRE_POST_exact_noise_injection_sigcat%d',sigcat);
figure('Color','white','Name',figname)
hold on
acc_diff = (acc(:,:,1)-acc(:,:,2))*100;
yavg = mean(acc_diff);
yerr = std(acc_diff)/sqrt(N);
xvec = noise_array;
patch([xvec,fliplr(xvec)],[yavg+yerr,fliplr(yavg-yerr)], ...
    [0.8,0.8,0.8],'EdgeColor','none');
plot(noise_array,yavg,'Color',[.5,.5,.5]/2,'LineWidth',1.5)
plot([min(noise_array) max(noise_array)],[0,0],'k:')
% cluster-based permutation test on PRE - POST, one-tailed (PRE > POST)
cfg = cfg_stat;
cfg.tail = 'right';
cfg.xvec = noise_array;
stat_diff = cluster_permutation(acc_diff,[],cfg);
% significant clusters (one-tailed: PRE > POST)
fprintf('<strong>accuracy difference %s - %s with computation noise</strong>\n',condnam{1},condnam{2})
for ic = 1:numel(stat_diff.clusters)
    c = stat_diff.clusters(ic);
    fprintf('\tnoise %.2f-%.2f: %s > %s, cluster tstat = %+.3f, dz = %+.2f, p = %.5f\n', ...
        c.tbeg,c.tend,condnam{1},condnam{2},c.tstat,mean(stat_diff.t(c.idx))/sqrt(N),c.p)
    if c.p < alphathresh
        plot([c.tbeg c.tend],[0 0],'Color',rgb(1,:),'LineWidth',2)
    end
end
xlim([min(noise_array) max(noise_array)])
xlabel('neural noise')
ylabel(sprintf('accuracy\ndifference (%%)'))
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% mean accuracy (PRE and POST)
figname = sprintf('PRE_POST_mean_accuracy_sigcat_%d',sigcat);
figure('Color','white','name',figname)
hold on;
% plot single RNN
for irnn = 1:N
    plot(noise_array,mean(acc(irnn,:,:),3),'color',[.8,.8,.8]);
end
yavg = mean(mean(acc,3));
yerr = std(mean(acc,3))/sqrt(N);
xvec = noise_array;
patch([xvec,fliplr(xvec)],[yavg+yerr,fliplr(yavg-yerr)], ...
    [0.8,0.8,0.8],'EdgeColor','none');
plot(noise_array,yavg,'Color','black','linewidth',1.5)
xlim([min(noise_array) max(noise_array)])
ylim([0.5,0.85])
xlabel('neural noise')
ylabel(sprintf('mean\naccuracy'))
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% histogram of human accuracy to add in illustrator
% behavioral data table organized per sequence
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');

participants = unique(dat.participant);
ndata = length(participants);

perf = nan(ndata,2);
for icond = 1:2
    for idata = 1:ndata
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        resp = dat.response(ifilt);
        seqcat = dat.seqcat(ifilt);
        perf(idata,icond) = mean(resp == seqcat);
    end
end

figname = 'histogram_perf_mean';
figure('color','white','name',figname)
histogram(mean(perf,2),'facecolor',[0.5 0.5 0.5],'edgecolor','w')
xlim([0.5,0.85]);
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% histogram of the computation noise matching human accuracy, to add in illustrator
% for each participant and each network: noise level at which the mean accuracy (PRE and POST) of the network
% equals the mean accuracy (PRE and POST) of the participant, and accuracy difference PRE - POST of the network
% at this noise level
perf_humans = mean(perf,2);
noise_dense = linspace(noise_array(1),noise_array(end),10000);
sigma_matched = nan(ndata,N);
diff_matched = nan(ndata,N);
for irnn = 1:N
    acc_dense = interp1(noise_array,mean(acc(irnn,:,:),3),noise_dense,'pchip');
    diff_dense = interp1(noise_array,acc_diff(irnn,:),noise_dense,'pchip');
    for idata = 1:ndata
        % closest accuracy over the whole noise range (the accuracy curve may not be monotonic)
        [~,imin] = min(abs(acc_dense-perf_humans(idata)));
        sigma_matched(idata,irnn) = noise_dense(imin);
        diff_matched(idata,irnn) = diff_dense(imin);
    end
end
fprintf('computation noise matching human accuracy: median = %.3f, interquartile range = %.3f-%.3f, 5th-95th percentiles = %.3f-%.3f\n', ...
    median(sigma_matched(:)),quantile(sigma_matched(:),[.25 .75]),quantile(sigma_matched(:),[.05 .95]))
x = mean(diff_matched,1); % mean over participants, for each network
fprintf('accuracy difference %s - %s at the noise matching human accuracy: %.2f ± %.2f %% (mean ± s.e.m.)\n', ...
    condnam{1},condnam{2},mean(x),std(x)/sqrt(N))
fprintf('accuracy difference %s - %s in humans: %.2f %%\n',condnam{1},condnam{2},100*mean(perf(:,1)-perf(:,2)))

figname = sprintf('histogram_computation_noise_sigcat%d',sigcat);
figure('color','white','name',figname)
histogram(sigma_matched(:),'facecolor',[0.5 0.5 0.5],'edgecolor','w')
xlim([0 1]);
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% coding precision on correct (solid) and error (dashed) sequences
icond = 1;
varnam = {'cumang','cumnul','cumllr'};
varlab = {'mean orientation','null belief','true belief'};
fprintf('<strong>coding precision on correct and error sequences with computation noise (mean over noise levels)</strong>\n')
for ivar = 1:numel(varnam)
    figname = sprintf('PRE_POST_noise_injection_%s_corerr_sigcat%d',varnam{ivar},sigcat);
    figure('Color','white','Name',figname)
    hold on
    y1 = rdec_cor.(varnam{ivar})(:,:,icond);
    y2 = rdec_err.(varnam{ivar})(:,:,icond);
    xvec = noise_array;
    % plot single RNN
    for irnn = 1:N
        plot(xvec,y1(irnn,:),'color',(1+rgb(icond,:))/2)
        plot(xvec,y2(irnn,:),'color',(1+rgb(icond,:))/2)
    end
    plot(xvec,mean(y1),'Color',rgb(icond,:),'LineWidth',1.5)
    plot(xvec,mean(y2),'Color',rgb(icond,:),'LineWidth',1.5,'Linestyle','--')
    xlim([min(noise_array) max(noise_array)])
    ylim([-0.1 1])
    ylabel(sprintf('coding precision (r)'))
    xlabel('neural noise')
    helper_ILL(gcf,pbar,HGT);
    if saveornot
        print(fullfile(figpath,figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',figname,figpath)
    end
    fprintf('\t%s %s: difference between correct and error trials: %4.3f ± %4.3f\n', ...
        condnam{icond},varlab{ivar},mean(mean(y1-y2,2)),std(mean(y1-y2,2))/sqrt(N))
end

%% Figure 7d: causal perturbation
load(fullfile(respath,'test_perturbation.mat'))
N = min(nrnn,size(acc0,1)); % networks 1:N
eff = structfun(@(x) x(1:N,:),eff,'UniformOutput',false);

% effects averaged over the perturbation strengths (networks x bars):
% POST mean orientation, PRE mean orientation, PRE true belief, PRE null belief
barnam = {'POST mean orient.','PRE mean orient.','PRE true belief','PRE null belief'};

cfg = cfg_fig;
cfg.plotx = true;
cfg.HGT = 2.5;
cfg.pbar = 3;
cfg.width = 0.6;
cfg.rgb = [rgb(2,:);rgb(1,:);rgb(1,:);rgb(1,:)];
cfg.xlabel = 'perturbation';
cfg.xticklabel = {'mean orient.','mean orient.','true belief','null belief'};

% accuracy drop
cfg.figname = sprintf('causal_perturbation_PRE_POST_bars_accuracy_drop_sigcat%d_percentage',sigcat);
cfg.x = eff.accuracy_drop;
cfg.ylim = [-25,0];
cfg.ylabel = sprintf('accuracy drop (%%)');
helper_plot_overall_bars(cfg)

% mean orientation coding change
cfg.figname = sprintf('causal_PRE_any_POST_cumang_bars_cumang_drop_sigcat%d',sigcat);
cfg.x = eff.mean_orient_coding_change;
cfg.ylim = [-.4,0];
cfg.ylabel = sprintf('mean orient. coding change (r)');
helper_plot_overall_bars(cfg)

% compression change
cfg.figname = sprintf('causal_PRE_any_POST_cumang_bars_compression_drop_sigcat%d',sigcat);
cfg.x = eff.compression_change;
cfg.ylim = [-.4,0];
cfg.ylabel = sprintf('compression change (r)');
helper_plot_overall_bars(cfg)

% compression change
cfg.figname = sprintf('causal_PRE_any_POST_cumang_bars_compression_drop_sigcat%d_nullbar',sigcat);
cfg.x = eff.compression_change;
cfg.ylim = [-.2,.2];
cfg.ylabel = sprintf('compression change (r)');
helper_plot_overall_bars(cfg)

fprintf('<strong>causal perturbation (mean ± s.e.m.)</strong>\n')
effnam = fieldnames(eff);
for ieff = 1:numel(effnam)
    x = eff.(effnam{ieff});
    fprintf('\t%s:\n',strrep(effnam{ieff},'_',' '))
    for ibar = 1:numel(barnam)
        fprintf('\t\t%s: %+.3f ± %.3f\n',barnam{ibar},mean(x(:,ibar)),std(x(:,ibar))/sqrt(N))
    end
end

% accuracy drop: POST mean orientation vs PRE mean orientation and vs PRE true belief
x = eff.accuracy_drop;
for ibar = 2:3
    [~,PVAL,~,STATS] = ttest(x(:,1),x(:,ibar));
    bf10 = bf.ttest('T',STATS.tstat,'N',N);
    if PVAL < 0.05, BFval = bf10; BFstr = 'BF10'; else, BFval = 1/bf10; BFstr = 'BF01'; end
    dz = mean(x(:,1)-x(:,ibar))/std(x(:,1)-x(:,ibar));
    fprintf('\tPOST-HOC ttest accuracy drop %s vs %s: t(%d) = %+.5f ; p = %.5f ; dz = %+.2f ; %s = %.5f\n', ...
        barnam{1},barnam{ibar},STATS.df,STATS.tstat,PVAL,dz,BFstr,BFval)
end

function r = partial_corr(X,Y,Z)
% partial correlation between each column of X and the same column of Y, controlling for the
% same column of Z (rows = networks): correlation between the residuals of the linear regressions
% (with intercept) of X and Y on Z
X = X-mean(X,1);
Y = Y-mean(Y,1);
Z = Z-mean(Z,1);
X = X-Z.*(sum(Z.*X,1)./sum(Z.^2,1));
Y = Y-Z.*(sum(Z.*Y,1)./sum(Z.^2,1));
r = sum(X.*Y,1)./sqrt(sum(X.^2,1).*sum(Y.^2,1));
end
