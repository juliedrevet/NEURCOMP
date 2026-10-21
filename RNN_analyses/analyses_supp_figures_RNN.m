%% analyses_supp_figures_RNN
% this script runs the statistics and plots of the supplementary figures of the RNNs:
%   Supp. Fig. 11a: accuracy of the computational model as a function of the noise on the category axis,
%       simulated with the best-fitting parameters of each participant in REF (BEH_analyses/results/results_fit.mat),
%       compared to participants' accuracy in PRE and REF
% and, from the results of the Python analyses of the RNNs trained in PRE and POST (exported in .mat files
% by the export cells of the notebooks):
%   Supp. Fig. 12: resilience to input noise (python/scripts/analyses_PRE_POST_SUPP.ipynb, test_input_noise.mat):
%       accuracy difference PRE - POST (one-tailed cluster-based permutation test), mean accuracy, coding
%       precision on correct and error sequences
%   Supp. Fig. 13: causal perturbations as a function of the perturbation strength
%       (python/scripts/analyses_PRE_POST.ipynb, test_perturbation.mat): accuracy, mean orientation coding
%       and compression
% arrays of Supp. Figs. 12-13 are networks x (noise levels or perturbation strengths) x conditions (1:PRE | 2:POST)

clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

sigcat = 15;
condnam = {'PRE','POST'};

figpath = './figs';
[~,~] = mkdir(figpath);
saveornot = false;

rgb = lines;

alphathresh = 0.05;

nrnn = Inf; % number of networks used (networks 1:nrnn), Inf: all networks

% cluster-based permutation tests
cfg_stat = [];
cfg_stat.alpha_form = 0.05;
cfg_stat.nperm = 1e4;
cfg_stat.seed = 1234;
cfg_stat.verbose = false;

% results of the Python analyses
respath = fullfile(cfg_path.rnnsim,'RNN_train_PRE_POST',sprintf('sigcat_%d',sigcat));

%% Noise on the category axis
% for each participant and each noise level sigcat (0:2:20°, s.d.), the computational model (inference imbalance
% gamma and inference noise siginf of REF, lapse epsi) is simulated on 1000 sequences of each length (4, 8, 12),
% the category axis being retrieved with Gaussian noise at each stimulus; same sequences and same noise samples
% for all sigcat (seed 100*participant)

% best-fitting parameters of the winning model (BEH_analyses/run_fit.m)
load(fullfile(cfg_path.results,'results_fit.mat'),'out_fit')
ndata = size(out_fit,1);
xavg = cell2mat(cellfun(@(x)getfield(x,'xavg'),out_fit,'UniformOutput',false));
gamma  = xavg(:,contains(out_fit{1}.xnam,'gamma'));
gamma  = gamma(:,3);  % REF
siginf = xavg(:,contains(out_fit{1}.xnam,'siginf'));
siginf = siginf(:,3); % REF
epsi   = xavg(:,contains(out_fit{1}.xnam,'epsi'));

% human accuracy (1:PRE | 2:POST | 3:REF)
load(fullfile(cfg_path.data,'data_table_seq.mat'),'dat');
participants = unique(dat.participant);
perf = nan(ndata,3);
for icond = 1:3
    for idata = 1:ndata
        ifilt = dat.participant == participants(idata) & dat.condtn == icond;
        perf(idata,icond) = mean(dat.response(ifilt) == dat.seqcat(ifilt));
    end
end

sigcat_array = 0:2:20; % s.d. (degrees) of the noise on the category axis
nsigcat = numel(sigcat_array);
nsmp = 5;   % simulations of each sequence
sigsel = 0; % no selection noise
acc_sim = nan(ndata,nsigcat);
for idata = 1:ndata
    fprintf('simulating participant %02d/%02d ...\n',idata,ndata)
    for isig = 1:nsigcat
        blck = generate_equalized_sequences(1000,[4,8,12],sigcat_array(isig),100*idata);
        llrhat = blck.llrhat;
        seqcat = blck.seqcat;
        nseq = numel(seqcat);
        x = zeros(nseq,nsmp);
        r = zeros(nseq,nsmp);
        for iseq = 1:nseq
            for itlt = 1:numel(llrhat{iseq})
                x(iseq,:) = exp(gamma(idata))*x(iseq,:) + llrhat{iseq}(itlt) + randn(1,nsmp)*siginf(idata);
            end
            % response lapse
            i0 = rand(1,nsmp) < epsi(idata);
            r(iseq,i0) = sign(rand(1,nnz(i0))-.5); % random response
            r(iseq,~i0) = sign(x(iseq,~i0)+sigsel*randn(1,nnz(~i0)));
        end
        r(r<0) = 2; % express response as 1 or 2
        acc_sim(idata,isig) = mean(mean(r == seqcat,1),2);
    end
end

figname = 'sigcat_simulation';
figure('Color','white','Name',figname)
hold on
for idata = 1:ndata
    plot(sigcat_array,acc_sim(idata,:),'Color',[.8 .8 .8],'LineWidth',1)
end
plot(sigcat_array,mean(acc_sim),'Color',[.5 .5 .5],'LineWidth',2)
plot(xlim,mean(perf(:,3))*[1 1],'Color',[.7 .7 .7],'LineWidth',2,'LineStyle',':')
plot(xlim,mean(perf(:,1))*[1 1],'Color',rgb(1,:),'LineWidth',2,'LineStyle',':')
ylim([0.69,0.85])
xlabel('noise on category axis (°)')
ylabel('accuracy')
helper_ILL(gcf,1.3,4);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

%% Resilience to input noise
load(fullfile(respath,'test_input_noise.mat'))
N = min(nrnn,size(acc,1));
acc = acc(1:N,:,:);
rdec_cor = structfun(@(x) x(1:N,:,:),rdec_cor,'UniformOutput',false);
rdec_err = structfun(@(x) x(1:N,:,:),rdec_err,'UniformOutput',false);
noise_array = input_noise_array; % s.d. (degrees) of the noise on the orientation of each stimulus

HGT = 4;
pbar = 1.5;

% accuracy difference PRE - POST
figname = sprintf('PRE_POST_input_noise_injection_sigcat%d',sigcat);
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
fprintf('<strong>accuracy difference %s - %s with input noise</strong>\n',condnam{1},condnam{2})
if isempty(stat_diff.clusters)
    fprintf('\tno cluster where %s > %s\n',condnam{1},condnam{2})
end
for ic = 1:numel(stat_diff.clusters)
    c = stat_diff.clusters(ic);
    fprintf('\tinput noise %.0f-%.0f: %s > %s, cluster tstat = %+.3f, dz = %+.2f, p = %.5f\n', ...
        c.tbeg,c.tend,condnam{1},condnam{2},c.tstat,mean(stat_diff.t(c.idx))/sqrt(N),c.p)
    if c.p < alphathresh
        plot([c.tbeg c.tend],[0 0],'Color',rgb(1,:),'LineWidth',2)
    end
end
xlim([min(noise_array) max(noise_array)])
ylim([-3,1])
xlabel('input noise (°)')
ylabel(sprintf('accuracy\ndifference (%%)'))
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% mean accuracy (PRE and POST)
% (histogram of human accuracy to add in illustrator: same as Figure 7c, histogram_perf_mean in analyses_figure_7.m)
figname = sprintf('PRE_POST_mean_accuracy_sigcat_%d_input_noise',sigcat);
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
xlabel('input noise (°)')
ylabel(sprintf('mean\naccuracy'))
helper_ILL(gcf,pbar,HGT);
if saveornot
    print(fullfile(figpath,figname),'-vector','-dpdf');
    fprintf('\nfigure %s saved in %s !\n',figname,figpath)
end

% coding precision on correct (solid) and error (dashed) sequences, in PRE
icond = 1;
varnam = {'cumang','cumnul','cumllr'};
varlab = {'mean orientation','null belief','true belief'};
fprintf('<strong>coding precision on correct and error sequences with input noise (mean over noise levels)</strong>\n')
for ivar = 1:numel(varnam)
    figname = sprintf('PRE_POST_noise_injection_%s_corerr_sigcat%d_input_noise',varnam{ivar},sigcat);
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
    xlabel('input noise (°)')
    helper_ILL(gcf,pbar,HGT);
    if saveornot
        print(fullfile(figpath,figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',figname,figpath)
    end
    fprintf('\t%s %s: difference between correct and error trials: %4.3f ± %4.3f\n', ...
        condnam{icond},varlab{ivar},mean(mean(y1-y2,2)),std(mean(y1-y2,2))/sqrt(N))
end

%% Causal perturbations
load(fullfile(respath,'test_perturbation.mat'))
N = min(nrnn,size(acc0,1));
acc0 = acc0(1:N,:);
acc_pert = structfun(@(x) x(1:N,:,:),acc_pert,'UniformOutput',false);
rdec_pert = structfun(@(s) structfun(@(x) x(1:N,:,:),s,'UniformOutput',false),rdec_pert,'UniformOutput',false);
% coding of the mean orientation under its own perturbation: cumang_matched = 1: coding of the perturbed component
% only (cumsumcos under the cumsumcos perturbation, cumsumsin under the cumsumsin perturbation)
if cumang_matched
    fprintf('mean orientation coding under its own perturbation: perturbed component only\n')
else
    fprintf('mean orientation coding under its own perturbation: both components\n')
end

HGT = 4;
pbar = 1.5;

% perturbation of the mean orientation, of the true belief and of the null belief
pertnam = {'cumang','cumllr','cumnul'};
pertlab = {'mean orientation','true belief','null belief'};
xvec = sel_noise_array;
for ipert = 1:numel(pertnam)
    p = pertnam{ipert};

    % accuracy, in % of the accuracy of the unperturbed network
    figname = sprintf('PRE_POST_causal_perturbation_%s_accuracy_sigcat%d',p,sigcat);
    figure('Color','white','Name',figname)
    hold on
    for icond = 1:2
        y = acc_pert.(p)(:,:,icond)./acc0(:,icond)*100;
        for irnn = 1:N
            plot(xvec,y(irnn,:),'color',(1+rgb(icond,:))/2)
        end
    end
    for icond = 1:2
        y = acc_pert.(p)(:,:,icond)./acc0(:,icond)*100;
        plot(xvec,mean(y),'Color',rgb(icond,:),'LineWidth',1.5)
    end
    plot([min(xvec) max(xvec)],[100 100],'k:')
    xlim([min(xvec) max(xvec)])
    ylim([50 100])
    ylabel('accuracy (% optimal)')
    xlabel('perturbation strength')
    helper_ILL(gcf,pbar,HGT);
    if saveornot
        print(fullfile(figpath,figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',figname,figpath)
    end

    % coding precision of the mean orientation, and compression
    resnam = {'cumang','compression'};
    reslab = {sprintf('mean orientation\ncoding precision (r)'),'compression (r)'};
    for ires = 1:numel(resnam)
        figname = sprintf('PRE_POST_causal_perturbation_%s_%s_sigcat%d',p,resnam{ires},sigcat);
        figure('Color','white','Name',figname)
        hold on
        for icond = 1:2
            y = rdec_pert.(p).(resnam{ires})(:,:,icond);
            for irnn = 1:N
                plot(xvec,y(irnn,:),'color',(1+rgb(icond,:))/2)
            end
        end
        for icond = 1:2
            y = rdec_pert.(p).(resnam{ires})(:,:,icond);
            plot(xvec,mean(y),'Color',rgb(icond,:),'LineWidth',1.5)
        end
        xlim([min(xvec) max(xvec)])
        ylim([0 1])
        ylabel(reslab{ires})
        xlabel('perturbation strength')
        helper_ILL(gcf,pbar,HGT);
        if saveornot
            print(fullfile(figpath,figname),'-vector','-dpdf');
            fprintf('\nfigure %s saved in %s !\n',figname,figpath)
        end
    end

    % values at the strongest perturbation (mean ± s.e.m.)
    fprintf('<strong>perturbation of the %s (strength %g)</strong>\n',pertlab{ipert},xvec(end))
    for icond = 1:2
        a = acc_pert.(p)(:,end,icond)./acc0(:,icond)*100;
        m = rdec_pert.(p).cumang(:,end,icond);
        c = rdec_pert.(p).compression(:,end,icond);
        fprintf('\t%s: accuracy %.1f ± %.1f %% optimal, mean orientation coding %.3f ± %.3f, compression %.3f ± %.3f\n', ...
            condnam{icond},mean(a),std(a)/sqrt(N),mean(m),std(m)/sqrt(N),mean(c),std(c)/sqrt(N))
    end
end
