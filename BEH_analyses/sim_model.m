function [out] = sim_model(cfg)
%  SIM_MODEL  Simulate model
%
%  This function simulates the model of fit_model.m:
%    * seqcat  => sequence category (1|2)
%    * tltllr  => card LLRs
%    * condtn  => condition (1:PRE|2:POST|3:REF)
%    * gamma1  => inference imbalance in PRE
%    * gamma2  => inference imbalance in POST
%    * gamma3  => inference imbalance in REF
%    * siginf1 => inference noise in PRE
%    * siginf2 => inference noise in POST
%    * siginf3 => inference noise in REF
%    * sigsel1 => selection noise in PRE
%    * sigsel2 => selection noise in POST
%    * sigsel3 => selection noise in REF
%    * epsi    => lapse rate (common to all conditions)
%    * nsmp    => number of simulations (default 1000)
%    * seed    => random seed (default from clock)

% check configuration parameters
if ~all(isfield(cfg,{'seqcat','tltllr','condtn',...
                     'gamma1','gamma2','gamma3',...
                     'siginf1','siginf2','siginf3',...
                     'sigsel1','sigsel2','sigsel3','epsi'}))
    error('Incomplete configuration structure!');
end
if ~isfield(cfg,'nsmp')
    cfg.nsmp = 1e3;
end
if ~isfield(cfg,'seed') % seed for random number generator
    cfg.seed = round(sum(100*clock)); % seed based on current time
end
% initialize random number generator
sd = RandStream('mt19937ar','Seed',cfg.seed);
RandStream.setGlobalStream(sd);

gamma  = [cfg.gamma1 cfg.gamma2 cfg.gamma3];
siginf = [cfg.siginf1 cfg.siginf2 cfg.siginf3];
sigsel = [cfg.sigsel1 cfg.sigsel2 cfg.sigsel3];
epsi = cfg.epsi;

seqcat = cfg.seqcat; % sequence category (1|2)
tltllr = cfg.tltllr;
condtn = cfg.condtn;
nsmp   = cfg.nsmp;

seqllr = cellfun(@sum,tltllr);

seqlen = cellfun(@length,tltllr);
nseq   = numel(seqcat);
ncond  = numel(unique(condtn));

x = zeros(nseq,nsmp);
r = zeros(nseq,nsmp);

for iseq = 1:nseq
    icond = condtn(iseq);
    for itlt = 1:seqlen(iseq)
        % update sequence (with inference imbalance and noise)
        x(iseq,:) = exp(gamma(icond))*x(iseq,:) + tltllr{iseq}(itlt) + randn(1,nsmp)*siginf(icond);
    end
    % response lapse
    i0 = rand(1,nsmp) < epsi;
    r(iseq,i0) = sign(rand(1,nnz(i0))-.5); % random response
    % with selection noise
    r(iseq,~i0) = sign(x(iseq,~i0)+sigsel(icond)*randn(1,nnz(~i0)));
end

r(r<0) = 2; % express response as 1 or 2

pcor_avg = nan(1,ncond);
pcor_std = nan(1,ncond);

pcor_opt = nan(1,ncond);
for icond = 1:ncond
    ifilt = condtn == icond;
    pcor_avg(icond) = mean(mean(r(ifilt,:) == seqcat(ifilt,:),1),2);
    pcor_std(icond) = std(mean(r(ifilt,:) == seqcat(ifilt,:),1),[],2);
    pcor_opt(icond) = mean(mean(sign(1.5-r(ifilt,:)) == sign(seqllr(ifilt,:)),1),2);
end

pcor_avg_all = mean(mean(r == seqcat,1),2);
pcor_std_all = std(mean(r == seqcat,1),[],2);

out = cfg;
out.r = r;
out.x = x;
out.pcor_avg = pcor_avg;
out.pcor_std = pcor_std;
out.pcor_avg_all = pcor_avg_all;
out.pcor_std_all = pcor_std_all;
out.pcor_opt = pcor_opt;

end