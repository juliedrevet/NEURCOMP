function out = fit_model_sensory(cfg)
% fit model to the participants' data
% parameters to be fitted:
%           * gamma  : inference imbalance (one per condition)
%           * siginf : inference noise (one per condition)
%           * sigsel : selection noise (one per condition)
%           * epsi   : lapse rate (common to all conditions)
% SENSORY INFORMATION RATHER THAN CATEGORIZED INFORMATION (no LOGLR)

% check configuration parameters
if ~all(isfield(cfg,{'seqang','catang','seqcat','condtn','resp'}))
    error('Incomplete configuration structure!');
end
if ~isfield(cfg,'verbose')
    cfg.verbose = false;
end
if ~isfield(cfg,'kappa')
    cfg.kappa = 0.5;
end
if ~isfield(cfg,'nsmp')
    cfg.nsmp = 1e3;
end
if ~isfield(cfg,'nres')
    cfg.nres = 1e3;
end
if ~isfield(cfg,'seed')
    cfg.seed = round(sum(100*clock));
end

% randsample), so the fit is reproducible only with a fixed seed
sd = RandStream('mt19937ar','Seed',cfg.seed);
RandStream.setGlobalStream(sd);


% get configuration parameters
seqcat = cfg.seqcat; % sequence category (1|2)
catang = cfg.catang;
seqang = cfg.seqang;
condtn = cfg.condtn;
resp   = cfg.resp;
kappa  = cfg.kappa;

truang = catang+pi/4;
% 
seqlen = cellfun(@length,seqang);
nseq   = length(seqlen); % number of sequences

nsmp = cfg.nsmp; % number of samples
nres = cfg.nres; % number of re-samples for ll bootstrapping
peps = 1e-3;     % infinitesimal response probability

verbose = cfg.verbose; % display level

% define model parameters
pnam = {}; % name
pmin = []; % minimum value
pmax = []; % maximum value
pfun = {}; % log-prior function
pini = []; % initial value
pplb = []; % plausible lower bound
ppub = []; % plausible upper bound
% define general parameters:
% 1/ inference imbalance gamma
pnam{1,1} = 'gamma1';
pmin(1,1) = -2;
pmax(1,1) = 2;
pfun{1,1} = @(x)normpdf(x,0,0.2);
pini(1,1) = normstat(0,0.2);
pplb(1,1) = norminv(0.15,0,0.2);
ppub(1,1) = norminv(0.85,0,0.2);
%
pnam{1,2} = 'gamma2';
pmin(1,2) = -2;
pmax(1,2) = 2;
pfun{1,2} = @(x)normpdf(x,0,0.2);
pini(1,2) = normstat(0,0.2);
pplb(1,2) = norminv(0.15,0,0.2);
ppub(1,2) = norminv(0.85,0,0.2);
%
pnam{1,3} = 'gamma3';
pmin(1,3) = -2;
pmax(1,3) = 2;
pfun{1,3} = @(x)normpdf(x,0,0.2);
pini(1,3) = normstat(0,0.2);
pplb(1,3) = norminv(0.15,0,0.2);
ppub(1,3) = norminv(0.85,0,0.2);
% 2/ inference noise
pnam{1,4} = 'siginf1';
pmin(1,4) = 0.001;
pmax(1,4) = 5; 
pfun{1,4} = @(x)exppdf(x,0.5);
pini(1,4) = expstat(0.5);
pplb(1,4) = expinv(0.15,0.5);
ppub(1,4) = expinv(0.85,0.5);
%
pnam{1,5} = 'siginf2';
pmin(1,5) = 0.001;
pmax(1,5) = 5;
pfun{1,5} = @(x)exppdf(x,0.5);
pini(1,5) = expstat(0.5);
pplb(1,5) = expinv(0.15,0.5);
ppub(1,5) = expinv(0.85,0.5);
%
pnam{1,6} = 'siginf3';
pmin(1,6) = 0.001;
pmax(1,6) = 5;
pfun{1,6} = @(x)exppdf(x,0.5);
pini(1,6) = expstat(0.5);
pplb(1,6) = expinv(0.15,0.5);
ppub(1,6) = expinv(0.85,0.5);

% 3/ selection noise
pnam{1,7} = 'sigsel1';
pmin(1,7) = 0.001;
pmax(1,7) = 5;
pfun{1,7} = @(x)exppdf(x,1);
pini(1,7) = expstat(1);
pplb(1,7) = expinv(0.15,1);
ppub(1,7) = expinv(0.85,1);
%
pnam{1,8} = 'sigsel2';
pmin(1,8) = 0.001;
pmax(1,8) = 5;
pfun{1,8} = @(x)exppdf(x,1);
pini(1,8) = expstat(1);
pplb(1,8) = expinv(0.15,1);
ppub(1,8) = expinv(0.85,1);
%
pnam{1,9} = 'sigsel3';
pmin(1,9) = 0.001;
pmax(1,9) = 5;
pfun{1,9} = @(x)exppdf(x,1);
pini(1,9) = expstat(1);
pplb(1,9) = expinv(0.15,1);
ppub(1,9) = expinv(0.85,1);

% 4/ selection lapses ~ beta(1,49)
pnam{1,10} = 'epsi';
pmin(1,10) = 1e-4; 
pmax(1,10) = 1; 
pfun{1,10} = @(x)betapdf(x,1,49);
pini(1,10) = betastat(1,49);
pplb(1,10) = betainv(0.15,1,49);
ppub(1,10) = betainv(0.85,1,49);

% check fixed parameters
if isfield(cfg,'gamma1')
    cfg.gamma1 = min(max(cfg.gamma1,0),2);
end
if isfield(cfg,'gamma2')
    cfg.gamma2 = min(max(cfg.gamma2,0),2);
end
if isfield(cfg,'gamma3')
    cfg.gamma3 = min(max(cfg.gamma3,0),2);
end
if isfield(cfg,'siginf1')
    cfg.siginf1 = max(cfg.siginf1,0);
end
if isfield(cfg,'siginf2')
    cfg.siginf2 = max(cfg.siginf2,0);
end
if isfield(cfg,'siginf3')
    cfg.siginf3 = max(cfg.siginf3,0);
end
if isfield(cfg,'sigsel1')
    cfg.sigsel1 = max(cfg.sigsel1,0);
end
if isfield(cfg,'sigsel2')
    cfg.sigsel2 = max(cfg.sigsel2,0);
end
if isfield(cfg,'sigsel3')
    cfg.sigsel3 = max(cfg.sigsel3,0);
end
if isfield(cfg,'epsi')
    cfg.epsi = max(cfg.epsi,0);
end

% define fixed parameters
npar = numel(pnam);
pfix = cell(1,npar);
for i = 1:npar
    if isfield(cfg,pnam{i})
        pfix{i} = cfg.(pnam{i});
    end
end

% define free parameters
ifit = cell(1,npar);
pfit_ini = [];
pfit_min = [];
pfit_max = [];
pfit_plb = [];
pfit_pub = [];
n = 1;
for ipar = 1:npar
    if isempty(pfix{ipar}) % free parameter
        ifit{ipar} = n;
        pfit_ini = cat(2,pfit_ini,pini(ipar));
        pfit_min = cat(2,pfit_min,pmin(ipar));
        pfit_max = cat(2,pfit_max,pmax(ipar));
        pfit_plb = cat(2,pfit_plb,pplb(ipar));
        pfit_pub = cat(2,pfit_pub,ppub(ipar));
        n = n+1;
    end
end
nfit = length(pfit_ini);

if nfit > 0
    
    % configure VBMC
    options = vbmc('defaults');
    options.MaxIter = 300; % maximum number of iterations
    options.MaxFunEvals = 500; % maximum number of function evaluations
    options.SpecifyTargetNoise = true; % noisy log-posterior function
    % set display level
    if verbose
        options.Display = 'iter';
    else
        options.Display = 'final';
    end
    
    % fit model using VBMC
    [vp,elbo,~,exitflag,output] = vbmc(@(x)fun(x), ...
        pfit_ini,pfit_min,pfit_max,pfit_plb,pfit_pub,options);
    
    % generate 10^6 samples from the variational posterior
    xsmp = vbmc_rnd(vp,1e6);
    
    xmap = vbmc_mode(vp); % posterior mode
    xavg = mean(xsmp,1); % posterior mean
    xstd = std(xsmp,[],1); % posterior s.d.
    xcov = cov(xsmp); % posterior covariance matrix
    xmed = median(xsmp,1); % posterior medians
    xiqr = quantile(xsmp,[0.25,0.75],1); % posterior interquartile ranges
    
    % create output structure with parameter values
    xhat = xavg; % use posterior average
    phat = getpval(xhat);
    out = cell2struct(phat(:),pnam(:));
    out.llh = getll(phat{:});
    
    % create substructures with parameter values
    phat_map = getpval(xmap); % posterior mode
    out.pmap = cell2struct(phat_map(:),pnam(:));
    phat_avg = getpval(xavg); % posterior mean
    out.pavg = cell2struct(phat_avg(:),pnam(:));

    % store VBMC output
    out.vp   = vp;
    out.elbo = elbo; % ELBO (expected lower bound on log-marginal likelihood)
    out.xnam = pnam(cellfun(@isempty,pfix)); % fitted parameters
    out.xmap = xmap; % posterior mode
    out.xavg = xavg; % posterior means
    out.xstd = xstd; % posterior s.d.
    out.xcov = xcov; % posterior covariance matrix
    out.xmed = xmed; % posterior medians
    out.xiqr = xiqr; % posterior interquartile ranges

    % store extra VBMC output
    out.exitflag = exitflag;
    out.output   = output;
    out.options  = options; 

else
    
    % create output structure with fixed parameter values
    phat = getpval([]);
    out = cell2struct(phat(:),pnam(:));
    
end
out.cfg = cfg;


    function [l,l_sd] = fun(p)
        % get parameter values
        pval = getpval(p);
        % get log-prior
        lp = 0;
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                lp = lp+log(pfun{k}(pval{k}));
            end
        end
        % get log-likelihood
        [ll,ll_sd] = getll(pval{:});
        % get log-posterior statistics
        l = sum(ll)+lp; % estimate
        l_sd = ll_sd; % bootstrap s.d. estimate
    end

    function [ll,ll_sd] = getll(varargin)
        % compute filtered response probability
        p = getp(varargin{:});
        % compute bootstrapped log-likelihood s.d.
        lres = nan(nres,1);
        for ires = 1:nres
            jres = randsample(nsmp,nsmp,true);
            pres = mean(p(:,jres),2);
            lres(ires) = ...
                sum(log(pres(resp==1)))+ ...
                sum(log(1-pres(resp==2)));
        end
        ll_sd = max(std(lres),1e-6);
        % compute log-likelihood
        p = mean(p,2);
        p = peps+(1-peps*2)*p;
        if any(p==0)
            warning('ll will be Inf')
        end
        ll = ...
            sum(log(p(resp==1)))+ ...
            sum(log(1-p(resp==2)));
    end

    function [pval] = getpval(p)
        pval = cell(1,npar); % parameter values
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                pval{k} = p(ifit{k});
            else % fixed parameter
                pval{k} = pfix{k};
            end
        end
    end


    function [p,x] = getp( gamma1,gamma2,gamma3,...
            siginf1,siginf2,siginf3,...
            sigsel1,sigsel2,sigsel3,...
            epsi)
        gamma_vec = [gamma1 gamma2 gamma3];
        siginf_vec = [siginf1 siginf2 siginf3];
        sigsel_vec = [sigsel1 sigsel2 sigsel3];
        
        x = zeros(nseq,nsmp);
        p = zeros(nseq,nsmp);

        for iseq = 1:nseq
            icond = condtn(iseq);
            
            z_seqang = kappa*exp(1i*2*seqang{iseq});
            z_truang = 2*exp(1i*2*truang(iseq));
            tru_real = real(z_truang);
            tru_imag = imag(z_truang);
            
            x_real = real(z_seqang);
            x_imag = imag(z_seqang);
            
            xtmp_real = 0;
            xtmp_imag = 0;
            for itlt = 1:seqlen(iseq)
                % update sequence (with inference imbalance and noise)
                xtmp_real = exp(gamma_vec(icond))*xtmp_real + x_real(itlt) + randn(1,nsmp)*siginf_vec(icond)/2;
                xtmp_imag = exp(gamma_vec(icond))*xtmp_imag + x_imag(itlt) + randn(1,nsmp)*siginf_vec(icond)/2;
            end
            x(iseq,:) = xtmp_real.*tru_real+xtmp_imag.*tru_imag;
            
            p(iseq,:) = normcdf(x(iseq,:),0,repmat(sigsel_vec(icond),1,nsmp));
            p(iseq,:) = (1-epsi)*p(iseq,:)+epsi*0.5;
        end
    end
end
