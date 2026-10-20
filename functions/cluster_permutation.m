function stat = cluster_permutation(Y,X,cfg)
% cluster_permutation  cluster-based permutation test on time courses
%   stat = cluster_permutation(Y,X,cfg) tests each time sample and corrects
%   for multiple comparisons across time with a cluster-based permutation
%   test (Maris & Oostenveld, 2007): adjacent samples whose t-value exceeds
%   the cluster-forming threshold are grouped into clusters, and the tstat of
%   each cluster (sum of its t-values) is compared to the null distribution
%   of the maximum absolute cluster tstat across permutations.
%
%   Y     participants x time samples (e.g. decoding time course of each participant)
%   X     participants x 1 variable to correlate with Y (Pearson correlation,
%         permutation of participants), or [] to test Y against 0 (one-sample
%         t-test, random sign flips of participants). For a paired test, pass
%         the difference between the two conditions as Y and X = [].
%   cfg   optional structure:
%         alpha_form  cluster-forming threshold (uncorrected p of each sample)  [0.05]
%         tail        'both', 'right' (positive effects) or 'left'              ['both']
%         nperm       number of permutations                                    [1e4]
%         minsize     minimum number of adjacent samples in a cluster           [2]
%         seed        seed of the random number stream (the global stream is   [1234]
%                     not modified)
%         xvec        time of each sample, for printing                         [1:nsmp]
%         verbose     print the clusters                                        [true]
%         Z           covariate for a partial correlation between X and Y       [none]
%                     (participants x 1, or participants x time samples): the linear effect
%                     of Z (with intercept) is removed from X and Y at each sample, the
%                     residuals of X are permuted, df = ndata-3
%   X can also be participants x time samples: one variable per sample, correlated
%   with the same sample of Y (the participants are permuted the same way for all samples)
%
%   stat  structure:
%         t           1 x nsmp t-values (computed from Pearson r for a correlation)
%         r           1 x nsmp Pearson r (correlation only)
%         df          degrees of freedom
%         tcrit       cluster-forming threshold on |t|
%         clusters    struct array, one element per cluster:
%                       idx   samples of the cluster
%                       tstat sum of the t-values of the cluster
%                       p     cluster-corrected p = (1 + #(null_tstat >= |tstat|))/(nperm + 1)
%                       tbeg  time of the first sample
%                       tend  time of the last sample
%         null_tstat  nperm x 1 maximum absolute cluster tstat of each permutation
%                     (0 when a permutation has no cluster)
% [adapted by claude]

if nargin < 3, cfg = []; end
if nargin < 2, X = []; end
if ~isfield(cfg,'alpha_form'), cfg.alpha_form = 0.05; end
if ~isfield(cfg,'tail'),       cfg.tail = 'both';     end
if ~isfield(cfg,'nperm'),      cfg.nperm = 1e4;       end
if ~isfield(cfg,'minsize'),    cfg.minsize = 2;       end
if ~isfield(cfg,'seed'),       cfg.seed = 1234;       end
if ~isfield(cfg,'verbose'),    cfg.verbose = true;    end

[ndata,nsmp] = size(Y);
if ~isfield(cfg,'xvec'), cfg.xvec = 1:nsmp; end
rs = RandStream('mt19937ar','Seed',cfg.seed); % own stream

iscorr = ~isempty(X);
if iscorr % correlation with X across participants
    if isvector(X), X = X(:); end % participants x 1, or participants x time samples
    if size(X,1) ~= ndata || ~ismember(size(X,2),[1 nsmp])
        error('X must have one value per participant (row of Y), or one column per sample of Y!')
    end
    df = ndata-2;
    if isfield(cfg,'Z') && ~isempty(cfg.Z) % partial correlation: residuals of X and Y after regression on Z
        Z = cfg.Z;
        if isvector(Z), Z = Z(:); end
        X = residuals(X,Z);
        Y = residuals(Y,Z);
        df = ndata-3;
    end
    if size(X,2) == 1 % same X for all samples
        r = corr(X,Y);
        t = r.*sqrt(df./(1-r.^2));
        % permutations: shuffle participants
        Xperm = zeros(ndata,cfg.nperm);
        for iperm = 1:cfg.nperm
            Xperm(:,iperm) = X(randperm(rs,ndata));
        end
        rperm = corr(Xperm,Y); % nperm x nsmp
    else % one X per sample: same permutation of the participants for all samples
        r = corr_paired(X,Y);
        t = r.*sqrt(df./(1-r.^2));
        rperm = zeros(cfg.nperm,nsmp);
        for iperm = 1:cfg.nperm
            rperm(iperm,:) = corr_paired(X(randperm(rs,ndata),:),Y);
        end
    end
    tperm = rperm.*sqrt(df./(1-rperm.^2));
else % one-sample t-test against 0
    df = ndata-1;
    t = mean(Y,1)./(std(Y,[],1)/sqrt(ndata));
    % permutations: flip the sign of participants
    S = 2*(rand(rs,ndata,cfg.nperm) > 0.5)-1; % ndata x nperm
    m = (S'*Y)/ndata; % nperm x nsmp
    v = (sum(Y.^2,1)-ndata*m.^2)/(ndata-1); % squares do not depend on the signs
    tperm = m./sqrt(v/ndata);
end

% cluster-forming threshold
switch cfg.tail
    case 'both'
        tcrit = tinv(1-cfg.alpha_form/2,df);
        sgn = [+1 -1];
    case 'right'
        tcrit = tinv(1-cfg.alpha_form,df);
        sgn = +1;
    case 'left'
        tcrit = tinv(1-cfg.alpha_form,df);
        sgn = -1;
    otherwise
        error('unknown tail ''%s'': use ''both'', ''right'' or ''left''!',cfg.tail)
end

% observed clusters
[cidx,ctstat] = find_clusters(t,tcrit,sgn,cfg.minsize);

% null distribution of the maximum absolute cluster tstat
null_tstat = zeros(cfg.nperm,1);
for iperm = 1:cfg.nperm
    [~,tstat_perm] = find_clusters(tperm(iperm,:),tcrit,sgn,cfg.minsize);
    if ~isempty(tstat_perm)
        null_tstat(iperm) = max(abs(tstat_perm));
    end
end

% cluster-corrected p-values
stat = [];
stat.t = t;
if iscorr
    stat.r = r;
end
stat.df = df;
stat.tcrit = tcrit;
stat.clusters = struct('idx',{},'tstat',{},'p',{},'tbeg',{},'tend',{});
for ic = 1:numel(cidx)
    stat.clusters(ic).idx  = cidx{ic};
    stat.clusters(ic).tstat = ctstat(ic);
    stat.clusters(ic).p    = (1+sum(null_tstat >= abs(ctstat(ic))))/(cfg.nperm+1);
    stat.clusters(ic).tbeg = cfg.xvec(cidx{ic}(1));
    stat.clusters(ic).tend = cfg.xvec(cidx{ic}(end));
end
stat.null_tstat = null_tstat;

if cfg.verbose
    % same printouts as the former code of the figure scripts
    fprintf('%d cluster(s) found\n',numel(cidx))
    for ic = 1:numel(cidx)
        fprintf('\tcluster tstat (sum) : %+5.3f on [%+2.3f %+2.3f] seconds\n',stat.clusters(ic).tstat,stat.clusters(ic).tbeg,stat.clusters(ic).tend);
    end
    if iscorr
        fprintf('<strong>Identifying significance cluster shuffling participants</strong>\n\t')
    else
        fprintf('<strong>Identifying significance cluster flipping the sign of participants</strong>\n\t')
    end
    fprintf('\tplease wait...')
    fprintf('done.\n')
    for ic = 1:numel(cidx)
        fprintf('<strong>Comparing cluster %d/%d :</strong>\n',ic,numel(cidx))
        fprintf('\tcluster-corrected p = %.5f\n',stat.clusters(ic).p)
    end
end
end

function A = residuals(A,Z)
% residuals of the linear regression (with intercept) of each column of A on the
% matching column of Z (A and Z: participants x 1 or participants x time samples)
ncol = max(size(A,2),size(Z,2));
A = A.*ones(1,ncol);
Z = Z.*ones(1,ncol);
A = A-mean(A,1);
Z = Z-mean(Z,1);
A = A-Z.*(sum(Z.*A,1)./sum(Z.^2,1));
end

function r = corr_paired(X,Y)
% Pearson correlation between each column of X and the same column of Y
X = X-mean(X,1);
Y = Y-mean(Y,1);
r = sum(X.*Y,1)./sqrt(sum(X.^2,1).*sum(Y.^2,1));
end

function [cidx,ctstat] = find_clusters(t,tcrit,sgn,minsize)
% clusters of at least minsize adjacent samples with sgn*t > tcrit, for each sign in sgn
cidx = {};
ctstat = [];
for s = sgn
    above = [false, s*t > tcrit, false];
    ibeg = find(diff(above) == 1);
    iend = find(diff(above) == -1)-1;
    for k = 1:numel(ibeg)
        if iend(k)-ibeg(k)+1 >= minsize
            cidx{end+1} = ibeg(k):iend(k); 
            ctstat(end+1) = sum(t(ibeg(k):iend(k))); 
        end
    end
end
end
