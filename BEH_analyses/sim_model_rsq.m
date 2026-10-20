function  out = sim_model_rsq(cfg)
% simulate model with inference noise on repeated trials with parametrized
% shared variance rsq between noise of the two trial presentations

% check configuration parameters
if ~all(isfield(cfg,{'seqllr','seqlen','siginf','rsq'}))
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

seqllr = cfg.seqllr;
seqlen = cfg.seqlen;
siginf = cfg.siginf;
rsq = cfg.rsq;

nseq  = numel(seqllr);
ncond = numel(siginf);
nsmp  = cfg.nsmp;

r = nan(nseq,nsmp,2,ncond);

x1 = randn(nseq,nsmp); % 1st presentation standard noise
x2 = randn(nseq,nsmp); % independent standard noise
x3 = sqrt(rsq)*x1+sqrt(1-rsq)*x2; % 2nd presentation standard noise: fraction rsq of its variance shared with the 1st

% scale by condition specific SIGINF
for icond = 1:ncond
    y1 = siginf(icond).*x1;
    y2 = siginf(icond).*x3;
    
    % accumulate evidence
    xevi1 = seqllr+sqrt(seqlen).*y1;
    r1 = sign(xevi1);
    
    xevi2 = seqllr+sqrt(seqlen).*y2;
    r2 = sign(xevi2);
    
    r(:,:,:,icond) = cat(3,r1,r2);
end

pcon_rsq = squeeze(mean(mean(r(:,:,1,:)==r(:,:,2,:),1),2));
pcon_rsq = reshape(pcon_rsq,size(siginf));

out = [];
out.rsq = rsq;
out.siginf = siginf;
out.pcon_rsq = pcon_rsq;
out.cfg = cfg;

end