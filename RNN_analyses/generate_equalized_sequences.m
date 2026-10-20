function out = generate_equalized_sequences(nseq,seqlen,sigcat,seed)
% generate_equalized_sequences  random sequences of the task, with the same number of sequences per length
%   out = generate_equalized_sequences(nseq,seqlen,sigcat,seed)
%   nseq    number of sequences per sequence length
%   seqlen  sequence lengths
%   sigcat  s.d. (deg) of the noise on the category axis retrieved at each stimulus
%   seed    seed of the random number generator (rng(seed,'twister'))
%   out     nseq*numel(seqlen) sequences: stimulus orientations (seqang), generative category (seqcat),
%           category axis (catang), evidence of each stimulus with the true (tltllr) or the noisy retrieved
%           (llrhat) category axis, evidence relative to the null categories (tltnul, nulhat), true and noisy
%           retrieved reference (reftru, refhat), evidence of the sequence (seqllr)

rng(seed,'twister');

nseqlen = numel(seqlen);

N = nseq*nseqlen;

out = [];
out.seqang = cell(N,1);
out.seqcat = nan(N,1);
out.catang = nan(N,1);
out.reftru = nan(N,1);
out.refhat = cell(N,1);
out.tltllr = cell(N,1);
out.seqllr = nan(N,1);
out.tltnul = cell(N,1);
out.llrhat = cell(N,1);
out.nulhat = cell(N,1);
out.seqtlt = cell(N,1);

rndvm = @(n,t,k)mod(vonmises_rnd(t*2,k,n)/2,pi);
kappa = 0.5;   % concentration of the von Mises distributions

sigcat_rad  = sigcat * pi / 180;

idx = 1;
% generate the sequences of each length
for iseqlen =  1:nseqlen
    ncue = seqlen(iseqlen);

    seqtlt  = rndvm([nseq, ncue], 0, kappa);
    catang  = rand(nseq,1) * pi;
    seqcat  = randi(2,nseq,1);
    reftru  = mod(catang - pi/4, pi);
    refhat  = reftru + randn(nseq,ncue) * sigcat_rad;

    offset = mod(catang+(seqcat-1)*pi/2,pi); % generative category depending on seqcat
    seqang = mod(seqtlt+offset,pi);

    tltllr = 2*kappa*cos(2*seqtlt);
    tltllr = tltllr.*(seqcat==1)-tltllr.*(seqcat==2); % sign depending on the generative category: ORANGE > 0; BLUE < 0
    tltnul = 2*kappa*(cos(2*(seqtlt-pi/4)).*(seqcat==1)+cos(2*(seqtlt+pi/4)).*(seqcat==2));
    seqllr = sum(tltllr,2);

    llrhat = 2*kappa*sin(2*(seqang-refhat));
    nulhat = -2*kappa*cos(2*(seqang-refhat));
    % fill preallocated outputs
    for i = 1:nseq
        out.seqang{idx} = seqang(i, :);
        out.seqcat(idx) = seqcat(i);
        out.reftru(idx) = reftru(i);
        out.refhat{idx} = refhat(i, :);
        out.seqllr(idx) = seqllr(i);
        out.seqtlt{idx} = seqtlt(i, :);
        out.tltllr{idx} = tltllr(i, :);
        out.catang(idx) = catang(i);
        out.tltnul{idx} = tltnul(i, :);
        out.llrhat{idx} = llrhat(i, :);
        out.nulhat{idx} = nulhat(i, :);
        idx = idx + 1;
    end
end

end
