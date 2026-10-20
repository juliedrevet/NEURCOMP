function dec = get_decoding_variables(cfg_path,isubj,varnam)
% get the variable to be decoded from the MEG data of one participant
%   varnam = 'absllr'  evidence strength |LLR|
%   varnam = 'conllr'  belief update along the true axis (conllr_particle) and along the null axis
%                      (nulllr_particle), from the same particles
%   varnam = 'seqang'  stimulus orientation, as cos(2*seqang) and sin(2*seqang)
% dec.C        stimulus x regressor
% dec.ikeep    stimulus x 1, stimuli decoded
% dec.filnam   names of the output files (SXX_decoding_<filnam>_condC_FFHz.mat)
% dec.ireg     regressors (columns of C) saved in each output file
% dec.condtn   condition of each stimulus
% dec.idxopt   optimal trials (response consistent with the sign of the sequence LLR)
% dec.idxsub   suboptimal trials
% [with claude]

kappa = 0.5;
npart = 1e4; % number of particles (conllr)

evt = read_events(cfg_path.megcont,isubj);
% LLR w.r.t. true categories
angllr = 2*kappa*sin(2*(evt.seqang-evt.catang));

switch varnam
    case 'absllr'
        C      = abs(angllr);
        ikeep  = true(size(C));
        filnam = {'absllr'};
        ireg   = {1};
    case 'conllr'
        % LLR w.r.t. null categories (category boundaries)
        angllr_null = 2*kappa*sin(2*(evt.seqang-(evt.catang+pi/4)));
        % best-fitting parameters of the participant
        load(fullfile(cfg_path.results,'results_fit.mat'),'out_fit');
        out_fit = out_fit{isubj+1};
        gamma  = out_fit.xavg(contains(out_fit.xnam,'gamma'));  % per condition
        siginf = out_fit.xavg(contains(out_fit.xnam,'siginf')); % per condition
        epsi   = out_fit.xavg(strcmp(out_fit.xnam,'epsi'));
        % belief before each stimulus on each axis (particle filter) and belief updates
        rng(100*isubj,'twister'); % same particles at each run, at the sensor and source levels
        fprintf('Simulating accumulated LLR (particle filter)...')
        [belief,belief_null] = particle_belief(angllr,angllr_null,evt.seqpos,evt.seqlen,evt.condtn,evt.repcat,gamma,siginf,epsi,npart);
        fprintf('done.\n')
        C      = [sign(belief).*angllr,sign(belief_null).*angllr_null]; % NaN for the first stimulus of each sequence
        ikeep  = evt.seqpos ~= 1;
        filnam = {'conllr_particle','nulllr_particle'};
        ireg   = {1,2};
    case 'seqang'
        C      = [cos(2*evt.seqang),sin(2*evt.seqang)];
        ikeep  = true(size(evt.seqang));
        filnam = {'seqang'};
        ireg   = {[1 2]};
    otherwise
        error('unknown variable %s!',varnam)
end

dec = [];
dec.C      = C;
dec.ikeep  = ikeep;
dec.filnam = filnam;
dec.ireg   = ireg;
dec.condtn = evt.condtn;
dec.idxopt = sign(1.5-evt.repcat) == sign(evt.seqllr);
dec.idxsub = sign(1.5-evt.repcat) ~= sign(evt.seqllr);

end

function [belief,belief_null] = particle_belief(llr,llr_null,seqpos,seqlen,condtn,repcat,gamma,siginf,epsi,npart)
% mean belief before each stimulus on the true and null axes (NaN for the
% first stimulus of each sequence): npart noisy accumulations with inference imbalance per
% sequence of the 2D stimulus vectors (true LLR, null LLR) with isotropic
% inference noise (independent on each coordinate), resampled according to
% the likelihood of the participant's choice given the evidence accumulated
% on the true axis after the last stimulus (with lapses, no selection noise)
belief      = nan(size(llr));
belief_null = nan(size(llr));
istart = find(seqpos == 1);
for iseq = 1:numel(istart)
    iang = istart(iseq)+(0:seqlen(istart(iseq))-1);
    icond = condtn(iang(1));
    x  = zeros(numel(iang),npart); % evidence accumulated after each stimulus - true axis
    xn = zeros(numel(iang),npart); % evidence accumulated after each stimulus - null axis
    xcur  = zeros(1,npart);
    xncur = zeros(1,npart);
    for k = 1:numel(iang)
        xcur  = exp(gamma(icond))*xcur  + llr(iang(k))      + siginf(icond)*randn(1,npart);
        xncur = exp(gamma(icond))*xncur + llr_null(iang(k)) + siginf(icond)*randn(1,npart);
        x(k,:)  = xcur;
        xn(k,:) = xncur;
    end
    % likelihood of the participant's choice (+1 = response 1, -1 = response 2)
    rtru = sign(1.5-repcat(iang(end)));
    w = epsi/2 + (1-epsi)*(sign(x(end,:)) == rtru);
    ikeep = randsample(npart,npart,true,w/sum(w));
    belief(iang(2:end))      = mean(x(1:end-1,ikeep),2);
    belief_null(iang(2:end)) = mean(xn(1:end-1,ikeep),2);
end
end
