function collect_decoding_source(varnam,freq,subjlist)
% group subject-level files of the source-level decoding
% decoding precision of each participant is projected from the participant's 
% cortex to the cortex of the default anatomy (as in Brainstorm), all participants can be averaged
%   varnam    decoded variable, as in the file names: 'conllr_particle', 'nulllr_particle', 'absllr' (8 Hz),
%             'seqang' (16 Hz)
%   freq      low-pass frequency of the epochs (Hz)
%   subjlist  participants    [default: 0:36]

if nargin < 3
    subjlist = 0:36;
end

addpath('..'); cfg_path = config_paths;
inpath  = './results/zdec_source/tmp';
outpath = './results/zdec_source';
[~,~] = mkdir(outpath);

ndata = numel(subjlist);
ncond = 3;
sfx = {'','_optsub'};
for k = 1:2
    zdec = [];
    for idata = 1:ndata
        isubj = subjlist(idata);
        fprintf('%s%s: S%02d\n',varnam,sfx{k},isubj)
        for icond = 1:ncond
            fname = fullfile(inpath,sprintf('S%02d_decoding_%s_cond%d_%02dHz%s.mat',isubj,varnam,icond,freq,sfx{k}));
            dat = load(fname,'zdec','tim_evt');
            % template vertex x time sample (x cos/sin) (x opt/sub)
            zdef = single(project_source_default(cfg_path.megcont,isubj,dat.zdec));
            sz = size(zdef);
            if isempty(zdec)
                zdec = nan([sz(1:2),ndata,ncond,sz(3:end)],'single');
                tim_evt = dat.tim_evt;
            elseif max(abs(dat.tim_evt-tim_evt)) > 1e-9
                error('not the same time samples in %s!',fname)
            end
            zdec(:,:,idata,icond,:) = reshape(zdef,[sz(1:2),1,1,prod(sz(3:end))]);
        end
    end
    fname = fullfile(outpath,sprintf('zdec_%s_%02dHz%s.mat',varnam,freq,sfx{k}));
    save(fname,'zdec','tim_evt','subjlist','-v7.3')
    fprintf('saved in %s\n',fname)
end

end
