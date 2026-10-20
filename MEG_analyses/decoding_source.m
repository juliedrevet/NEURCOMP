function decoding_source(isubj,varnam,freq,nblock)
% decoding_source  decoding at the source level, for one participant with a searchlight:
%   * all MEG channels (magnetometers and gradiometers) projected on the participant's cortical sources
%     with the inverse kernel (SXX_source_model.mat)
%   * in each condition, at each vertex and time sample, ridge decoding from the real and
%     imaginary parts of the vertex and its nearest neighbours on the inflated cortex
%   isubj   participant number
%   varnam  'absllr', 'conllr' (conllr and nulllr) or 'seqang'
%   freq    low-pass frequency of the epochs (Hz)    [default: 16 for seqang, 8 otherwise]
%   nblock  number of time samples read at once (e.g. 20 if the memory is short)    [default: all]
% results saved in ./results/zdec_source/tmp:
%   SXX_decoding_<name>_condC_FFHz.mat         zdec: vertex x time sample (x cos/sin for seqang)
%   SXX_decoding_<name>_condC_FFHz_optsub.mat  zdec: vertex x time sample (x cos/sin) x opt/sub

if nargin < 3 || isempty(freq)
    if strcmp(varnam,'seqang')
        freq = 16;
    else
        freq = 8;
    end
end

smp2decode = 1:200;
ntim = numel(smp2decode);
if nargin < 4 || isempty(nblock)
    nblock = ntim;
end

addpath('..'); cfg_path = config_paths;
outpath = './results/zdec_source/tmp';
[~,~] = mkdir(outpath);

nfold   = 10;  % number of cross-validation folds
freg    = 500; % ridge regularization factor
nsearch = 10;  % number of vertices of the searchlight

% variables to decode
dec = get_decoding_variables(cfg_path,isubj,varnam);

epochdir = fullfile(cfg_path.megepo,sprintf('S%02d/epochs_%02dHz_bin',isubj,freq));

% source model: inverse kernel matched to the channels of the epochs, searchlight neighbours
load(fullfile(epochdir,'cfg_cond1.mat'),'cfg_cond');
[K,ichan,iv_search] = load_source_model(cfg_path.megcont,isubj,cfg_cond.label,nsearch);
fprintf('%d vertices, %d channels\n',size(K,1),numel(ichan))

% loop over conditions
condnam = {'PRE','POST','REF'};
for icond = 1:3
    fprintf('*** S%02d %s: decoding %s (%s) ***\n',isubj,varnam,condnam{icond},datestr(now,'HH:MM'))

    % variables to decode, stimuli decoded and optimal/suboptimal trials in this condition
    C      = dec.C(dec.condtn == icond,:);
    ifilt  = dec.ikeep(dec.condtn == icond);
    idxopt = dec.idxopt(dec.condtn == icond);
    idxsub = dec.idxsub(dec.condtn == icond);

    % load epoch configuration structure
    load(fullfile(epochdir,sprintf('cfg_cond%d.mat',icond)),'cfg_cond');
    nepo = numel(cfg_cond.epochs_real);
    if nepo ~= 1440
        fprintf('Warning: not the usual number of epochs for this subject and condition!\n')
    end
    if size(C,1) ~= nepo
        error('not the same number of epochs available!')
    end

    tic
    nreg = size(C,2);
    zdec_all        = nan(size(K,1),ntim,nreg);   % vertex x time sample x regressor
    zdec_optsub_all = nan(size(K,1),ntim,nreg,2); % vertex x time sample x regressor x opt/sub
    for ib = 1:ceil(ntim/nblock)
        iblk = (ib-1)*nblock+1:min(ib*nblock,ntim); % time samples of the block (positions in smp2decode)
        % load epochs: channel x time sample x epoch, all MEG channels of the kernel
        [dat_real,dat_imag] = read_epochs(epochdir,icond,cfg_cond,ichan,smp2decode(iblk),ifilt);

        % decoding at each vertex and time sample
        [zdec_all(:,iblk,:),zdec_optsub_all(:,iblk,:,:)] = decode_source_searchlight(dat_real,dat_imag,K,iv_search,C(ifilt,:),idxopt(ifilt),idxsub(ifilt),nfold,freg);
        clearvars dat_real dat_imag
    end
    toc

    % save
    tim_evt = cfg_cond.tim_evt(smp2decode); % decoded time samples only
    for ifil = 1:numel(dec.filnam)
        zdec = zdec_all(:,:,dec.ireg{ifil});
        fname = fullfile(outpath,sprintf('S%02d_decoding_%s_cond%d_%02dHz.mat',isubj,dec.filnam{ifil},icond,freq));
        save(fname,'zdec','tim_evt')
        fprintf('Decoding matrix saved:\n%s\n',fname)
        zdec = zdec_optsub_all(:,:,dec.ireg{ifil},:);
        if numel(dec.ireg{ifil}) == 1
            zdec = reshape(zdec,size(zdec,1),size(zdec,2),2); % vertex x time sample x opt/sub
        end
        fname = fullfile(outpath,sprintf('S%02d_decoding_%s_cond%d_%02dHz_optsub.mat',isubj,dec.filnam{ifil},icond,freq));
        save(fname,'zdec','tim_evt')
        fprintf('Decoding matrix saved:\n%s\n',fname)
    end
end

end
