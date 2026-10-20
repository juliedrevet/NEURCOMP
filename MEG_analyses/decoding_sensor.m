function decoding_sensor(isubj,varnam,freq,nblock)
% decoding_sensor  decoding at the sensor level, for one participant: in each condition and at each time
% sample, ridge decoding of the variable(s) from the real and imaginary parts of the planar gradiometers
%   isubj   participant number
%   varnam  'absllr', 'conllr' (conllr and nulllr) or 'seqang' (get_decoding_variables.m)
%   freq    low-pass frequency of the epochs (Hz)    [default: 16 for seqang, 8 otherwise]
%   nblock  number of time samples read at once (e.g. 20 if the memory is short)    [default: all]
% results saved in ./results/zdec_sensor/tmp:
%   SXX_decoding_<name>_condC_FFHz.mat         zdec: 1 x time sample (x cos/sin for seqang)
%   SXX_decoding_<name>_condC_FFHz_optsub.mat  zdec: 1 x time sample (x cos/sin) x opt/sub
% [with claude]

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
outpath = './results/zdec_sensor/tmp';
[~,~] = mkdir(outpath);

nfold = 10;  % number of cross-validation folds
freg  = 500; % ridge regularization factor
% define Fisher transform for correlation coefficient
r2z = @(r)0.5*log((1+r)./(1-r));

% variables to decode
dec  = get_decoding_variables(cfg_path,isubj,varnam);
nreg = size(dec.C,2);

epochdir = fullfile(cfg_path.megepo,sprintf('S%02d/epochs_%02dHz_bin',isubj,freq));

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

    % planar gradiometers (odd:MEG*2 even:MEG*3)
    channum = cellfun(@(s)str2num(s(4:end)),cfg_cond.label);
    igrd = find(mod(channum,10) ~= 1);

    tic
    C      = C(ifilt,:);
    idxopt = idxopt(ifilt);
    idxsub = idxsub(ifilt);
    zcell = cell(1,ntim);
    ocell = cell(1,ntim);
    for ib = 1:ceil(ntim/nblock)
        iblk = (ib-1)*nblock+1:min(ib*nblock,ntim); % time samples of the block (positions in smp2decode)
        % load epochs: epoch x channel x time sample
        [dat_real,dat_imag] = read_epochs(epochdir,icond,cfg_cond,igrd,smp2decode(iblk),ifilt);
        dat_real = permute(dat_real,[3 1 2]);
        dat_imag = permute(dat_imag,[3 1 2]);

        % decoding at each time sample
        zblk = cell(1,numel(iblk));
        oblk = cell(1,numel(iblk));
        parfor k = 1:numel(iblk)
            data = [dat_real(:,:,k),dat_imag(:,:,k)]; % epoch x (grad real, grad imag)
            z = nan(1,nreg);
            o = nan(1,nreg,2);
            for ireg = 1:nreg
                out_dec = decode_ridg(C(:,ireg),data,nfold,freg,true);
                z(ireg) = r2z(out_dec.rdec);
                % correlation on optimal/suboptimal trials
                o(1,ireg,1) = r2z(corr(out_dec.Ctru(idxopt),out_dec.Chat(idxopt),'type','pearson'));
                o(1,ireg,2) = r2z(corr(out_dec.Ctru(idxsub),out_dec.Chat(idxsub),'type','pearson'));
            end
            zblk{k} = z;
            oblk{k} = o;
        end
        zcell(iblk) = zblk;
        ocell(iblk) = oblk;
        clearvars dat_real dat_imag
    end
    toc
    zdec_all        = reshape(cat(1,zcell{:}),[1 ntim nreg]);   % 1 x time sample x regressor
    zdec_optsub_all = reshape(cat(1,ocell{:}),[1 ntim nreg 2]); % 1 x time sample x regressor x opt/sub

    % save
    tim_evt = cfg_cond.tim_evt(smp2decode);
    for ifil = 1:numel(dec.filnam)
        zdec = zdec_all(:,:,dec.ireg{ifil});
        fname = fullfile(outpath,sprintf('S%02d_decoding_%s_cond%d_%02dHz.mat',isubj,dec.filnam{ifil},icond,freq));
        save(fname,'zdec','tim_evt')
        fprintf('Decoding matrix saved:\n%s\n',fname)
        zdec = zdec_optsub_all(:,:,dec.ireg{ifil},:);
        if numel(dec.ireg{ifil}) == 1
            zdec = reshape(zdec,size(zdec,1),size(zdec,2),2); % 1 x time sample x opt/sub
        end
        fname = fullfile(outpath,sprintf('S%02d_decoding_%s_cond%d_%02dHz_optsub.mat',isubj,dec.filnam{ifil},icond,freq));
        save(fname,'zdec','tim_evt')
        fprintf('Decoding matrix saved:\n%s\n',fname)
    end
end

end
