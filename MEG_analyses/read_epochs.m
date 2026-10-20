function [dat_real,dat_imag] = read_epochs(epochdir,icond,cfg_cond,ichan,ismp,iepo)
% read_epochs  real and imaginary parts (Hilbert) of the epochs of one condition, from the binary files of
% epoching_meg_from_continuous.m: channel x time sample x epoch, for the channels ichan, the time samples
% ismp and the epochs iepo
%   epochdir  folder of the epochs (cfg_path.megepo/SXX/epochs_FFHz_bin)
%   cfg_cond  epoch configuration structure of the condition (cfg_condC.mat)
%   ismp      ':' or all the time samples: whole files read at once (fastest)
%             block of consecutive time samples: only this block is read, one read per epoch (less memory)

nchan = numel(cfg_cond.label);
nsmp  = numel(cfg_cond.tim_evt);
nepo  = numel(cfg_cond.epochs_real);
if numel(cfg_cond.epochs_imag) ~= nepo
    error('missing real or imag file!')
end
switch cfg_cond.type
    case {'double','float64'}
        nbyte = 8;
    case {'single','float32'}
        nbyte = 4;
    otherwise
        error('unknown data type %s!',cfg_cond.type)
end
whole = ischar(ismp) || isequal(ismp(:)',1:nsmp);
if ~whole && any(diff(ismp) ~= 1)
    error('the time samples of a block must be consecutive!')
end
if islogical(iepo)
    iepo = find(iepo);
end

part = {'real','imag'};
dat = cell(1,2);
for ipart = 1:2
    d = dir(fullfile(epochdir,sprintf('*cond%d_hilbert_%s.bin',icond,part{ipart})));
    if d.bytes ~= nchan*nsmp*nepo*nbyte
        error('unexpected size of binary file!')
    end
    fileID = fopen(fullfile(d.folder,d.name));
    if whole
        % whole binary file at once (channel x time sample x epoch)
        A = fread(fileID,[nchan nsmp*nepo],cfg_cond.type);
        A = reshape(A,nchan,nsmp,nepo);
        dat{ipart} = A(ichan,ismp,iepo);
        clearvars A
    else
        % block of time samples: for each epoch, the block is one piece of the file (channel x time sample)
        B = nan(numel(ichan),numel(ismp),numel(iepo));
        for k = 1:numel(iepo)
            fseek(fileID,((iepo(k)-1)*nsmp+ismp(1)-1)*nchan*nbyte,'bof');
            a = fread(fileID,[nchan numel(ismp)],cfg_cond.type);
            B(:,:,k) = a(ichan,:);
        end
        dat{ipart} = B;
    end
    fclose(fileID);
end
[dat_real,dat_imag] = dat{:};

end
