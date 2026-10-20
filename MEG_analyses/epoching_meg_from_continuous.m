% epoching_meg_from_continuous
% epoching from the preprocessed continuous data (shared data), instead of the raw .fif files
% (previous epoching script): same epochs, same binary files, read as before by the decoding scripts.
%
% Input: shared data of the participant in cfg_path.megcont/SXX (prepare/export_meg_continuous.m)
%   SXX_meg_continuous.mat  dat (1 x 6 cell, runs: channel x time sample, 200 Hz, first sample at t = 0),
%                           fsample, label: MEG data already preprocessed (blink and heart beat components
%                           removed, high-pass 0.5 Hz, low-pass 40 Hz, downsampled to 200 Hz)
%   SXX_events.tsv          events (read_events.m): run, onset (trigger-video offset corrected)
%                           and condition of each stimulus
% Output (as before): one BIN file per condition (real and imag) in cfg_path.megepo/SXX/epochs_%02dHz_bin
% and cfg_cond1/2/3.mat; each epoch:
%   * type: double
%   * 200 Hz (200 time steps between -200 and +800ms)
%   * low-pass filter (freq), linear trend removed and Hilbert transform on the epoch with a padding
%     of 200 ms on each side (-400 to +1000 ms), then padding removed
%   * stimulus onset on the nearest sample (200 Hz grid: up to 2.5 ms of difference)

function epoching_meg_from_continuous(isubj,freq)
%   isubj  participant number
%   freq   low-pass frequency (Hz)    [default: 8]

if nargin < 2
    freq = 8;
end

clearvars -except isubj freq
close all
clc

addpath('..'); cfg_path = config_paths;
addpath(cfg_path.fieldtrip);
ft_defaults

contpath = cfg_path.megcont; % shared data

% preprocessed continuous data and stimulus events of the participant
fname_cont = fullfile(contpath,sprintf('S%02d/S%02d_meg_continuous.mat',isubj,isubj));
fprintf('Loading %s...',fname_cont)
cont = load(fname_cont,'dat','fsample','label');
fprintf('done.\n')
fsample = cont.fsample;
evt = read_events(contpath,isubj);

% preprocessing configuration (in seconds, converted to samples at fsample)
twin     = [-0.200,+0.800]; % time window of interest (s)
tpad     = 0.200; % padding window (s)
flowpass = freq; % low-pass frequency (Hz)
iwin     = round(twin*fsample); % time window of interest (samples)
ipad     = round(tpad*fsample); % padding window (samples)
tim_evt  = (iwin(1):(iwin(2)-1))/fsample;

epochdir = fullfile(cfg_path.megepo,sprintf('S%02d/epochs_%02dHz_bin',isubj,freq));

% add epoch folder
if ~exist(epochdir,'dir')
    mkdir(epochdir);
end

cfg_dat = [];
cfg_dat.participant = isubj;
cfg_dat.flowpass = flowpass;
cfg_dat.contfilename = fname_cont;
cfg_dat.contpath = contpath;
cfg_dat.fsample = fsample;
cfg_dat.iwin = iwin;
cfg_dat.ipad = ipad;
cfg_dat.type = 'double';
cfg_dat.label = cont.label;
cfg_dat.process = {'continuous data: components removed, high-pass 0.5 Hz, low-pass 40 Hz, 200 Hz'};

fname_real = cell(1,3);
fname_imag = cell(1,3);
epochs_real = cell(1,3); % list of epochs per condition
epochs_imag = cell(1,3);
for icond = 1:3
    fname_real{icond} = fullfile(epochdir,sprintf('S%02d_cond%d_hilbert_real.bin',isubj,icond));
    fname_imag{icond} = fullfile(epochdir,sprintf('S%02d_cond%d_hilbert_imag.bin',isubj,icond));
    % create (!overwrites!) binary files
    fileID = fopen(fname_real{icond},'w'); fclose(fileID);
    fileID = fopen(fname_imag{icond},'w'); fclose(fileID);
end

for irun = 1:6
    dat_run = double(cont.dat{irun});
    irow    = find(evt.run == irun); % stimuli of the run
    onset   = evt.onset(irow);
    condtn  = evt.condtn(irow);
    nevt    = numel(irow);

    fprintf('\nEpoching run %d...\n',irun)
    for ievt = 1:nevt
        if mod(ievt,30)==0
            fprintf('\t...%03d/%03d\n',ievt,nevt)
        end
        % stimulus onset: nearest sample (first sample of the run at t = 0)
        i0 = round(onset(ievt)*fsample)+1;
        % extract data segment (with padding)
        dat_tmp = dat_run(:,i0+(iwin(1)-ipad:iwin(2)-1+ipad));
        % low-pass filter data segment
        dat_tmp = ft_preproc_lowpassfilter(dat_tmp,fsample,flowpass,6,'but');
        % remove linear trend in the data
        dat_tmp = ft_preproc_polyremoval(dat_tmp,1);
        % apply Hilbert transform to extract real and imaginary components
        dat_evt_real = ft_preproc_hilbert(dat_tmp,'real');
        dat_evt_imag = ft_preproc_hilbert(dat_tmp,'imag');
        % remove padding
        dat_evt_real = dat_evt_real(:,1+ipad:end-ipad);
        dat_evt_imag = dat_evt_imag(:,1+ipad:end-ipad);

        % append to the binary files of the condition
        icond = condtn(ievt);
        fileID = fopen(fname_real{icond},'a');
        fwrite(fileID,dat_evt_real,cfg_dat.type);
        fclose(fileID);
        fileID = fopen(fname_imag{icond},'a');
        fwrite(fileID,dat_evt_imag,cfg_dat.type);
        fclose(fileID);
        epochs_real{icond} = [epochs_real{icond};{sprintf('S%02d_run%02d_epoch%03d_real_cond%d',isubj,irun,ievt,icond)}];
        epochs_imag{icond} = [epochs_imag{icond};{sprintf('S%02d_run%02d_epoch%03d_imag_cond%d',isubj,irun,ievt,icond)}];
    end
    fprintf('done.\n')
end
cfg_dat.process = [cfg_dat.process; ...
    {sprintf('epoch %+d to %+d ms with padding %d ms',round(1000*twin(1)),round(1000*twin(2)),round(1000*tpad))}; ...
    {'ft_preproc_lowpassfilter 6 but'};{'ft_preproc_polyremoval 1st order'};{'ft_preproc_hilbert real and imag'};{'unpad'}];
cfg_dat.tim_evt = tim_evt;

% save configuration structure per condition
for icond = 1:3
    cfg_cond = cfg_dat;
    cfg_cond.epochs_real = epochs_real{icond};
    cfg_cond.epochs_imag = epochs_imag{icond};
    save(fullfile(epochdir,sprintf('cfg_cond%d.mat',icond)),'cfg_cond')
end

end
