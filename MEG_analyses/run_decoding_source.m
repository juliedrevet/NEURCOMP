%% run_decoding_source
% source-level decoding of all participants (decoding_source.m) LONG!
% Results: results/zdec_source/tmp/SXX_decoding_<variable>_condC_FFHz.mat (zdec: vertex x
% time sample, on the participant's cortex) and _optsub.mat; then group files projected on the default
% anatomy (collect_decoding_source.m, once all participants are done): results/zdec_source/zdec_<variable>_FFHz.mat
% and _optsub.mat
% Time (12 workers): about 55 min per participant for conllr and nulllr

clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

subjlist = 0:6;
% parpool(8); % fewer workers if the memory is short (default: one worker per core)
nblock = []; % time samples read at once: [] = all (fastest), e.g. 20 if the memory is short

% check that the source models and the epochs of all participants exist (16 Hz: seqang, 8 Hz: absllr and conllr)
missing = {};
for isubj = subjlist
    fname = fullfile(cfg_path.megcont,sprintf('S%02d/S%02d_source_model.mat',isubj,isubj));
    if ~exist(fname,'file')
        missing{end+1} = fname; %#ok<SAGROW>
    end
    for freq = [16 8]
        epochdir = fullfile(cfg_path.megepo,sprintf('S%02d/epochs_%02dHz_bin',isubj,freq));
        for icond = 1:3
            fpat = {sprintf('cfg_cond%d.mat',icond),sprintf('*cond%d_hilbert_real.bin',icond),sprintf('*cond%d_hilbert_imag.bin',icond)};
            for i = 1:numel(fpat)
                if isempty(dir(fullfile(epochdir,fpat{i})))
                    missing{end+1} = fullfile(epochdir,fpat{i});
                end
            end
        end
    end
end
if ~isempty(missing)
    fprintf('%s\n',missing{:})
    error('%d files missing (listed above): source models from Zenodo, epochs from run_epoching.m!',numel(missing))
end

%% belief updates along the true axis (conllr) and the null axis (nulllr), particle filter, 8 Hz
for isubj = subjlist
    fprintf('\n##### S%02d: conllr and nulllr (%s) #####\n',isubj,datestr(now))
    decoding_source(isubj,'conllr',8,nblock)
end
% group files on the default anatomy, once all participants are done
collect_decoding_source('conllr_particle',8,0:36)
collect_decoding_source('nulllr_particle',8,0:36)

%% stimulus orientation, 16 Hz
for isubj = subjlist
    fprintf('\n##### S%02d: seqang (%s) #####\n',isubj,datestr(now))
    decoding_source(isubj,'seqang',16,nblock)
end
collect_decoding_source('seqang',16,0:36)

%% evidence strength, 8 Hz
for isubj = subjlist
    fprintf('\n##### S%02d: absllr (%s) #####\n',isubj,datestr(now))
    decoding_source(isubj,'absllr',8,nblock)
end
collect_decoding_source('absllr',8,0:36)
