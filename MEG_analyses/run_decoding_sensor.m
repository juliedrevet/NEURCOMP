%% run_decoding_sensor
% sensor-level decoding of all participants (decoding_sensor.m)
% Results: results/zdec_sensor/tmp/SXX_decoding_<variable>_condC_FFHz.mat and _optsub.mat; then group files
% (collect_decoding_sensor.m): results/zdec_sensor/zdec_<variable>_FFHz.mat

clearvars
close all
clc

addpath('..'); cfg_path = config_paths;

subjlist = 0:36;
% parpool(8); % fewer workers if the memory is short (default: one worker per core)
nblock = []; % time samples read at once: [] = all (fastest), e.g. 20 if the memory is short

% check that the epochs of all participants exist (16 Hz: seqang, 8 Hz: absllr and conllr)
missing = {};
for isubj = subjlist
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
    error('%d epoch files missing (listed above): run run_epoching.m first!',numel(missing))
end

%% stimulus orientation, 16 Hz
for isubj = subjlist
    fprintf('\n##### S%02d: seqang (%s) #####\n',isubj,datestr(now))
    decoding_sensor(isubj,'seqang',16,nblock)
end
% group files, once all participants are done
collect_decoding_sensor('seqang',16,0:36)

%% evidence strength, 8 Hz
for isubj = subjlist
    fprintf('\n##### S%02d: absllr (%s) #####\n',isubj,datestr(now))
    decoding_sensor(isubj,'absllr',8,nblock)
end
collect_decoding_sensor('absllr',8,0:36)

%% belief updates along the true axis (conllr) and the null axis (nulllr), particle filter, 8 Hz
for isubj = subjlist
    fprintf('\n##### S%02d: conllr and nulllr (%s) #####\n',isubj,datestr(now))
    decoding_sensor(isubj,'conllr',8,nblock)
end
collect_decoding_sensor('conllr_particle',8,0:36)
collect_decoding_sensor('nulllr_particle',8,0:36)
