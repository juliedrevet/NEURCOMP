function collect_decoding_sensor(varnam,freq,subjlist)
% group subject-level files of the sensor-level decoding
%   varnam    decoded variable, as in the file names: 'conllr_particle', 'nulllr_particle', 'absllr' (8 Hz),
%             'seqang' (16 Hz)
%   freq      low-pass frequency of the epochs (Hz)
%   subjlist  participants    [default: 0:36]

if nargin < 3
    subjlist = 0:36;
end

inpath  = './results/zdec_sensor/tmp';
outpath = './results/zdec_sensor';
[~,~] = mkdir(outpath);

ndata = numel(subjlist);
ncond = 3;
outnam = ['zdec_' strrep(varnam,'_particle','')];

% decoding precision: participant x condition (x cos/sin)
grp = [];
for idata = 1:ndata
    for icond = 1:ncond
        fname = fullfile(inpath,sprintf('S%02d_decoding_%s_cond%d_%02dHz.mat',subjlist(idata),varnam,icond,freq));
        load(fname,'zdec','tim_evt');
        nreg = size(zdec,3); % 2 for seqang (cos/sin), 1 otherwise
        if isempty(grp)
            grp = cell(ndata,ncond,nreg);
        end
        for ireg = 1:nreg
            grp{idata,icond,ireg} = zdec(:,:,ireg);
        end
    end
end
out = [];
out.(outnam) = grp;
out.tim_evt = tim_evt;
fname = fullfile(outpath,sprintf('zdec_%s_%02dHz.mat',varnam,freq));
save(fname,'-struct','out')
fprintf('saved in %s\n',fname)

% decoding precision on correct/error (optima/suboptimal) trials: participant x condition (x cos/sin) x opt/sub
if nreg == 1
    grp = cell(ndata,ncond,2);
else
    grp = cell(ndata,ncond,nreg,2);
end
for idata = 1:ndata
    for icond = 1:ncond
        fname = fullfile(inpath,sprintf('S%02d_decoding_%s_cond%d_%02dHz_optsub.mat',subjlist(idata),varnam,icond,freq));
        load(fname,'zdec','tim_evt');
        for iopt = 1:2
            if nreg == 1
                grp{idata,icond,iopt} = zdec(1,:,iopt);
            else
                for ireg = 1:nreg
                    grp{idata,icond,ireg,iopt} = zdec(1,:,ireg,iopt);
                end
            end
        end
    end
end
out = [];
out.(outnam) = grp;
out.tim_evt = tim_evt;
fname = fullfile(outpath,sprintf('zdec_%s_%02dHz_optsub.mat',varnam,freq));
save(fname,'-struct','out')
fprintf('saved in %s\n',fname)

end
