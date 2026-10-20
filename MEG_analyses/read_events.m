function evt = read_events(contpath,isubj)
% stimulus events of one participant from SXX_events.tsv
%   contpath  folder of the MEG data (cfg_path.megcont)
%   isubj     participant number
% evt: run, onset (s), sample (200 Hz), trigger_sample (1000 Hz), seqidx, condtn, seqlen, seqpos, seqcat,
%      seqtlt, seqllr, seqang, catang, catmap, keypos (2 x stimulus), repcat, repcor, repkey, reptim

fname = fullfile(contpath,sprintf('S%02d/S%02d_events.tsv',isubj,isubj));
opts = detectImportOptions(fname,'FileType','text','Delimiter','\t');
numvar = setdiff(opts.VariableNames,{'event'});
opts = setvartype(opts,numvar,'double');
opts = setvartype(opts,'event','char');
opts = setvaropts(opts,numvar,'TreatAsMissing','n/a');
T = readtable(fname,opts);

istim = strcmp(T.event,'stimulus');
var2keep = {'run','onset','sample','trigger_sample','seqidx','condtn','seqlen','seqpos','seqcat','seqtlt', ...
    'seqllr','seqang','catang','catmap','repcat','repcor','repkey','reptim'};
evt = struct();
for k = 1:numel(var2keep)
    evt.(var2keep{k}) = T.(var2keep{k})(istim);
end
evt.keypos = [T.keypos1(istim),T.keypos2(istim)]';

end
