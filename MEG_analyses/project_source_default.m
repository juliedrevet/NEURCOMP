function zdef = project_source_default(contpath,isubj,zdec)
% projection of a source map of one participant (e.g. decoding precision of decoding_*_source.m)
% from the participant's cortex to the cortex of the default anatomy (default_cortex.mat), with the
% interpolation matrix computed by Brainstorm (Wdefault in SXX_source_model.mat):
% identical to the projection of the sources on the default anatomy in Brainstorm
%   contpath  folder of the MEG data (cfg_path.megcont)
%   isubj     participant number
%   zdec      vertex x ... (participant's cortex, any number of other dimensions)
% zdef        template vertex x ... (default anatomy)
% [with claude] 

sm = load(fullfile(contpath,sprintf('S%02d/S%02d_source_model.mat',isubj,isubj)),'Wdefault');
sz = size(zdec);
if sz(1) ~= size(sm.Wdefault,2)
    error('not the same number of vertices in the map and in the participant''s cortex!')
end
zdef = reshape(sm.Wdefault*reshape(double(zdec),sz(1),[]),[size(sm.Wdefault,1),sz(2:end)]);

end
