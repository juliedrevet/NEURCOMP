function [K,ichan,iv_search] = load_source_model(contpath,isubj,label,nsearch)
% source model of one participant for the source-level decoding:
% inverse kernel and searchlight neighbours, from SXX_source_model.mat
%   contpath     folder of the MEG data (cfg_path.megcont)
%   isubj        participant number
%   label        channel labels of the epochs (cfg_cond.label)
%   nsearch      number of vertices of each searchlight
% K          vertex x channel inverse kernel, columns in the order of label(ichan)
% ichan      channels of the epochs used (those of the kernel)
% iv_search  vertex x nsearch, vertices of the searchlight of each vertex
% [with claude]

sm = load(fullfile(contpath,sprintf('S%02d/S%02d_source_model.mat',isubj,isubj)),'ImagingKernel','label','Vinflated');

% kernel columns matched to the channels of the epochs by name
[isinkern,ikern] = ismember(label,sm.label);
if ~all(ismember(sm.label,label))
    error('channels of the kernel missing in the epochs!')
end
ichan = find(isinkern);
K = sm.ImagingKernel(:,ikern(isinkern));

% searchlight: nearest vertices on the inflated cortex (the vertex itself included)
V = sm.Vinflated;
nvert = size(V,1);
if size(K,1) ~= nvert
    error('not the same number of vertices in the kernel and in the inflated cortex!')
end
iv_search = zeros(nvert,nsearch);
for iv = 1:nvert
    dist = sqrt(sum(bsxfun(@minus,V,V(iv,:)).^2,2));
    [~,ivnew] = mink(dist,nsearch);
    iv_search(iv,:) = ivnew;
end

end
