function [zdec,zdec_optsub] = decode_source_searchlight(dat_real,dat_imag,K,iv_search,C,idxopt,idxsub,nfold,freg)
% ridge decoding at the source level with a searchlight: 
% at each time sample, the epochs are projected on the sources with the inverse 
% kernel, then at each vertex the variable(s) are decoded from the real and 
% imaginary parts of the vertices of its searchlight
%   dat_real, dat_imag  channel x time sample x epoch, real and imaginary parts of the epochs (Hilbert)
%   K                   vertex x channel, inverse kernel (columns in the order of the channels of dat)
%   iv_search           vertex x nsearch, vertices of the searchlight of each vertex
%   C                   epoch x nvar, variables to decode (one decoder per column)
%   idxopt, idxsub      epoch x 1 (logical), optimal and suboptimal trials
%   nfold, freg         number of cross-validation folds and ridge regularization factor
% zdec         vertex x time sample x nvar, decoding precision (Fisher z of the correlation)
% zdec_optsub  vertex x time sample x nvar x 2 (optimal, suboptimal trials)

r2z = @(r)0.5*log((1+r)./(1-r));
nsmp  = size(dat_real,2);
nvert = size(K,1);
nvar  = size(C,2);

% one cell per time sample, so that each worker only receives the data of its time sample
X = cell(1,nsmp);
for ismp = 1:nsmp
    X{ismp} = cat(3,permute(dat_real(:,ismp,:),[1 3 2]),permute(dat_imag(:,ismp,:),[1 3 2])); % channel x epoch x (real,imag)
end

zcell = cell(1,nsmp);
ocell = cell(1,nsmp);
parfor ismp = 1:nsmp
    % project on the sources: epoch x vertex
    src_real = (K*X{ismp}(:,:,1))';
    src_imag = (K*X{ismp}(:,:,2))';
    z = nan(nvert,nvar);
    o = nan(nvert,nvar,2);
    for iv = 1:nvert
        % searchlight: epoch x (real, imag of the nsearch vertices)
        data = [src_real(:,iv_search(iv,:)),src_imag(:,iv_search(iv,:))];
        for ivar = 1:nvar
            out_dec = decode_ridg(C(:,ivar),data,nfold,freg,true);
            z(iv,ivar)   = r2z(out_dec.rdec);
            o(iv,ivar,1) = r2z(corr(out_dec.Ctru(idxopt),out_dec.Chat(idxopt),'type','pearson'));
            o(iv,ivar,2) = r2z(corr(out_dec.Ctru(idxsub),out_dec.Chat(idxsub),'type','pearson'));
        end
    end
    zcell{ismp} = z;
    ocell{ismp} = o;
end
zdec        = permute(cat(3,zcell{:}),[1 3 2]);   % vertex x time sample x nvar
zdec_optsub = permute(cat(4,ocell{:}),[1 4 2 3]); % vertex x time sample x nvar x 2

end
