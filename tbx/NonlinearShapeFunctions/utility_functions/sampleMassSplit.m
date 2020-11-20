function massValues_new = sampleMassSplit(massValues,x_old,x_new,varargin)

%Discrete Mass Type Sampling
%Values of the old discribution are split between neighboring sample points of the new disctribution.

dim = get_option(varargin,'dim',1);

permuteVec = [1 2 3];
permuteVec(1) = dim;
permuteVec(dim) = 1;

smat_discrete = sampleMat(x_new,x_old);

massValues_permute = permute(massValues,permuteVec);
sz = size(massValues_permute);
massValues_permute_reshape = reshape(massValues_permute,sz(1),[],1);
massValues_new_permute_reshape = smat_discrete*massValues_permute_reshape;
if issparse(massValues_new_permute_reshape), massValues_new_permute_reshape = full(massValues_new_permute_reshape); end
massValues_new_permute = reshape(massValues_new_permute_reshape,[numel(x_new),sz(2:end)]);
massValues_new = permute(massValues_new_permute,permuteVec);

end