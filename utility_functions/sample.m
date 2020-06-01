function Asample = sample(A,eval_idx,dim)

%linearly sample matrix A at the indices eval_idx along the dimension 'dim'
%example:
%A = [1 2 3;2 4 6];
%>> sample(A,[1.1 2.5 3],2)
%    1.1    2.5    3.0
%    2.2    5.0    6.0

Andim = numel(size(A));
Adiff = diff(A,1,dim);

%---- extract the integer part of the evaluation indices
eval_idx_floor  = floor(eval_idx);
indices_fl  = repmat({':'},1,Andim); indices_fl{dim}  = eval_idx_floor;
%---- modified integer part to handle special case eval_idx(i_) = size(A,dim)
eval_idx_floor_ = min(eval_idx_floor,size(A,dim)-1);
indices_fl_ = indices_fl;    indices_fl_{dim} = eval_idx_floor_;
%---- extract the decimal part of the evaluation indices
indices_dc = eval_idx - eval_idx_floor;
dimensions_indices_dc = ones(1,Andim); dimensions_indices_dc(dim) = length(eval_idx);

indices_dc_rs = reshape(indices_dc,dimensions_indices_dc);

%linear sampling of A matrix
Asample = A(indices_fl{1:Andim}) + bsxfun(@times,Adiff(indices_fl_{1:Andim}),indices_dc_rs);

end
