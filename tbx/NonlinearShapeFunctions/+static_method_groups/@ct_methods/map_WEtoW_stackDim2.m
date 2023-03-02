function M_W_flat = map_WEtoW_stackDim2(M_WE_flat,Rs_W_WE_tr,Rv_W_WE_tr,TD)
%mapping performing the transform M_W = R_W_WE*M_WE*R_W_WE_tr
%where M_WE consists of a group of 9x1xns flattened matrices concatinated along dim 2
%rotation R_W_WE must consist of combined 90 degree rotations
%    i.e. all entries of R_W_WE are 1 or -1

% smw = Rs_W_WE_tr; pmw = Rv_W_WE_tr; ns = size(M_WE_flat,3);
% M_W_flat = M_WE_flat;
% for ii_ = 1:3
%     for jj_ = 1:3
%         M_W_flat(ii_+3*(jj_-1),:,:) = smw(ii_)*smw(jj_)*M_WE_flat(pmw(ii_)+3*(pmw(jj_)-1),:,:);
%     end
% end
% 
% if ~isempty(TD)
%     M_W = reshape(M_W_flat,3,[],ns);
%     for I = 1:size(M_W_flat,2)
%         M_W(:,(1:3)+3*(I-1),:) = utility_functions.MultiProd_(M_W(:,(1:3)+3*(I-1),:),TD);
%     end
%     M_W_flat = reshape(M_W,size(M_W_flat));
% end

if 0
    M_W_flat = static_method_groups.mex('map_WEtoW_stackDim2',M_WE_flat,Rs_W_WE_tr,Rv_W_WE_tr,TD);
else
    M_W_flat = static_method_groups.map_WEtoW_stackDim2(M_WE_flat,Rs_W_WE_tr,Rv_W_WE_tr,TD);
end

end
