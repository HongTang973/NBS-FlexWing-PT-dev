function M_W = map_WEtoW(M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD)
%M_W = MultiProd_(R_W_WE,MultiProd_(M_WE,R_W_WE_tr,[1 2]),[1 2]);%(3)x(3)x(ns)

%if the rotation R_W_WE is modulo 90 degrees use following code
%much quicker than MultiProd_ line
% smw = Rs_W_WE_tr; pmw = Rv_W_WE_tr;
% M_W = M_WE;
% for ii_ = 1:3
%     for jj_ = 1:3
%         M_W(ii_,jj_,:) = smw(ii_)*smw(jj_)*M_WE(pmw(ii_),pmw(jj_),:);
%     end
% end
% if ~isempty(TD)
%     M_W = utility_functions.MultiProd_(M_W,TD);
% end

if 0
    M_W = static_method_groups.mex('map_WEtoW', M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD);

else
    M_W = static_method_groups.map_WEtoW(M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD);
end
end