function M_G = map_AtoG(M_A,R_G_A)
M_G = utility_functions.mult_Anm1_Bmpz(R_G_A,M_A);
end

% function M_G = map_AtoG(M_A,R_G_A)
% 
% if 1
%     M_G = static_method_groups.ct_mex('map_AtoG', M_A, R_G_A);
% else
%     M_G = static_method_groups.map_AtoG(M_A, R_G_A);
% end
% end