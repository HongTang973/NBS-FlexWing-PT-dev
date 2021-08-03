function Gamma = get_Gamma(Gamma_Integration_Function,del_s,E,root_idx,Gamma_root)
Gamma = Gamma_Integration_Function(del_s,E(:,2,:),root_idx) + Gamma_root;
end