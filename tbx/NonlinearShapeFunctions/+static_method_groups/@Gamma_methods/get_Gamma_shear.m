function Gamma = get_Gamma_shear(Gamma_Integration_Function,del_s,E,TauVec,root_idx,Gamma_root)
Gamma_Integrand = utility_functions.MultiProd_(E,TauVec);
Gamma = Gamma_Integration_Function(del_s,Gamma_Integrand,root_idx) + Gamma_root;
end