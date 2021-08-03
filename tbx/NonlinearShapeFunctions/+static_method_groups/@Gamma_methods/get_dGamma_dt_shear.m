function dGamma_dt = get_dGamma_dt_shear(Gamma_Integration_Function,del_s,E,dE_dt_33ns,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,dGamma_dt_root_G)
if Gamma_approx_lvl == 0
    dGamma_dt_Integrand = utility_functions.MultiProd_(dE_dt_33ns,TauVec)+MultiProd_(E,dTauVec_dt);
elseif Gamma_approx_lvl == 1
    dGamma_dt_Integrand = utility_functions.MultiProd_(dE_dt_33ns,TauVec);
end
dGamma_dt = Gamma_Integration_Function(del_s,dGamma_dt_Integrand,root_idx) + dGamma_dt_root_G;
end