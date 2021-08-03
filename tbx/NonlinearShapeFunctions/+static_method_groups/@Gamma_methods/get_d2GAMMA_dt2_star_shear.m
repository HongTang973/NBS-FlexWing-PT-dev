function d2Gamma_dt2_star = get_d2GAMMA_dt2_star_shear(Gamma_Integration_Function,del_s,dE_dt,d2E_dt2_star,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,d2Gamma_dt2_root_G_star)
if Gamma_approx_lvl == 0
    d2Gamma_dt2_star_Integrand = utility_functions.MultiProd_(d2E_dt2_star,TauVec)+2*MultiProd_(dE_dt,dTauVec_dt);
elseif Gamma_approx_lvl == 1
    d2Gamma_dt2_star_Integrand = utility_functions.MultiProd_(d2E_dt2_star,TauVec);
end
d2Gamma_dt2_star = Gamma_Integration_Function(del_s,d2Gamma_dt2_star_Integrand,root_idx) + d2Gamma_dt2_root_G_star;
end