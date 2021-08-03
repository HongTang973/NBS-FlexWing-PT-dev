function dGamma_dt = get_dGamma_dt(Gamma_Integration_Function,del_s,dE_dt_33ns,root_idx,dGamma_dt_root_G)
dGamma_dt = Gamma_Integration_Function(del_s,dE_dt_33ns(:,2,:),root_idx) + dGamma_dt_root_G;
end