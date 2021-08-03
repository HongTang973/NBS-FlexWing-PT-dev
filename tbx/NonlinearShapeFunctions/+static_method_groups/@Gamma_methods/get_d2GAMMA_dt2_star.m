function d2Gamma_dt2_star = get_d2GAMMA_dt2_star(Gamma_Integration_Function,del_s,d2E_dt2_star,root_idx,d2Gamma_dt2_root_G_star)
d2Gamma_dt2_star = Gamma_Integration_Function(del_s,d2E_dt2_star(:,2,:),root_idx) + d2Gamma_dt2_root_G_star;
end