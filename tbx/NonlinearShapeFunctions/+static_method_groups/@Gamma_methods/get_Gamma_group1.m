function [Gamma_G,dGamma_dt_G,dGamma_dq_G_Dim3x1xnsxnq2nd,dGamma_m_dq_G_Dim3x1xnsxnq2nd,d2Gamma_dt2_G_star,d2Gamma_m_dt2_G_star] = get_Gamma_group1(...
    Gamma_Integration_Function,del_s,ns,root_idx,...
    E_W,dE_dt_G,d2E_dt2_G_star,E_G,dE_dq_G_Dim3x3xnsxnq2nd,...
    tau_x,tau_y,tau_z,...
    dtau_x_dt,dtau_y_dt,dtau_z_dt,...
    dtau_x_dq,dtau_y_dq,dtau_z_dq,...
    massOffset_I,FLAG_massOffset,...                               %dR_G_W_dqr456_931,FLAG_free_free,...
    FLAG_shear,Gamma_approx_lvl,...
    Gamma_root_G,dGamma_dt_root_G,dGamma_dq_root_G,d2Gamma_dt2_root_G_star)

ex_G = E_G(:,1,:);
ey_G = E_G(:,2,:);
ez_G = E_G(:,3,:);

dex_dq_G_Dim3x1xnsxnq2nd = dE_dq_G_Dim3x3xnsxnq2nd(:,1,:,:);
dey_dq_G_Dim3x1xnsxnq2nd = dE_dq_G_Dim3x3xnsxnq2nd(:,2,:,:);
dez_dq_G_Dim3x1xnsxnq2nd = dE_dq_G_Dim3x3xnsxnq2nd(:,3,:,:);

if FLAG_shear == false
    Gamma_G = static_method_groups.Gamma_methods.get_Gamma(Gamma_Integration_Function,del_s,E_G,root_idx,Gamma_root_G);
    dGamma_dt_G = static_method_groups.Gamma_methods.get_dGamma_dt(Gamma_Integration_Function,del_s,dE_dt_G,root_idx,dGamma_dt_root_G);
    d2Gamma_dt2_G_star = static_method_groups.Gamma_methods.get_d2GAMMA_dt2_star(Gamma_Integration_Function,del_s,d2E_dt2_G_star,root_idx,d2Gamma_dt2_root_G_star);
    dGamma_dq_G_Dim3x1xnsxnq2nd = static_method_groups.Gamma_methods.get_dGamma_dq(Gamma_Integration_Function,del_s,dey_dq_G_Dim3x1xnsxnq2nd,root_idx,dGamma_dq_root_G);
    
else
    TauVec = [tau_x ; 1 + tau_y ; tau_z];
    dTauVec_dt = [dtau_x_dt ; dtau_y_dt ; dtau_z_dt];
    
    Gamma_G = static_method_groups.Gamma_methods.get_Gamma_shear(Gamma_Integration_Function,del_s,E_W,TauVec,root_idx,Gamma_root_G);
    dGamma_dt_G = static_method_groups.Gamma_methods.get_dGamma_dt_shear(Gamma_Integration_Function,del_s,E_G,dE_dt_G,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,dGamma_dt_root_G);
    d2Gamma_dt2_G_star = static_method_groups.Gamma_methods.get_d2GAMMA_dt2_star_shear(Gamma_Integration_Function,del_s,dE_dt_G,d2E_dt2_G_star,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,d2Gamma_dt2_root_G_star);
    
    dGamma_dq_G_Dim3x1xnsxnq2nd = static_method_groups.Gamma_methods.get_dGamma_dq_shear(Gamma_Integration_Function,del_s,...
        ex_G,ey_G,ez_G,...
        dex_dq_G_Dim3x1xnsxnq2nd,dey_dq_G_Dim3x1xnsxnq2nd,dez_dq_G_Dim3x1xnsxnq2nd,...
        tau_x,tau_y,tau_z,...
        dtau_x_dq,dtau_y_dq,dtau_z_dq,...
        root_idx,dGamma_dq_root_G);
    
end

if ~FLAG_massOffset
    dGamma_m_dq_G_Dim3x1xnsxnq2nd = dGamma_dq_G_Dim3x1xnsxnq2nd;
    d2Gamma_m_dt2_G_star = d2Gamma_dt2_G_star;
    %             dGamma_m_dt_W = dGamma_dt_W;
    %             dGamma_m_dt_G = dGamma_dt_G;
else
    %            Gamma_m_W = get_Gamma_W_offset(Gamma_W,E_W,massOffset_I);
    %            dGamma_m_dt_W = get_dGamma_dt_offset(dGamma_dt_W,dE_dt_W,massOffset_I);
    %            dGamma_m_dt_G = get_dGamma_dt_offset(dGamma_dt_G,dE_dt_G,massOffset_I);
    dGamma_m_dq_G_Dim3x1xnsxnq2nd = static_method_groups.Gamma_methods.get_dGamma_dq_offset(dGamma_dq_G_Dim3x1xnsxnq2nd,dE_dq_G_Dim3x3xnsxnq2nd,massOffset_I);
    d2Gamma_m_dt2_G_star = static_method_groups.Gamma_methods.get_d2Gamma_dt2_star_offset(d2Gamma_dt2_G_star,d2E_dt2_G_star,massOffset_I);
end

end