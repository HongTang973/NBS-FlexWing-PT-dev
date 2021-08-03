classdef Gamma_methods < handle
    %A collection of methods responsible for the calculation of Gamma and it's derivatives
    
    
    methods (Static = true)
        
        %method to return group 1 Gamma variables
        [Gamma_G,dGamma_dt_G,dGamma_dq_G_Dim3x1xnsxnq2nd,dGamma_m_dq_G_Dim3x1xnsxnq2nd,d2Gamma_dt2_G_star,d2Gamma_m_dt2_G_star] = get_Gamma_group1(...
            Gamma_Integration_Function,del_s,ns,root_idx,...
            E_W,dE_dt_G,d2E_dt2_G_star,E_G,dE_dq_G_Dim3x3xnsxnq2nd,...
            tau_x,tau_y,tau_z,...
            dtau_x_dt,dtau_y_dt,dtau_z_dt,...
            dtau_x_dq,dtau_y_dq,dtau_z_dq,...
            massOffset_I,FLAG_massOffset,...                               %dR_G_W_dqr456_931,FLAG_free_free,...
            FLAG_shear,Gamma_approx_lvl,...
            Gamma_root_G,dGamma_dt_root_G,dGamma_dq_root_G,d2Gamma_dt2_root_G_star)
        
    end
    
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %GAMMA
        %----------------------------------------------------------------------
        Gamma = get_Gamma(Gamma_Integration_Function,del_s,E,root_idx,Gamma_root)
        
        %~~~~~~~~~~~~~~~~~~~~~~~
        Gamma = get_Gamma_shear(Gamma_Integration_Function,del_s,E,TauVec,root_idx,Gamma_root)
        
        %----------------------------------------------------------------------
        
        %dGAMMA_dt
        %----------------------------------------------------------------------
        dGamma_dt = get_dGamma_dt(Gamma_Integration_Function,del_s,dE_dt_33ns,root_idx,dGamma_dt_root_G)
        
        %~~~~~~~~~~~~~~~~~~~~~~~
        dGamma_dt = get_dGamma_dt_shear(Gamma_Integration_Function,del_s,E,dE_dt_33ns,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,dGamma_dt_root_G)
        
        %----------------------------------------------------------------------
        
        %d2GAMMA_dt2_star
        %----------------------------------------------------------------------
        d2Gamma_dt2_star = get_d2GAMMA_dt2_star(Gamma_Integration_Function,del_s,d2E_dt2_star,root_idx,d2Gamma_dt2_root_G_star)
        
        %~~~~~~~~~~~~~~~~~~~~~~~
        d2Gamma_dt2_star = get_d2GAMMA_dt2_star_shear(Gamma_Integration_Function,del_s,dE_dt,d2E_dt2_star,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,d2Gamma_dt2_root_G_star)
        
        %----------------------------------------------------------------------
        
        %dGAMMA_dq
        %----------------------------------------------------------------------
        dGamma_dq = get_dGamma_dq(Gamma_Integration_Function,del_s,dey_dq,root_idx,dGamma_dq_root_G)
        
        %~~~~~~~~~~~~~~~~~~~~~~~
        dGamma_dq = get_dGamma_dq_shear(Gamma_Integration_Function,del_s,ex,ey,ez,dex_dq,dey_dq,dez_dq,tau_x,tau_y,tau_z,dtau_x_dq,dtau_y_dq,dtau_z_dq,root_idx,dGamma_dq_root_G)
        
        
        %----------------------------------------------------------------------
        
        %Gamma_offset
        %----------------------------------------------------------------------
        Gamma_offset = get_Gamma_W_offset(Gamma,E,massOffset_I)
        
        %----------------------------------------------------------------------
        
        %dGamma_dt_offset
        %----------------------------------------------------------------------
        dGamma_dt_offset = get_dGamma_dt_offset(dGamma_dt,dE_dt,massOffset_I)
        
        %----------------------------------------------------------------------
        
        %dGamma_dq_offset
        %----------------------------------------------------------------------
        dGamma_dq_offset = get_dGamma_dq_offset(dGamma_dq,dE_dq,massOffset_I)
        
        %----------------------------------------------------------------------
        
        
        d2Gamma_dt2_star_offset = get_d2Gamma_dt2_star_offset(d2Gamma_dt2_star,d2E_dt2_star,massOffset_I)
        
        
    end
    
end

% dGamma_dq_G
% dGamma_dt_W
%
% Gamma_m_W = Gamma_W + MultiProd(E_W,massOffset_I);
%     dGamma_m_dq_G = dGamma_dq_G + bsxfun(@times,massOffset_am,dex_dq_G)+ bsxfun(@times,massOffset_bm,dey_dq_G) + bsxfun(@times,massOffset_cm,dez_dq_G);
%     dGamma_m_dt_W
