classdef Gamma_methods < handle
    %A collection of methods responsible for the calculation of Gamma and it's derivatives
    
    
methods (Static = true)
    
    %method to return group 1 Gamma variables
    function [Gamma_G,dGamma_dt_G,dGamma_dq_G_Dim3x1xnsxnq2nd,dGamma_m_dq_G_Dim3x1xnsxnq2nd,d2Gamma_dt2_G_star,d2Gamma_m_dt2_G_star] = get_Gamma_group1(...
            Gamma_Integration_Function,del_s,ns,root_idx,...
            E_W,dE_dt_G,d2E_dt2_G_star,E_G,dE_dq_G_Dim3x3xnsxnq2nd,...
            tau_x,tau_y,tau_z,...
            dtau_x_dt,dtau_y_dt,dtau_z_dt,...
            dtau_x_dqs,dtau_y_dqs,dtau_z_dqs,...
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
            Gamma_G = Gamma_methods.get_Gamma(Gamma_Integration_Function,del_s,E_G,root_idx,Gamma_root_G);
            dGamma_dt_G = Gamma_methods.get_dGamma_dt(Gamma_Integration_Function,del_s,dE_dt_G,root_idx,dGamma_dt_root_G);
            d2Gamma_dt2_G_star = Gamma_methods.get_d2GAMMA_dt2_star(Gamma_Integration_Function,del_s,d2E_dt2_G_star,root_idx,d2Gamma_dt2_root_G_star);
            dGamma_dq_G_Dim3x1xnsxnq2nd = Gamma_methods.get_dGamma_dq(Gamma_Integration_Function,del_s,dey_dq_G_Dim3x1xnsxnq2nd,root_idx,dGamma_dq_root_G);
            
        else
            TauVec = [tau_x ; 1 + tau_y ; tau_z];
            dTauVec_dt = [dtau_x_dt ; dtau_y_dt ; dtau_z_dt];
            
            Gamma_G = Gamma_methods.get_Gamma_shear(Gamma_Integration_Function,del_s,E_W,TauVec,root_idx,Gamma_root_G);
            dGamma_dt_G = Gamma_methods.get_dGamma_dt_shear(Gamma_Integration_Function,del_s,E_G,dE_dt_G,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,dGamma_dt_root_G);
            d2Gamma_dt2_G_star = Gamma_methods.get_d2GAMMA_dt2_star_shear(Gamma_Integration_Function,del_s,dE_dt_G,d2E_dt2_G_star,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,d2Gamma_dt2_root_G_star);
            
            dGamma_dq_G_Dim3x1xnsxnq2nd = Gamma_methods.get_dGamma_dq_shear(Gamma_Integration_Function,del_s,...
                ex_G,ey_G,ez_G,...
                dex_dq_G_Dim3x1xnsxnq2nd,dey_dq_G_Dim3x1xnsxnq2nd,dez_dq_G_Dim3x1xnsxnq2nd,...
                tau_x,tau_y,tau_z,...
                dtau_x_dqs,dtau_y_dqs,dtau_z_dqs,...
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
            dGamma_m_dq_G_Dim3x1xnsxnq2nd = Gamma_methods.get_dGamma_dq_offset(dGamma_dq_G_Dim3x1xnsxnq2nd,dE_dq_G_Dim3x3xnsxnq2nd,massOffset_I);
            d2Gamma_m_dt2_G_star = Gamma_methods.get_d2Gamma_dt2_star_offset(d2Gamma_dt2_G_star,d2E_dt2_G_star,massOffset_I);
        end
        
    end
    
end


methods (Static = true, Access = private) %Group 1 Static Methods
    
    %GAMMA
    %----------------------------------------------------------------------
    function Gamma = get_Gamma(Gamma_Integration_Function,del_s,E,root_idx,Gamma_root)
        Gamma = Gamma_Integration_Function(del_s,E(:,2,:),root_idx) + Gamma_root;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~
    function Gamma = get_Gamma_shear(Gamma_Integration_Function,del_s,E,TauVec,root_idx,Gamma_root)
        Gamma_Integrand = MultiProd_(E,TauVec);
        Gamma = Gamma_Integration_Function(del_s,Gamma_Integrand,root_idx) + Gamma_root;
    end
    %----------------------------------------------------------------------
    
    %dGAMMA_dt
    %----------------------------------------------------------------------
    function dGamma_dt = get_dGamma_dt(Gamma_Integration_Function,del_s,dE_dt_33ns,root_idx,dGamma_dt_root_G)
        dGamma_dt = Gamma_Integration_Function(del_s,dE_dt_33ns(:,2,:),root_idx) + dGamma_dt_root_G;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~
    function dGamma_dt = get_dGamma_dt_shear(Gamma_Integration_Function,del_s,E,dE_dt_33ns,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,dGamma_dt_root_G)
        if Gamma_approx_lvl == 0
            dGamma_dt_Integrand = MultiProd_(dE_dt_33ns,TauVec)+MultiProd_(E,dTauVec_dt);
        elseif Gamma_approx_lvl == 1
            dGamma_dt_Integrand = MultiProd_(dE_dt_33ns,TauVec);
        end
        dGamma_dt = Gamma_Integration_Function(del_s,dGamma_dt_Integrand,root_idx) + dGamma_dt_root_G;
    end
    %----------------------------------------------------------------------
    
    %d2GAMMA_dt2_star
    %----------------------------------------------------------------------
    function d2Gamma_dt2_star = get_d2GAMMA_dt2_star(Gamma_Integration_Function,del_s,d2E_dt2_star,root_idx,d2Gamma_dt2_root_G_star)
        d2Gamma_dt2_star = Gamma_Integration_Function(del_s,d2E_dt2_star(:,2,:),root_idx) + d2Gamma_dt2_root_G_star;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~
    function d2Gamma_dt2_star = get_d2GAMMA_dt2_star_shear(Gamma_Integration_Function,del_s,dE_dt,d2E_dt2_star,TauVec,dTauVec_dt,root_idx,Gamma_approx_lvl,d2Gamma_dt2_root_G_star)
        if Gamma_approx_lvl == 0
            d2Gamma_dt2_star_Integrand = MultiProd_(d2E_dt2_star,TauVec)+2*MultiProd_(dE_dt,dTauVec_dt);
        elseif Gamma_approx_lvl == 1
            d2Gamma_dt2_star_Integrand = MultiProd_(d2E_dt2_star,TauVec);
        end
        d2Gamma_dt2_star = Gamma_Integration_Function(del_s,d2Gamma_dt2_star_Integrand,root_idx) + d2Gamma_dt2_root_G_star;
    end
    %----------------------------------------------------------------------
    
    %dGAMMA_dq
    %----------------------------------------------------------------------
    function dGamma_dq = get_dGamma_dq(Gamma_Integration_Function,del_s,dey_dq,root_idx,dGamma_dq_root_G)
        dGamma_dq = bsxfun(@plus, Gamma_Integration_Function(del_s,dey_dq,root_idx) , dGamma_dq_root_G);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~
    function dGamma_dq = get_dGamma_dq_shear(Gamma_Integration_Function,del_s,ex,ey,ez,dex_dq,dey_dq,dez_dq,tau_x,tau_y,tau_z,dtau_x_dq,dtau_y_dq,dtau_z_dq,root_idx,dGamma_dq_root_G)
        
%         dGamma_dq_Integrand = ...
%             bsxfun(@times,dex_dq,tau_x)+...
%             bsxfun(@times,dey_dq,tau_y)+...
%             bsxfun(@times,dez_dq,tau_z)+ MultiProd_(E,dTauVec_dq);
        
        if ~isempty(dex_dq) && ~isempty(dey_dq) && ~isempty(dez_dq)
            dGamma_dq_Integrand_pt1 = ...
                bsxfun(@times,dex_dq,tau_x)+...
                bsxfun(@times,dey_dq,1+tau_y)+...
                bsxfun(@times,dez_dq,tau_z);
        else
            dGamma_dq_Integrand_pt1 = 0;
        end
    
        if ~isempty(dtau_x_dq) && ~isempty(dtau_y_dq) && ~isempty(dtau_z_dq)
            dGamma_dq_Integrand_pt2 = [...
                bsxfun(@times,ex,dtau_x_dq),...
                bsxfun(@times,ey,dtau_y_dq),...
                bsxfun(@times,ez,dtau_z_dq)];
        else
            dGamma_dq_Integrand_pt2 = 0;
        end
        
        dGamma_dq_Integrand = dGamma_dq_Integrand_pt1 + dGamma_dq_Integrand_pt2;
        
        dGamma_dq = Gamma_Integration_Function(del_s,dGamma_dq_Integrand,root_idx) + dGamma_dq_root_G;
    end
    %----------------------------------------------------------------------
    
    %Gamma_offset
    %----------------------------------------------------------------------
    function Gamma_offset = get_Gamma_W_offset(Gamma,E,massOffset_I)
        Gamma_offset = Gamma + MultiProd_(E,massOffset_I);
    end
    %----------------------------------------------------------------------
    
    %dGamma_dt_offset
    %----------------------------------------------------------------------
    function dGamma_dt_offset = get_dGamma_dt_offset(dGamma_dt,dE_dt,massOffset_I)
        dGamma_dt_offset = dGamma_dt + MultiProd_(dE_dt,massOffset_I);
    end
    %----------------------------------------------------------------------
    
    %dGamma_dq_offset
    %----------------------------------------------------------------------
    function dGamma_dq_offset = get_dGamma_dq_offset(dGamma_dq,dE_dq,massOffset_I)
        dGamma_dq_offset = dGamma_dq + MultiProd_(dE_dq,massOffset_I);
    end
    %----------------------------------------------------------------------
    
    
    function d2Gamma_dt2_star_offset = get_d2Gamma_dt2_star_offset(d2Gamma_dt2_star,d2E_dt2_star,massOffset_I)
        d2Gamma_dt2_star_offset = d2Gamma_dt2_star + MultiProd_(d2E_dt2_star,massOffset_I);
    end
    
    
end

end





% dGamma_dq_G
% dGamma_dt_W
% 
% Gamma_m_W = Gamma_W + MultiProd(E_W,massOffset_I);
%     dGamma_m_dq_G = dGamma_dq_G + bsxfun(@times,massOffset_am,dex_dq_G)+ bsxfun(@times,massOffset_bm,dey_dq_G) + bsxfun(@times,massOffset_cm,dez_dq_G);
%     dGamma_m_dt_W
