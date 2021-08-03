classdef varTheta_methods < handle
    %A collection of methods responsible for the calculation of varTheta and it's derivatives
    
    methods (Static = true)
        
        %method to return group 1 varTheta variables
        [dvarTheta_dqa_G_3x1xnsxnqa,dvarTheta_dt_G,d2varTheta_dt2_G_star] = get_varTheta_group1(...
                ...
                ss, cs, ey_WE, dE_dt_W,...
                ...
                nszrs,ns,...
                R_W_WE,R_G_W,...
                dvarTheta_dt_root_G,dvarTheta_dt_root_G_skew,...
                ...
                dth_dt,dsi_dt,dph_dt,...
                B_th_tr,B_si_tr,B_ph_tr,...
                ...
                d2varTheta_dt2_root_G_star)
            
  
    end
    
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %dvarTheta
        %----------------------------------------------------------------------
        [dvarTheta_dq_G,dvarTheta_dt_G] = get_dvarTheta(dthVec_G,dsiVec_G,dphVec_G,B_th_tr,B_si_tr,B_ph_tr,dzeta_a_dt_tr,dvarTheta_dt_root_G)
           
            %---------------------------

        %----------------------------------------------------------------------
        
    end
    
end

