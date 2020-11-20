classdef dPi_dq_Kinetic_Translation < handle
    
    methods (Static = true)
        
        %method to return group 1 varTheta variables
        function [dW_dq_Kinetic_Translation_ddqComponent,dW_dq_Kinetic_Translation_remainder] = get_dPi_dq_Kinetic_Translation(...
                ms, dGamma_m_dq_G_Dim3xnq2ndxns, dGamma_m_dq_G_Dimnq2ndx3xns,...
                d2Gamma_m_dt2_G_star)
            
            %[dvarTheta_dqa_G,dvarTheta_dqr456_G,dvarTheta_dt_G,d2varTheta_dt2_G_star]
            
            %---------- common terms -----------
            
            dGamma_m_dq_tr_ms = bsxfun(@times, dGamma_m_dq_G_Dimnq2ndx3xns, ms);
            
            %-----------ddqComponent------------
            
            dPi_dq_Kinetic_Translation_ddqComponent = ...
                MultiProd_(dGamma_m_dq_tr_ms,dGamma_m_dq_G_Dim3xnq2ndxns);
            
            %-----------remaining_pt------------
            
            dPi_dq_Kinetic_Translation_remainder = ...
                MultiProd_(dGamma_m_dq_tr_ms,d2Gamma_m_dt2_G_star);
            
            %-----------------------------------
            
            
            %===================================
            dW_dq_Kinetic_Translation_ddqComponent = sum(dPi_dq_Kinetic_Translation_ddqComponent,3);
            dW_dq_Kinetic_Translation_remainder = sum(dPi_dq_Kinetic_Translation_remainder,3);
            %===================================
            
        end
        
    end
    
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %..
        %----------------------------------------------------------------------
        
        %----------------------------------------------------------------------
        
    end
    
end

