classdef dPi_dq_Kinetic_Rotation < handle
    
    methods (Static = true)
        
        %method to return group 1 varTheta variables
        function [dW_dq_Kinetic_Rotation_ddqComponent,dW_dq_Kinetic_Rotation_remainder] = get_dPi_dq_Kinetic_Rotation(...
                E_G, E_G_tr, I_varTheta,...
                dvarTheta_dq_G_Dim3xnq2ndxns, dvarTheta_dq_G_Dimnq2ndx3xns,...
                dvarTheta_dt_G,d2varTheta_dt2_G_star)
            
            %[dvarTheta_dqa_G,dvarTheta_dqr456_G,dvarTheta_dt_G,d2varTheta_dt2_G_star]
            
            %---------- common terms -----------
            
            E_G_IvarTheta = MultiProd_(E_G,I_varTheta);
            I_varTheta_G = MultiProd_(E_G_IvarTheta,E_G_tr);
            
            dvarThdqtr_EG_IvarThI_EGtr = MultiProd_(dvarTheta_dq_G_Dimnq2ndx3xns,I_varTheta_G);
            
            %-----------ddqComponent------------
            
            dPi_dq_Kinetic_Rotation_ddqComponent = ...
                MultiProd_(dvarThdqtr_EG_IvarThI_EGtr,dvarTheta_dq_G_Dim3xnq2ndxns);
            
            %-----------remaining_pt------------
            
            dPi_dq_Kinetic_Rotation_remainder_pt1 = ...
                MultiProd_(dvarTheta_dq_G_Dimnq2ndx3xns,MultiProd_(getSkewMat(dvarTheta_dt_G),MultiProd_(I_varTheta_G,dvarTheta_dt_G)));
            
            dPi_dq_Kinetic_Rotation_remainder_pt2 = ...
                MultiProd_(dvarThdqtr_EG_IvarThI_EGtr,d2varTheta_dt2_G_star);
            
            %------------
            dPi_dq_Kinetic_Rotation_remainder = ...
                dPi_dq_Kinetic_Rotation_remainder_pt1 ...
               +dPi_dq_Kinetic_Rotation_remainder_pt2;
            %-----------------------------------
            
            
            %===================================
            dW_dq_Kinetic_Rotation_ddqComponent = sum(dPi_dq_Kinetic_Rotation_ddqComponent,3);
            dW_dq_Kinetic_Rotation_remainder = sum(dPi_dq_Kinetic_Rotation_remainder,3);
            %===================================
            
        end
        
    end
    
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %..
        %----------------------------------------------------------------------
        
        %----------------------------------------------------------------------
        
    end
    
end

