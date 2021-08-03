classdef dPi_dq_Kinetic_Rotation < handle
    
    methods (Static = true)
        
        %method to return group 1 varTheta variables
        [dW_dq_Kinetic_Rotation_ddqComponent,dW_dq_Kinetic_Rotation_remainder] = get_dPi_dq_Kinetic_Rotation(...
            E_G, E_G_tr, I_varTheta,...
            dvarTheta_dq_G_Dim3xnq2ndxns, dvarTheta_dq_G_Dimnq2ndx3xns,...
            dvarTheta_dt_G,d2varTheta_dt2_G_star)

    end
    
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %..
        %----------------------------------------------------------------------
        
        %----------------------------------------------------------------------
        
    end
    
end



