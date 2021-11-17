function out = custom_free_states(SimObject,Q,t) %#ok<INUSD>
    
%     %get the index of the rotational state (1 state in this case)
%     rR_idx  = SimObject.StateInfo{'Index','qRigidR'}{:};
%     %get the index of the rotational velocity state (1 state in this case)
%     drR_idx  = SimObject.StateInfo{'Index','dqRigidR'}{:};
%     
%     %get the index of the rotational state (1 state in this case)
%     rT_idx  = SimObject.StateInfo{'Index','qRigidT'}{:};
%     %get the index of the rotational velocity state (1 state in this case)
%     drT_idx  = SimObject.StateInfo{'Index','dqRigidT'}{:};
    
    %get the value of the rotational state
    qR_custom = Q(SimObject.rR_idx);
    %get the value of the rotational velocity state
    dqR_custom = Q(SimObject.drR_idx);

    %get the value of the rotational state
    qT_custom = Q(SimObject.rT_idx);
    %get the value of the rotational velocity state
    dqT_custom = Q(SimObject.drT_idx);
    %---------------------------------
    %in this function we wish to populate and return the following
    %kinematic quantities, relating them where appropriate to the custom rotational state
    
    %position in global system - constant
    out.rBarA_G = [0;0;qT_custom];
    %linear velocity in global system - zero
    out.drBarA_dt_G = [0;0;dqT_custom];
    
    %rotation in global system - custom state mapped to x component
    out.Beta_G = [0;qR_custom;0];
    %rotation matrix of above rotation
    out.R_G_A = utility_functions.r_matrix(out.Beta_G);
    %rotational velocity in global system - custom rate state mapped to x component
    out.Omega_G = [0;dqR_custom;0];
    
    %the derivative of the global to hub rotation matrix with respect to the rotational state
    skewY = utility_functions.getSkewMat([0;1;0]);
    out.dR_G_A_dqr_Dim3x3x1xnqr = utility_functions.MultiProd_(cat(4,zeros(3,3,1,1),skewY) , out.R_G_A);
    %the derivative of the global rotation vector with respect to the rotational state
    out.dvarTheta_dqr_G_Dim3x1x1xnqr = reshape([zeros(3,1) [0;1;0]],3,1,1,2);
    %the derivative of the global position vector with respect to the rotational state
    out.drBarA_G_dqr_G_Dim3x1x1xnqr = reshape([[0;0;1] zeros(3,1)],3,1,1,2);
    
end

