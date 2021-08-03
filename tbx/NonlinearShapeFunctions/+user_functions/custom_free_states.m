function out = custom_free_states(SimObject,Q,t) %#ok<INUSD>
    
    %get the index of the rotational state (1 state in this case)
    rR_idx  = SimObject.StateInfo{'Index','qRigidR'}{:};
    %get the index of the rotational velocity state (1 state in this case)
    drR_idx  = SimObject.StateInfo{'Index','qRigidR'}{:};
    
    %get the value of the rotational state
    q_custom = Q(rR_idx);
    %get the value of the rotational velocity state
    dq_custom = Q(drR_idx);

    %---------------------------------
    %in this function we wish to populate and return the following
    %kinematic quantities, relating them where appropriate to the custom rotational state
    
    %position in global system - constant
    out.rBarA_G = [0;0;0];
    %linear velocity in global system - zero
    out.drBarA_dt_G = [0;0;0];
    
    %rotation in global system - custom state mapped to x component
    out.Beta_G = [q_custom;0;0];
    %rotation matrix of above rotation
    out.R_G_A = r_matrix(out.Beta_G);
    %rotational velocity in global system - custom rate state mapped to x component
    out.Omega_G = [dq_custom;0;0];
    
    %the derivative of the global to hub rotation matrix with respect to the rotational state
    out.dR_G_A_dqr_Dim3x3x1xnqr = getSkewMat([1;0;0]) * out.R_G_A;
    %the derivative of the global rotation vector with respect to the rotational state
    out.dvarTheta_dqr_G_Dim3x1x1xnqr = reshape([1;0;0],3,1,1,[]);
    %the derivative of the global position vector with respect to the rotational state (zero)
    out.drBarA_G_dqr_G_Dim3x1x1xnqr = reshape([0;0;0],3,1,1,[]);

end

