function [dW_dq_Kinetic_Rotation_ddqComponent,dW_dq_Kinetic_Rotation_remainder] = get_dPi_dq_Kinetic_Rotation(...
    E_G, E_G_tr, I_varTheta,...
    dvarTheta_dq_G_Dim3xnq2ndxns, dvarTheta_dq_G_Dimnq2ndx3xns,...
    dvarTheta_dt_G,d2varTheta_dt2_G_star)

%[dvarTheta_dqa_G,dvarTheta_dqr456_G,dvarTheta_dt_G,d2varTheta_dt2_G_star]

%---------- common terms -----------

E_G_IvarTheta = utility_functions.MultiProd_(E_G,I_varTheta);
I_varTheta_G = utility_functions.MultiProd_(E_G_IvarTheta,E_G_tr);

dvarThdqtr_EG_IvarThI_EGtr = utility_functions.MultiProd_(dvarTheta_dq_G_Dimnq2ndx3xns,I_varTheta_G);

%-----------ddqComponent------------

dPi_dq_Kinetic_Rotation_ddqComponent = ...
    utility_functions.MultiProd_(dvarThdqtr_EG_IvarThI_EGtr,dvarTheta_dq_G_Dim3xnq2ndxns);

%-----------remaining_pt------------

dPi_dq_Kinetic_Rotation_remainder_pt1 = ...
    utility_functions.MultiProd_(dvarTheta_dq_G_Dimnq2ndx3xns,...
    utility_functions.MultiProd_(utility_functions.getSkewMat(dvarTheta_dt_G),...
    utility_functions.MultiProd_(I_varTheta_G,dvarTheta_dt_G)));

dPi_dq_Kinetic_Rotation_remainder_pt2 = ...
    utility_functions.MultiProd_(dvarThdqtr_EG_IvarThI_EGtr,d2varTheta_dt2_G_star);

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