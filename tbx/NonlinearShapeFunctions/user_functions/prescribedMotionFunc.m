function prescribed_motion = prescribedMotionFunc(SimObject,t)
    
    omega_ = 1; %angular velocity
    Omega_G = [0;0;1].*omega_; %Angular velocity vector
    beta_ = t*Omega_G; %Angle
    Beta_G = [0;0;1].*beta_; %Angle vector
    dOmega_dt_G = [0;0;0]; %Angular acceleration vector
    
    OmegaSkew_G = getSkewMat(Omega_G); dOmega_dtSkew_G = getSkewMat(dOmega_dt_G);
    
    R_G_A = r_matrix(Beta_G);
    dR_G_A_dt = OmegaSkew_G*R_G_A;
    d2R_G_A_dt2 = dOmega_dtSkew_G*R_G_A + OmegaSkew_G*OmegaSkew_G*R_G_A;
    
    prescribed_motion.Omega_G = Omega_G;
    prescribed_motion.OmegaSkew_G = getSkewMat(Omega_G);
    prescribed_motion.dOmega_dt_G = dOmega_dt_G;
    prescribed_motion.R_G_A = R_G_A;
    prescribed_motion.dR_G_A_dt = dR_G_A_dt;
    prescribed_motion.d2R_G_A_dt2 = d2R_G_A_dt2;
    prescribed_motion.rBarA_G = [0;0;0]; %Translation
    prescribed_motion.drBarA_dt_G = [0;0;0]; %Translational velocity
    prescribed_motion.d2rBarA_dt2_G = [0;0;0]; %Translational acceleration
end