function [dvarTheta_dq_G,dvarTheta_dt_G] = get_dvarTheta(dthVec_G,dsiVec_G,dphVec_G,B_th_tr,B_si_tr,B_ph_tr,dzeta_a_dt_tr,dvarTheta_dt_root_G)
dvarTheta_dq_G = [B_th_tr.*dthVec_G, B_si_tr.*dsiVec_G, B_ph_tr.*dphVec_G];
%---------------------------
dvarTheta_dt_G = ...
    sum([dthVec_G,dsiVec_G,dphVec_G].*[dzeta_a_dt_tr;dzeta_a_dt_tr;dzeta_a_dt_tr],2) + dvarTheta_dt_root_G;
end
