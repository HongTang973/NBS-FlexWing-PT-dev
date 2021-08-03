function [KAPPA_I,dKAPPA_dt_I,dKAPPA_dqa_I_tr_Dim1x3xnsxnqa] = get_kappa_group1(...
    ...
    st,ct,sp,cp,dth_ds,dsi_ds,dph_ds,nszrs,R_W_WE,B_th,B_si,B_ph,dB_th,dB_si,dB_ph,st_sp,st_cp,ct_sp,ct_cp,yidx2,TD,...
    dzeta_a_dt_tr,d2th_dsdt,d2si_dsdt,d2ph_dsdt)

R_W_WE_tr = R_W_WE.';
TD_tr = permute(TD,[2 1 3]);

%-----------------------------------
kappa_x_IE = dsi_ds.*ct_sp + dth_ds.*cp;                       %(1)x(1)x(ns)
tau_IE     = dph_ds - dsi_ds.*st;                              %(1)x(1)x(ns)
kappa_z_IE =-dsi_ds.*ct_cp + dth_ds.*sp;                       %(1)x(1)x(ns)

KAPPA_IE = [kappa_x_IE ; tau_IE ; kappa_z_IE];
KAPPA_I = utility_functions.MultiProd_(TD_tr, utility_functions.mult_Anm1_Bmpz(R_W_WE,KAPPA_IE) );
%-----------------------------------

[dKAPPA_dt_I,dKAPPA_dqa_I_tr_Dimnqax3xns] = static_method_groups.kappa_methods.get_dKAPPA_I(st,ct,sp,cp,dth_ds,dsi_ds,nszrs,R_W_WE_tr,B_th,B_si,B_ph,dB_th,dB_si,dB_ph,dzeta_a_dt_tr,d2th_dsdt,d2si_dsdt,d2ph_dsdt,st_sp,st_cp,ct_sp,ct_cp,yidx2,TD);
dKAPPA_dqa_I_tr_Dim1x3xnsxnqa = permute(dKAPPA_dqa_I_tr_Dimnqax3xns,[4 2 3 1]);

end