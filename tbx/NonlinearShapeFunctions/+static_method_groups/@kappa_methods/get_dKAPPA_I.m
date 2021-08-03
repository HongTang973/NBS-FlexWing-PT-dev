function [dKAPPA_dt_I,dKAPPA_dq_I_tr] = get_dKAPPA_I(st,ct,sp,cp,dth_ds,dsi_ds,nszrs,R_W_WE_tr,B_th,B_si,B_ph,dB_th,dB_si,dB_ph,dzeta_a_dt_tr,d2th_dsdt,d2si_dsdt,d2ph_dsdt,st_sp,st_cp,ct_sp,ct_cp,yidx2,TD)
dxi_dzeta_IE = [...                                            %e.g. dkappa_z/dth = dxi_dzeta(1,3)
    -dsi_ds.*st_sp              ,-dsi_ds.*ct , dsi_ds.*st_cp             ;
    nszrs                       , nszrs      , nszrs                     ;
    dsi_ds.*ct_cp - dth_ds.*sp , nszrs      , dsi_ds.*ct_sp + dth_ds.*cp];
%%%%%%%%%
dxi_dzetads_IE = [...
    cp     , nszrs   , sp     ;
    ct_sp ,-st      ,-ct_cp ;
    nszrs  , nszrs+1 , nszrs ];
%-----------------------------------
dxi_dzeta_I = utility_functions.MultiProd_(utility_functions.mult_Anmz_Bmp1(dxi_dzeta_IE,R_W_WE_tr),TD);
dxi_dzetads_I = utility_functions.MultiProd_(utility_functions.mult_Anmz_Bmp1(dxi_dzetads_IE,R_W_WE_tr),TD);
%---------------------------------------

dKAPPA_dq_I_tr = bsxfun(@times,[B_th;B_si;B_ph],dxi_dzeta_I(yidx2,:,:))+bsxfun(@times,[dB_th;dB_si;dB_ph],dxi_dzetads_I(yidx2,:,:));

dKAPPA_dt_I_tr = utility_functions.MultiProd_(dzeta_a_dt_tr,dxi_dzeta_I)+...
    utility_functions.MultiProd_([d2th_dsdt,d2si_dsdt,d2ph_dsdt],dxi_dzetads_I);
dKAPPA_dt_I = reshape(dKAPPA_dt_I_tr,3,1,[]);

end