function E_WE = get_E_WE(st,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp)
ey_WE = [ct_ss ; ct_cs ; st];
ex_WE = [cs_cp + st_ss_sp;- ss_cp + st_cs_sp;-ct_sp];
ez_WE = [cs_sp - st_ss_cp;- ss_sp - st_cs_cp; ct_cp];
E_WE = [ex_WE ey_WE ez_WE];
end