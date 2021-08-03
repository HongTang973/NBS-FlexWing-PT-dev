function [dE_dt_W_91ns,dE_dqa_G_9qans] = get_dE(ct,st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,nszrs,ns,Rs_W_WE_tr,Rv_W_WE_tr,R_G_W,TD,Ba_tr,dzeta_a_dt_tr,yidx2)

dE_dth_WE_ = [...
    ct_ss_sp ,-st_ss ,-ct_ss_cp;
    ct_cs_sp ,-st_cs ,-ct_cs_cp;
    st_sp    , ct    ,-st_cp  ];
%%%%%%%%%
dE_dsi_WE_ = [...
    -ss_cp + st_cs_sp , ct_cs ,-ss_sp - st_cs_cp;
    -cs_cp - st_ss_sp ,-ct_ss ,-cs_sp + st_ss_cp;
    nszrs            , nszrs ,           nszrs];
%%%%%%%%%
dE_dph_WE_ = [...
    -cs_sp + st_ss_cp , nszrs , cs_cp + st_ss_sp;
    ss_sp + st_cs_cp , nszrs ,-ss_cp + st_cs_sp;
    -ct_cp            , nszrs ,-ct_sp          ];

dE_dzeta_WE = reshape([dE_dth_WE_,dE_dsi_WE_,dE_dph_WE_],9,3,ns);
dE_dzeta_W  = reshape(static_method_groups.ct_methods.map_WEtoW_stackDim2(dE_dzeta_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD),3,[],ns);
dE_dzeta_G  = static_method_groups.ct_methods.map_WtoG(dE_dzeta_W,R_G_W);

dE_dzeta_G_93 = reshape(dE_dzeta_G,9,3,ns);
dE_dqa_G_9qans = bsxfun(@times,Ba_tr,dE_dzeta_G_93(:,yidx2,:));

dE_dzeta_W_93 = reshape(dE_dzeta_W,9,3,ns);
dE_dt_W_91ns = sum(bsxfun(@times,dzeta_a_dt_tr,dE_dzeta_W_93),2);
end