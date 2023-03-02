classdef E_methods < handle
    %A collection of methods responsible for the calculation of Gamma and it's derivatives
    
    methods (Static = true)
        
        %method to return group 1 Gamma variables
        [E_WE,E_W,E_G,dE_dqa_G_3x3xnsxnqa,dE_dqe_G_3x3xnsxnqe,dE_dt_W] = get_E_group1(...
            ...
            st, ct,...
            st_ss,st_cs,st_sp,st_cp,...
            ct_ss,ct_cs,ct_sp,ct_cp,...
            ss_sp,ss_cp,cs_sp,cs_cp,...
            st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,...
            ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,...
            ...
            Rs_W_WE_tr,Rv_W_WE_tr,TD,nszrs,ns,...
            R_A_W,R_G_W,...
            ...
            dzeta_a_dt_tr,Ba_tr,yidx2,...
            ...
            dR_G_W_dqe_3x3x1xnqe,FLAG_free_free)
        
        %st,ss,sp,ct,cs,cp,st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp
        
        
        %d2E_dzeta2_W = E_methods.get_d2E_dzeta2_W(st,nszrs,Rs_W_WE_tr,Rv_W_WE_tr,st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,TD);
        
        
        %         +=
        %
        %
        %
        %         d2varTheta_dt2_W_star = ...
        %                 bsxfun(@times,dth_dt.*dsi_dt,dthdsiVec_W) +...
        %                 bsxfun(@times,dph_dt,dey_dt_W);
        %
        %             d2varTheta_dt2_G_star = ...
        %                 mult_Anm1_Bmpz(OmegaSkew_G,dvarTheta_dt_G)...
        %                +ct_methods.map_WtoG(d2varTheta_dt2_W_star,R_G_W);
        %
        %             if ~FLAG_free_free
        %                 d2varTheta_dt2_G_star = bsxfun(@plus,dOmega_dt_G,d2varTheta_dt2_G_star);
        %             end
        %
        %
        
        
    end
    
    
    
    
    methods (Static = true)%, Access = private) %Group 1 Static Methods
        
        %E_WE
        %----------------------------------------------------------------------
        E_WE = get_E_WE(st,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp)
        
        %----------------------------------------------------------------------
        
        %dE
        %----------------------------------------------------------------------
        [dE_dt_W_91ns,dE_dqa_G_9qans] = get_dE(ct,st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,nszrs,ns,Rs_W_WE_tr,Rv_W_WE_tr,R_G_W,TD,Ba_tr,dzeta_a_dt_tr,yidx2)
        
        
        %----------------------------------------------------------------------
        
        %     %d2E_dzeta2
        %     %----------------------------------------------------------------------
        %     function d2E_dzeta2_W = get_d2E_dzeta2_W(st,nszrs,Rs_W_WE_tr,Rv_W_WE_tr,st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,TD)
        %         DDex1_WE = [
        %             -st_ss_sp, ct_cs_sp       , ct_ss_cp       ,...
        %              ct_cs_sp,-cs_cp-st_ss_sp, ss_sp+st_cs_cp,...
        %              ct_ss_cp, ss_sp+st_cs_cp,-cs_cp-st_ss_sp];                   %(1)x(9)x(ns)
        %
        %         DDex2_WE = [
        %             -st_cs_sp,-ct_ss_sp       , ct_cs_cp       ,...
        %             -ct_ss_sp, ss_cp-st_cs_sp, cs_sp-st_ss_cp,...
        %              ct_cs_cp, cs_sp-st_ss_cp, ss_cp-st_cs_sp];                   %(1)x(9)x(ns)
        %
        %         DDex3_WE = [
        %              ct_sp, nszrs, st_cp,...
        %               nszrs, nszrs,  nszrs,...
        %              st_cp, nszrs, ct_sp];                                        %(1)x(9)x(ns)
        %
        %         DDey1_WE = [
        %             -ct_ss,-st_cs, nszrs,...
        %             -st_cs,-ct_ss, nszrs,...
        %               nszrs,  nszrs, nszrs];                                      %(1)x(9)x(ns)
        %
        %         DDey2_WE = [
        %             -ct_cs, st_ss, nszrs,...
        %              st_ss,-ct_cs, nszrs,...
        %               nszrs,  nszrs, nszrs];                                      %(1)x(9)x(ns)
        %
        %         DDey3_WE = [
        %                -st, nszrs, nszrs,...
        %              nszrs, nszrs, nszrs,...
        %              nszrs, nszrs, nszrs];                                        %(1)x(9)x(ns)
        %
        %         DDez1_WE = [
        %              st_ss_cp,-ct_cs_cp       , ct_ss_sp       ,...
        %             -ct_cs_cp,-cs_sp+st_ss_cp,-ss_cp+st_cs_sp,...
        %              ct_ss_sp,-ss_cp+st_cs_sp,-cs_sp+st_ss_cp];                   %(1)x(9)x(ns)
        %
        %         DDez2_WE = [
        %              st_cs_cp, ct_ss_cp       , ct_cs_sp       ,...
        %              ct_ss_cp, ss_sp+st_cs_cp,-cs_cp-st_ss_sp,...
        %              ct_cs_sp,-cs_cp-st_ss_sp, ss_sp+st_cs_cp];                   %(1)x(9)x(ns)
        %
        %         DDez3_WE = [
        %             -ct_cp, nszrs, st_sp,...
        %               nszrs, nszrs,  nszrs,...
        %              st_sp, nszrs,-ct_cp];                                        %(1)x(9)x(ns)
        %
        %         d2E_dzeta2_WE = [DDex1_WE;DDex2_WE;DDex3_WE;DDey1_WE;DDey2_WE;DDey3_WE;DDez1_WE;DDez2_WE;DDez3_WE]; %(9)x(9)x(ns)
        %
        %         d2E_dzeta2_W = ct_methods.map_WEtoW_stackDim2(d2E_dzeta2_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD); %(9)x(9)x(ns)
        %
        %     end
        
        
    end
    
end





% dGamma_dq_G
% dGamma_dt_W
%
% Gamma_m_W = Gamma_W + MultiProd(E_W,massOffset_I);
%     dGamma_m_dq_G = dGamma_dq_G + bsxfun(@times,massOffset_am,dex_dq_G)+ bsxfun(@times,massOffset_bm,dey_dq_G) + bsxfun(@times,massOffset_cm,dez_dq_G);
%     dGamma_m_dt_W
