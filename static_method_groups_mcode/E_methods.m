classdef E_methods < handle
    %A collection of methods responsible for the calculation of Gamma and it's derivatives
    
    
methods (Static = true)
    
    %method to return group 1 Gamma variables
    function [E_WE,E_W,E_G,dE_dqa_G_3x3xnsxnqa,dE_dqe_G_3x3xnsxnqe,dE_dt_W] = get_E_group1(...
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
        
        E_WE = E_methods.get_E_WE(st,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp);
        
        E_W = ct_methods.map_WEtoW(E_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD);
        
        E_G = ct_methods.map_AtoG(E_W,R_G_W);
        
        [dE_dt_W_91ns,dE_dqa_G_9qans] = E_methods.get_dE(ct,st_ss,st_cs,st_sp,st_cp,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp,ct_ss_sp,ct_ss_cp,ct_cs_sp,ct_cs_cp,nszrs,ns,Rs_W_WE_tr,Rv_W_WE_tr,R_G_W,TD,Ba_tr,dzeta_a_dt_tr,yidx2);
        dE_dt_W = reshape(dE_dt_W_91ns,3,3,ns);
        dE_dqa_G_3x3xnsxnqa = reshape(permute(dE_dqa_G_9qans,[1 4 3 2]),3,3,ns,[]);
        
        if ~isempty(dR_G_W_dqe_3x3x1xnqe)
            dE_dqe_G_3x3xnsxnqe = MultiProd_(dR_G_W_dqe_3x3x1xnqe,E_W);
            %dE_dqe_G = reshape([dE_dqe(1:3,:,:),dE_dqe(4:6,:,:),dE_dqe(7:9,:,:)],9,3,ns);
%             dE_dqe_G = [...
%                 reshape(dE_dqe_G_3nq_3(:,1,:),3,[],ns) ;...
%                 reshape(dE_dqe_G_3nq_3(:,2,:),3,[],ns) ;...
%                 reshape(dE_dqe_G_3nq_3(:,3,:),3,[],ns)];
        else
            dE_dqe_G_3x3xnsxnqe = zeros(3,3,ns,0);
        end
        
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
    
end


methods (Static = true, Access = private) %Group 1 Static Methods
    
    %E_WE
    %----------------------------------------------------------------------
    function E_WE = get_E_WE(st,ct_ss,ct_cs,ct_sp,ct_cp,ss_sp,ss_cp,cs_sp,cs_cp,st_ss_sp,st_ss_cp,st_cs_sp,st_cs_cp)
        ey_WE = [ct_ss ; ct_cs ; st];
        ex_WE = [cs_cp + st_ss_sp;- ss_cp + st_cs_sp;-ct_sp];
        ez_WE = [cs_sp - st_ss_cp;- ss_sp - st_cs_cp; ct_cp];
        E_WE = [ex_WE ey_WE ez_WE];
    end
    %----------------------------------------------------------------------
    
    %dE
    %----------------------------------------------------------------------
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
        dE_dzeta_W  = reshape(ct_methods.map_WEtoW_stackDim2(dE_dzeta_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD),3,[],ns);
        dE_dzeta_G  = ct_methods.map_WtoG(dE_dzeta_W,R_G_W);
        
        dE_dzeta_G_93 = reshape(dE_dzeta_G,9,3,ns);
        dE_dqa_G_9qans = bsxfun(@times,Ba_tr,dE_dzeta_G_93(:,yidx2,:));
        
        dE_dzeta_W_93 = reshape(dE_dzeta_W,9,3,ns);
        dE_dt_W_91ns = sum(bsxfun(@times,dzeta_a_dt_tr,dE_dzeta_W_93),2);
    end
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
