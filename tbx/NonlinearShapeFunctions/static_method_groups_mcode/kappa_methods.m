classdef kappa_methods < handle
    %A collection of methods responsible for the calculation of varTheta and it's derivatives
    
    methods (Static = true)
        
        %method to return group 1 varTheta variables
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
            KAPPA_I = MultiProd_(TD_tr, mult_Anm1_Bmpz(R_W_WE,KAPPA_IE) );
            %-----------------------------------
            
            [dKAPPA_dt_I,dKAPPA_dqa_I_tr_Dimnqax3xns] = kappa_methods.get_dKAPPA_I(st,ct,sp,cp,dth_ds,dsi_ds,nszrs,R_W_WE_tr,B_th,B_si,B_ph,dB_th,dB_si,dB_ph,dzeta_a_dt_tr,d2th_dsdt,d2si_dsdt,d2ph_dsdt,st_sp,st_cp,ct_sp,ct_cp,yidx2,TD);
            dKAPPA_dqa_I_tr_Dim1x3xnsxnqa = permute(dKAPPA_dqa_I_tr_Dimnqax3xns,[4 2 3 1]);
            
        end
        
    end
    
    
    methods (Static = true, Access = private) %Group 1 Static Methods
        
        %dKAPPA_I
        %----------------------------------------------------------------------
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
            dxi_dzeta_I = MultiProd_(mult_Anmz_Bmp1(dxi_dzeta_IE,R_W_WE_tr),TD);
            dxi_dzetads_I = MultiProd_(mult_Anmz_Bmp1(dxi_dzetads_IE,R_W_WE_tr),TD);
            %---------------------------------------
            
            dKAPPA_dq_I_tr = bsxfun(@times,[B_th;B_si;B_ph],dxi_dzeta_I(yidx2,:,:))+bsxfun(@times,[dB_th;dB_si;dB_ph],dxi_dzetads_I(yidx2,:,:));
            
            dKAPPA_dt_I_tr = MultiProd_(dzeta_a_dt_tr,dxi_dzeta_I)+...
                MultiProd_([d2th_dsdt,d2si_dsdt,d2ph_dsdt],dxi_dzetads_I);
            dKAPPA_dt_I = reshape(dKAPPA_dt_I_tr,3,1,[]);

        end
        
%         function val = spinCross(vec11,vec12,vec21,vec22,vec31,vec32)
%             val = 1/2*(MultiProd_(getSkewMat(vec11),vec12)+MultiProd_(getSkewMat(vec21),vec22)+MultiProd_(getSkewMat(vec31),vec32));
%         end
        
    end
    
end

