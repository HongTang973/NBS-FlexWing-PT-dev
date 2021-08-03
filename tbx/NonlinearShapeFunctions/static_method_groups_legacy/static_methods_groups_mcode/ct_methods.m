classdef ct_methods < handle
    %coordinate transform methods
    
methods (Static = true)
        
    function M_W = map_WEtoW(M_WE,Rs_W_WE_tr,Rv_W_WE_tr,TD)
        %M_W = MultiProd_(R_W_WE,MultiProd_(M_WE,R_W_WE_tr,[1 2]),[1 2]);%(3)x(3)x(ns)
        
        %if the rotation R_W_WE is modulo 90 degrees use following code
        %much quicker than MultiProd_ line
        smw = Rs_W_WE_tr; pmw = Rv_W_WE_tr;
        M_W = M_WE;
        for ii_ = 1:3
            for jj_ = 1:3
                M_W(ii_,jj_,:) = smw(ii_)*smw(jj_)*M_WE(pmw(ii_),pmw(jj_),:);
            end
        end
        if ~isempty(TD)
            M_W = MultiProd_(M_W,TD);
        end
    end
    
    function M_G = map_WtoG(M_W,R_G_W)
        M_G = mult_Anm1_Bmpz(R_G_W,M_W);
    end

    function M_A = map_WtoA(M_W,R_A_W)
        M_A = mult_Anm1_Bmpz(R_A_W,M_W);
    end

    function M_G = map_AtoG(M_A,R_G_A)
        M_G = mult_Anm1_Bmpz(R_G_A,M_A);
    end
    
    function M_W_flat = map_WEtoW_stackDim2(M_WE_flat,Rs_W_WE_tr,Rv_W_WE_tr,TD)
        %mapping performing the transform M_W = R_W_WE*M_WE*R_W_WE_tr
        %where M_WE consists of a group of 9x1xns flattened matrices concatinated along dim 2
        %rotation R_W_WE must consist of combined 90 degree rotations
        %    i.e. all entries of R_W_WE are 1 or -1
        
        smw = Rs_W_WE_tr; pmw = Rv_W_WE_tr; ns = size(M_WE_flat,3);
        M_W_flat = M_WE_flat;
        for ii_ = 1:3
            for jj_ = 1:3
                M_W_flat(ii_+3*(jj_-1),:,:) = smw(ii_)*smw(jj_)*M_WE_flat(pmw(ii_)+3*(pmw(jj_)-1),:,:);
            end
        end

        if ~isempty(TD)
            M_W = reshape(M_W_flat,3,[],ns);
            for I = 1:size(M_W_flat,2)
                M_W(:,(1:3)+3*(I-1),:) = MultiProd_(M_W(:,(1:3)+3*(I-1),:),TD);
            end
            M_W_flat = reshape(M_W,size(M_W_flat));
        end
    end
        
end
    
end
