function AB = mult_Anmz_Bmp1(A,B)
%multiplication of each 2D slice (:,:,i) of matrix A by the 2D matrix B
    B_ = permute(A,[2 1 3]);
    A_ = B.';
    szA_ = size(A_);
    szB_ = size(B_);
    if length(szB_)==3, szB_3=szB_(3); else, szB_3=1; end
    AB_ = reshape(A_*reshape(B_,szB_(1),szB_(2)*szB_3,1),[szA_(1) szB_(2) szB_3]);
    AB = permute(AB_,[2 1 3]);
end

