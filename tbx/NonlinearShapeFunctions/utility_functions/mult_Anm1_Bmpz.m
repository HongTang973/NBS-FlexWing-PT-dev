function AB = mult_Anm1_Bmpz(A,B)
%multiplication of the 2D matrix A by each 2D slice (:,:,i) of matrix B
    szA = size(A);
    szB = size(B);
    AB = reshape(A*reshape(B,szB(1),prod(szB(2:end))),[szA(1) szB(2:end)]);
end

