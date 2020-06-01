function AB = MultiProd_(A,B,~)
%returns the product of the multidimensional arrays A and B
%product obtained by multiplying together all of the (dim1,dim2) 2d sub-matrices of A and B
%automatic array expansion of unit length dimensions
%examples
%    MultiProd_(rand(2,3,4),rand(3,1,4))
%    MultiProd_(rand(2,3,1),rand(3,1,4))
%    MultiProd_(rand(2,3,1,2),rand(3,1,4,2))

if isempty(A)||isempty(B)
    
    %handling of special 0 length dimension case
    %warning: no checks performed on validity of A/B size
    szA = size(A); szB = size(B);
    nszA = numel(szA); nszB = numel(szB);
    ndim = max(nszA,nszB);
    szA(nszA+1:ndim+1) = 1; szB(nszB+1:ndim+1) = 1;
    ABdim = [szA(1),szB(2),min(szA(3:end),szB(3:end))];
    AB = zeros(ABdim);
    
else
    
    %if the mtimesx.mex exists then use it
    %else use multiprod.m
    global mult3d_mex
    if mult3d_mex
        AB = mtimesx(A,B);
    else
        AB = multiprod(A,B);
    end
    
end

end
