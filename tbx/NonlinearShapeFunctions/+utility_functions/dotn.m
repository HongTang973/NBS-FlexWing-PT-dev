function AB = dotn(A,B,dim)
%return the dot product of two matrices along the specified
%dimension 'dim'
AB = sum(A.*B,dim);
end
