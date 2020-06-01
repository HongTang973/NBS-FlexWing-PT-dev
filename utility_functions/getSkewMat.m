function skewMat = getSkewMat(a)

%return the skew symmetic matrix corresponding to the cross product
%operation of the vector 'a'
zrs(1,1,size(a,3)) = 0;
skewMat = [zrs -a(3,1,:) a(2,1,:);a(3,1,:) zrs -a(1,1,:);-a(2,1,:) a(1,1,:) zrs];

end
