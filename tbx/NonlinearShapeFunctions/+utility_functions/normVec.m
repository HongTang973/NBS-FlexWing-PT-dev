function [normVecs,mag] = normVec(v,dim)

%---normalise vectors along the dimension 'dim'
mag = sum((v.^2),dim).^0.5;
normVecs = v.*(1./max(mag,1e-10));

end
