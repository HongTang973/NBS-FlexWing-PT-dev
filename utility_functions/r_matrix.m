function R = r_matrix(d,del)
%returns the rotation matrix R for rotation of magnitude 'del' radians about the
%arbitrary vector d

if nargin == 1
    del = norm(d);
    d = d./max(del,1e-8);
end

sd = sin(del);
cd = cos(del);

R11 = d(1)^2+(1-d(1)^2)*cd;
R22 = d(2)^2+(1-d(2)^2)*cd;
R33 = d(3)^2+(1-d(3)^2)*cd;

R12 = d(1)*d(2)*(1-cd)-d(3)*sd;
R21 = d(1)*d(2)*(1-cd)+d(3)*sd;

R13 = d(1)*d(3)*(1-cd)+d(2)*sd;
R31 = d(1)*d(3)*(1-cd)-d(2)*sd;

R23 = d(2)*d(3)*(1-cd)-d(1)*sd;
R32 = d(2)*d(3)*(1-cd)+d(1)*sd;

R = [ R11 R12 R13 ; ...
    R21 R22 R23 ; ...
    R31 R32 R33 ];

%equivalent to the following expression
%R = eye(3) + sin(del)*getSkewMat(d) + (1-cos(del))*getSkewMat(d)*getSkewMat(d); slower to execute!
end