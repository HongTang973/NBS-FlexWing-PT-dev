function Rinfo = rTransform_projection_to_Rmat(sweep_deg,dihedral_deg,alphaFlag,alpha_deg)
%takes sweep, dihedral, and alpha angles in degrees and returns global
%rotation matrix and equivalent on-axis twist

%two forms of alpha referencing may be used depending on 'alphaFlag'
%alphaFlag = 'alphaXZglobal': alpha refers to orientation of ex projection in global XZ frame
%alphaFlag = 'alphaXZdihedral'; alpha refers to orientation of ex projection in XZ' dihedral frame

%calculate <ey>
%--------------------------------------------------------------------------
ey = [...
    sind(sweep_deg);
    cosd(sweep_deg);
    tand(dihedral_deg)*cosd(sweep_deg)]./(1+tand(dihedral_deg)^2*cosd(sweep_deg)^2)^0.5;

%calculate <ex>
%--------------------------------------------------------------------------
if isequal(alphaFlag,'alphaXZglobal')
    alpha_xzg = alpha_deg;
    ex_temp = [cosd(alpha_xzg);0;-sind(alpha_xzg)]; %start with ex projection in XZ plane
    
    %ex.ey = 0    ==>    ex1*ey1 + ex2*ey2 + ex3*ey3 = 0
    % ==>    ey = (-ex1*ey1 - ex3*ey3)/ey2
    ex_temp(2) = (-ex_temp(1)*ey(1) - ex_temp(3)*ey(3))/ey(2); %create y component of ex to ensure orthogonality with ey
    ex = ex_temp./norm(ex_temp); %normalise to obtain ex vector

elseif isequal(alphaFlag,'alphaXZdihedral')
    alpha_xzd = alpha_deg;
    
    rmat_dihedral = [...%rotation matrix from global to dihedral frame
        1  0               0             ;
        0  cosd(dihedral_deg) -sind(dihedral_deg);
        0  sind(dihedral_deg)  cosd(dihedral_deg)];
    
    eyPrime = rmat_dihedral.'*ey; %solve same problem as above in dihedral frame
    
    exPrime_temp = [cosd(alpha_xzd);0;-sind(alpha_xzd)];
    
    exPrime_temp(2) = (-exPrime_temp(1)*eyPrime(1) - exPrime_temp(3)*eyPrime(3))/eyPrime(2);
    exPrime = exPrime_temp./norm(exPrime_temp);
    
    ex = rmat_dihedral*exPrime; %transform back to global frame to get final ex vector

end

%calculate <ez>
%--------------------------------------------------------------------------
ez = cross(ex,ey);

%package output rotation quantities
%--------------------------------------------------------------------------
Rinfo.rotationMatrix = [ex ey ez]; %Transformation Matrix
Rinfo.axialRotation = asind(dot(ey,cross([cosd(sweep_deg);-sind(sweep_deg);0],ex))); %Equivalent axial rotation about ey

end

