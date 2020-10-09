function Rinfo = rTransform_Rmat_to_projection(Rmat)
%takes global rotation matrix and returns sweep, dihedral and alpha angles
%in degrees
%alpha is returned in both the global and dihedral systems such that:
    %Rinfo.alphaXZglobal': alpha refers to orientation of ex projection in global XZ frame
    %Rinfo.alphaXZdihedral'; alpha refers to orientation of ex projection in XZ' dihedral frame

ex = Rmat(:,1);
ey = Rmat(:,2);
ez = Rmat(:,3);

%calculate <sweep>
%--------------------------------------------------------------------------
planeNorm = [0;0;-1];
Rinfo.sweep_deg = get_RotationInPlane(planeNorm,[0;1;0],ey)*180/pi;

%calculate <dihedral>
%--------------------------------------------------------------------------
planeNorm = [1;0;0];
Rinfo.dihedral_deg = get_RotationInPlane(planeNorm,[0;1;0],ey)*180/pi;

%calculate <alphaXZglobal>
%--------------------------------------------------------------------------
planeNorm = [0;1;0];
Rinfo.alphaXZglobal_deg = get_RotationInPlane(planeNorm,[1;0;0],ex)*180/pi;

%calculate <alphaXZdihedral>
%--------------------------------------------------------------------------
planeNorm = [0;cosd(Rinfo.dihedral_deg);sind(Rinfo.dihedral_deg)];
Rinfo.alphaXZdihedral_deg = get_RotationInPlane(planeNorm,[1;0;0],ex)*180/pi;

end

