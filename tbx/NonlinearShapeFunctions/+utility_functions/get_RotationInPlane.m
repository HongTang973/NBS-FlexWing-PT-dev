function theta = get_RotationInPlane(planeNorm,Vector1,Vector2)

    %Calculates the rotation from 'Vector1' to 'Vector2' in the plane
    %defined by the vector normal 'planeNorm'
    
    %theta is returned in the range -pi < theta <= pi

    planeNorm = planeNorm(:);
    Vector1 = Vector1(:);
    Vector2 = Vector2(:);

    uVec1_Prj = utility_functions.normVec( Vector1 - utility_functions.dotn(Vector1,planeNorm,1)*planeNorm , 1);
    uVec2_Prj = utility_functions.normVec( Vector2 - utility_functions.dotn(Vector2,planeNorm,1)*planeNorm , 1);
    
    signArg = sign(utility_functions.dotn(utility_functions.crossn(uVec1_Prj,uVec2_Prj,1),planeNorm,1));
    if signArg == 0, signArg = 1; end
    
    theta = acos(utility_functions.dotn(uVec1_Prj,uVec2_Prj,1))*signArg;

end

