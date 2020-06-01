function T = T_matrix(d,del)

    %returns the tangency matrix T for rotation of magnitude 'del' radians about the
    %arbitrary vector d
    
    %omega = T*d/dt[del*d]
    
    %d = d./norm(d); %-- unit vector d must be supplied if commented out
    
    if del<=1e-6
        DELTAskew = [0 -d(3) d(2);d(3) 0 -d(1);-d(2) d(1) 0]*del;
        T = eye(3) - 1/2*DELTAskew + 1/6*DELTAskew*DELTAskew;
    else
    
    sdod = sin(del)/del;
    cdm1od = (cos(del)-1)/del;
    
    T11 = d(1)^2+(1-d(1)^2)*sdod;
    T22 = d(2)^2+(1-d(2)^2)*sdod;
    T33 = d(3)^2+(1-d(3)^2)*sdod;
    
    T12 = d(1)*d(2)*(1-sdod)-d(3)*cdm1od;
    T21 = d(1)*d(2)*(1-sdod)+d(3)*cdm1od;
    
    T13 = d(1)*d(3)*(1-sdod)+d(2)*cdm1od;
    T31 = d(1)*d(3)*(1-sdod)-d(2)*cdm1od;
    
    T23 = d(2)*d(3)*(1-sdod)-d(1)*cdm1od;
    T32 = d(2)*d(3)*(1-sdod)+d(1)*cdm1od;
    
    T = [ T11 T12 T13 ; ...
          T21 T22 T23 ; ...
          T31 T32 T33 ];
    
    %equivalent to the following expression
    %dskew = [0 -d(3) d(2);d(3) 0 -d(1);-d(2) d(1) 0];
    %T = eye(3) + (cos(del)-1)*dskew/del + (1-sin(del)/del)*dskew*dskew    slower to execute!

end