function data_new = interp1dim(data,s_old,s_new,dim,varargin)
    %resample data along the indicated dimension 'dim'
    
    method = get_option(varargin,'method','linear'); %choose from 'nearest', 'next', 'previous', 'linear','spline','pchip', 'makima', or 'cubic'
    extrapolation = get_option(varargin,'extrapolation','extrap'); %choose 'extrap' or provide a scalar value
    
    permuteVec = [1 2 3];
    permuteVec(1) = dim;
    permuteVec(dim) = 1;
    
    data_permute = permute(data,permuteVec);
    
    data_new_permute = interp1(s_old(:),data_permute,s_new(:),method,extrapolation);
    
    data_new = permute(data_new_permute,permuteVec);
end