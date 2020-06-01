function int_vec_out = integrate2_delxConst(x,y,fullValOnly,dim)
    %requires x to have odd number of elements (even intervals)
    %each pair of x intervals must have equal width
    %---int_y(x)_dx, returns cumulative output, simpsons quadratic form
    %dim is the dimension along which y is integrated
    %if fullValOnly is true then only the full integral value over the entire x range is returned
    %x: length (2n+1) vector, y: matrix with dimension 'dim' of length (2n+1)
    
    if nargin < 4
        dim = 2;
    end
    if nargin < 3
        fullValOnly = false;
    end
    
    nx = numel(x);
    szy = size(y);
    ndim_y = numel(szy);
    
    assert(nx == szy(dim),['x and y must have equal length along the specified dimension ''dim'' (dim = ' num2str(dim) ').'])
    
    idx = repmat({':'},1,ndim_y);
    idx1 = idx; idx1{dim} = 1:2:nx-2;
    idx2 = idx; idx2{dim} = 2:2:nx-1;
    idx3 = idx; idx3{dim} = 3:2:nx;
    idx_cumsum = idx; idx_cumsum{dim} = 2:nx;
    
    dx = x(2)-x(1);
    y1 = y(idx1{:});
    y2 = y(idx2{:});
    y3 = y(idx3{:});
    
    int_vec = zeros(szy);
    szdint = szy; szdint(dim) = szdint(dim)-1;
    dint = zeros(szdint);
    
    dint(idx1{:}) = (5*y1+8*y2-y3)*dx/12;
    dint(idx2{:}) = (-y1+8*y2+5*y3)*dx/12;
    
    int_vec(idx_cumsum{:}) = cumsum(dint,dim);
    
    if fullValOnly
        idx_out = idx; idx_out{dim} = nx;
        int_vec_out = int_vec(idx_out{:});
    else
        int_vec_out = int_vec;
    end

end
