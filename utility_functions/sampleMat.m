    function sampleMat = sampleMat(x_old,x_new)
        %generate linear sampling matrix from x_old to x_new system
        %e.g. wish to resample matrix A
        %A with rows as functions f(x_old)
        %returns sampleMat such that:
        %Asample = A*sampleMat
        %where: rows of Asample : f(x_new)
        
        Lx_old = length(x_old);
        Lx_new = length(x_new);
        sampleMat = sparse(Lx_old,Lx_new);
        
        xnew_idx = 1;
        xold_idx = 2;
        while (xold_idx <= Lx_old) && (xnew_idx <= Lx_new)
            if x_new(xnew_idx) <= x_old(xold_idx)
                frac = (x_old(xold_idx)-x_new(xnew_idx))/(x_old(xold_idx)-x_old(xold_idx-1));
                sampleMat(xold_idx-1:xold_idx,xnew_idx) = [frac;1-frac];
                %Asample(:,xnew_idx) = A(:,xold_idx-1)*frac+A(:,xold_idx)*(1-frac);
                xnew_idx = xnew_idx+1;
            else
                xold_idx = xold_idx+1;
            end
        end
    end