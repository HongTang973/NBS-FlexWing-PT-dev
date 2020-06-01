function bool = eq_tol(A,B,varargin)
%determines if two quantities are equal within a specified tolerance
%measure of distance between quantities is based upon 'norm()'
%if relTol is supplied then this acts on the first argument to the function

relTol = get_option(varargin,'relTol',inf);
absTol = get_option(varargin,'absTol',inf);

nmA = norm(A);
abs_normDiff = abs(nmA-norm(B));

bool = abs_normDiff < min(absTol,nmA*relTol);

end

