function val = crossndim(a,b,dim)

    %returns the cross product of the multidimensional arrays 'a' and 'b'
    %performed along the dimension 'dim'

    sza = size(a); ndim_a = numel(sza);
    idx1 = repmat({':'},1,ndim_a); idx1{dim} = 1;
    idx2 = repmat({':'},1,ndim_a); idx2{dim} = 2;
    idx3 = repmat({':'},1,ndim_a); idx3{dim} = 3;
    
    val = zeros(sza, 'like', b);
    val(idx1{:}) = a(idx2{:}).*b(idx3{:})-a(idx3{:}).*b(idx2{:});
    val(idx2{:}) = a(idx3{:}).*b(idx1{:})-a(idx1{:}).*b(idx3{:});
    val(idx3{:}) = a(idx1{:}).*b(idx2{:})-a(idx2{:}).*b(idx1{:});
    
end