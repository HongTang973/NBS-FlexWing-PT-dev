function [ClosestValue,index] = closestValue(Array,Values)
    %return the indices and values of [Array] that are closest to the values in [Values]

    ColVec = Array(:);
    RowVec_request = reshape(Values,1,[]);
    [~,index] = min(abs(bsxfun(@plus,ColVec,-RowVec_request)));
    ClosestValue = Array(index);
end

