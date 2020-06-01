function int_val = int_midpointShooting(del_x,y,root_idx)
    %rectangular mid-point integration along the 3rd dimension
    %del_x: step size vector TODO review this function for non-uniform del_x
    %root_idx: index of root point, integration zero about this point
    int_val = y; %allocate space
    y(:,:,3:2:end,:) = y(:,:,2:2:end,:);
    y_ = y(:,:,2:end,:);
    int_val(:,:,2:end,:) = cumsum(bsxfun(@times,del_x,y_),3);
    int_val(:,:,1,:) = y(:,:,1,:)*0;
    if root_idx~=1
        int_val = bsxfun(@plus, int_val, -int_val(:,:,root_idx,:));
    end
end
