clc;
propCell = properties(P2);
for i_ = 1:numel(propCell)
    prop = propCell{i_};
    if ~isequal(P2.(prop),P3.(prop))
        disp(prop);
    end
end