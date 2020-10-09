function option_value = get_option(VARARGIN,option_flag,default_value)
%checks for optional_flag among the input arguments in VARARGIN
%if flag is present, return the following argument as its value
%if flag is not present, return the default value
[member,index] = find(strcmpi(option_flag,VARARGIN));
if member
    option_value = VARARGIN{index+1};
else
    option_value = default_value;
end
end
