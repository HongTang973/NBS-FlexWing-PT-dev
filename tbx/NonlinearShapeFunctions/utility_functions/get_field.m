function fieldVal = get_field(Container,address)                               %#ok
%returns the field value contained at 'O.address'

if ~isempty(address) && address(1) ~= '.', address = ['.' address]; end
try
    eval(['Container' address ';'])
catch
    error(['get_field method return_field(arg,''' address '''): arg' address ' not a valid field address.']);
end

eval(['fieldVal = Container' address ';'])

end

