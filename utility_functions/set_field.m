function Container = set_field(Container,address,val)                                  %#ok
%sets the field 'O.address' to the value 'val' and returns 'O'
%need to edit this function to not return error if field does not exist
%i.e. something like QOI = [], QOI.Gamma_G = '-' should be allowed
if ~isempty(address) && address(1) ~= '.', address = ['.' address]; end
%        try
%            eval(['O' address ';'])
%        catch
%            error(['Static method return_field(O,''' address '''): O' address ' not a valid field address.']);
%        end

eval(['Container' address '=val;'])

end

