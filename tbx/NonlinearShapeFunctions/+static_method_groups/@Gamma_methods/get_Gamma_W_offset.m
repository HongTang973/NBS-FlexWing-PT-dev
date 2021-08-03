function Gamma_offset = get_Gamma_W_offset(Gamma,E,massOffset_I)
Gamma_offset = Gamma + utility_functions.MultiProd_(E,massOffset_I);
end