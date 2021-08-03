function dGamma_dt_offset = get_dGamma_dt_offset(dGamma_dt,dE_dt,massOffset_I)
dGamma_dt_offset = dGamma_dt + utility_functions.MultiProd_(dE_dt,massOffset_I);
end