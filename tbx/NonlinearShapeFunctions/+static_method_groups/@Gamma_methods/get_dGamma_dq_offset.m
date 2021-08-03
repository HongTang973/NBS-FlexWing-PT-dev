function dGamma_dq_offset = get_dGamma_dq_offset(dGamma_dq,dE_dq,massOffset_I)
dGamma_dq_offset = dGamma_dq + utility_functions.MultiProd_(dE_dq,massOffset_I);
end