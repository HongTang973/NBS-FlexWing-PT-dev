function dGamma_dq = get_dGamma_dq(Gamma_Integration_Function,del_s,dey_dq,root_idx,dGamma_dq_root_G)
dGamma_dq = bsxfun(@plus, Gamma_Integration_Function(del_s,dey_dq,root_idx) , dGamma_dq_root_G);
end