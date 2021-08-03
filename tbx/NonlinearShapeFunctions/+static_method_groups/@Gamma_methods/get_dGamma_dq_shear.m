function dGamma_dq = get_dGamma_dq_shear(Gamma_Integration_Function,del_s,ex,ey,ez,dex_dq,dey_dq,dez_dq,tau_x,tau_y,tau_z,dtau_x_dq,dtau_y_dq,dtau_z_dq,root_idx,dGamma_dq_root_G)

if ~isempty(dex_dq) && ~isempty(dey_dq) && ~isempty(dez_dq)
    dGamma_dq_Integrand_pt1 = ...
        bsxfun(@times,dex_dq,tau_x)+...
        bsxfun(@times,dey_dq,1+tau_y)+...
        bsxfun(@times,dez_dq,tau_z);
else
    dGamma_dq_Integrand_pt1 = 0;
end

if ~isempty(dtau_x_dq) && ~isempty(dtau_y_dq) && ~isempty(dtau_z_dq)
    dGamma_dq_Integrand_pt2 = ...
        bsxfun(@times,ex,dtau_x_dq)+...
        bsxfun(@times,ey,dtau_y_dq)+...
        bsxfun(@times,ez,dtau_z_dq);
else
    dGamma_dq_Integrand_pt2 = 0;
end

dGamma_dq_Integrand = dGamma_dq_Integrand_pt1 + dGamma_dq_Integrand_pt2;

dGamma_dq = Gamma_Integration_Function(del_s,dGamma_dq_Integrand,root_idx) + dGamma_dq_root_G;
end