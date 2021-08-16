function [dQaero,Qaero,Fqc,Mqc,Drag,alpha] = aero_stripTheory_usteady_LeishmanIndicial(Qaero,rho,Vinf,V3qrt,xAp,yAp,zAp,Omega,chord,width,AIC,C_D0,aeroCoeff2D,qsteady)
%% coded implementation of Leishman's Indicial Response Method  - C.Howcroft
% [use publish button to view latex comments]

% INPUT ARGUMENTS
% Qaero [nxv] - aerodynamic state vector
% rho [1] - air density
% Vinf [3x1] - free stream flow vector (magnitude = flight speed)
% V3qrt [3xv] - velocity of the three quarter chord point
% xAp [3xv] - unit vector in the strip chordwise direction from leading to trailing edge
% yAp [3xv] - unit vector perpendicular to xAp and normal to the aerofoil cross sectional plane, pitch up is positive in this direction
% zAp [3xv] - out of plane unit vector perpendicular to xAp and yAp
% Omega [3xv] - rotational velocity of the aerodynamic strip
% chord [1xv] - chord length of the strip in the xAp direction
% width [1xv] - width of the strip in the yAp direction
% AIC [1xv] - aerodynamic influence coefficients by which the aerodynamic forces are scaled
% C_D0 [1xv] - coefficient of parasitic drag, calculated externally by empirical function and supplied by user (set C_D0=0 to neglect this effect)
% qsteady [true/(false or unspecified)] - quasi steady flag, set to true to bypass unsteady terms

% Notes:
% - a single arbitrary inertial reference coordinate system must be selected in which to specify the [3x...] input quantities
%   output quantities are returned in this same reference system
% - this aerodynamic function is vectorised by stacking multiple inputs along dimension v
%   v may be along any dimension except the 1st (dim v = 2 in these example comments)
%   the output is vectorised along the same dimension

% OUTPUT
% dQaero [nxv] - time derivative of aerodynamic state vector
% Qaero [nxv] - if qsteady is true then Qaero is assigned the quasi steady state value, else Qaero is simply passed through the function unchanged
% Fqc [3xv] - aerodynamic force vector generated at the 1/4 chord location excluding parasitic drag
% Mqc [3xv] - aerodynamic moment vector generated at the 1/4 chord location
% Drag [3xv] - parasitic drag vector

%% Define Aerodynamic Parameters
%Calculate apparent flow vector at 3/4 chord
Vflow = Vinf-V3qrt;
velocity = sum(Vflow.^2,1).^0.5;

%Take components of this apparent flow in the strip coordinate system
vx3qrt = sum(Vflow.*xAp,1);
%vy3qrt = Vflow.'*yAp; Not Used
vz3qrt = sum(Vflow.*zAp,1);

%neglect linear and angular acceleration contributions to non-circulatory terms
dvzmid_dt = 0;
d2alpha_dt2 = 0;

% Dynamic Pressure
Pdyn = 0.5*rho*velocity.^2;


%Angle of attack of strip
alpha = atan(vz3qrt./vx3qrt);
% alpha = -atan2(vz3qrt,vx3qrt);

% R_Ap_A = permute([xAp, yAp, zAp],[2 1 3]);
% X_A = [1;0;0];
% X_Ap = MultiProd_(R_Ap_A,X_A);
% phi = atan(X_Ap(3,:,:)./ X_Ap(1,:,:));
% phi_ = squeeze(phi*180/pi);

%d[alpha]/dt
dalpha_dt = sum(Omega.*yAp,1);

%Half chord
b = chord/2;

%% Define Wagners Function
%%%
% Exponential approximation of Wagners Function psi(t)
%
% $$\psi(t) = 1 - \sum_{i=1}^{N} a_i e^{-b_i s}\;\;, \qquad s = \frac{v_x t}{b}$$
% for R.T. Jones approximation
aiCoeffs = [0.165 ; 0.335];
biCoeffs = [0.041 ; 0.32 ];

%initial value of Wagners Function
psi_0 = 1 - sum(aiCoeffs);

%% Aerodynamic State ODE
%%%
% Bypass unsteady terms if qsteady flag is 'true'
%
% $$\dot\xi_i = 0\;\;,\qquad$$
% $$\xi_i = \frac{b\alpha}{b_i v_x}$$
if nargin > 12 && qsteady
    dalpha_dt = 0;
    Qaero = b.*alpha./vx3qrt.*(1./biCoeffs);
    dQaero = Qaero*0;
else
    aeroCoeff2D = [];
    %%%
    % else get state derivative from ode
    %
    % $$\dot\xi_i + b_i\frac{v_x}{b}\xi_i = \alpha$$
    dQaero = alpha - (vx3qrt./b.*biCoeffs.*Qaero);
end

if isempty(aeroCoeff2D)
    %% Lift Coefficient
    %%%
    % Circulatory Component
    %
    % $$C_{Lc} = 2\pi\alpha\psi(0) + 2\pi\frac{v_x}{b}\sum_{i=1}^{N} a_i b_i \xi_i$$
    C_Lc = 2*pi*alpha*psi_0 + 2*pi*vx3qrt./b.*sum(aiCoeffs.*biCoeffs.*Qaero,1);

    %%%
    % Non-circulatory Component
    %
    % $$C_{Li} = \pi\frac{b}{v_x^2}( v_x \dot{\alpha} - \dot{v}_{zMid})$$
    C_Li = pi*b./vx3qrt.^2.*( vx3qrt.*dalpha_dt - dvzmid_dt );

    %%%
    % Add circulatory and non-circulatory components to obtain lift coefficient
    CL = C_Lc + C_Li;

    %% Moment Coefficient
    %%%
    % Circulatory Component
    C_Mc = 0; %zero when taken at quarter chord

    %%%
    % Non-circulatory Component
    %
    % $$C_{Li} = \pi\frac{b^2}{v_x^2}(-v_x \dot{\alpha} + \frac{1}{2}( \dot{v}_{zMid} - \frac{b}{4}\ddot{\alpha} ))$$
    C_Mi = pi*b.^2./vx3qrt.^2.*(-vx3qrt.*dalpha_dt + 0.5*(dvzmid_dt - b/4*d2alpha_dt2) );

    %%%
    % Add circulatory and non-circulatory components to obtain lift coefficient
    CM = C_Mc + C_Mi;

    %% Create Aerodynamic Force Vectors
    %%%
    % Lift
    Fqc  = Pdyn.*chord.*width.*AIC.*CL.*zAp; %Include circulatory component of drag by casting this aero force normal to the strip (ie. in the zAp direction).

    %%%
    % Moment
    Mqc  = Pdyn.*chord.*width.*AIC.*CM.*yAp;

    %%%
    % Parasitic Drag
    Drag = Pdyn.*chord.*width.*AIC.*C_D0.*xAp;

else

    [CL, CD, CM] = interpCoeff(alpha, aeroCoeff2D);

    CL = permute(CL,[3 2 1]);
    CD = permute(CD,[3 2 1]);
    CM = permute(CM,[3 2 1]);

    ca = cos(alpha);
    sa = sin(alpha);

    LIFT = Pdyn.*chord.*width.*AIC.*CL;
    DRAG = Pdyn.*chord.*width.*AIC.*CD;

    Fqc  = LIFT.*ca + DRAG.*sa.*zAp;
    Mqc  = Pdyn.*chord.*width.*AIC.*CM.*yAp;
    Drag = DRAG.*ca - LIFT.*sa.*xAp;

end

    function [cl, cd, cm, alpha] = interpCoeff(alpha, aeroCoeff2D)

        alpha_in_file = aeroCoeff2D(:,1,:);
        cl_in_file = aeroCoeff2D(:,2,:);
        cd_in_file = aeroCoeff2D(:,3,:);
        cm_in_file = aeroCoeff2D(:,4,:);

        Nseg = length(alpha);
        IDSeg = (1:Nseg)';
        Nrow_cl = size(cl_in_file,1);
        alpha = squeeze(alpha*180/pi);

        % --- Interpolate_LiftDrag
        L = alpha_in_file(2)-alpha_in_file(1);
        if all(isreal(alpha))
            Id          = ceil(alpha/L - alpha_in_file(1)/L);
            LinId      = Id + (IDSeg-1)*Nrow_cl; % Id + (find(NcId)-1)*Nrow_cl; % LinId = sub2ind(size(cl_in_file),Id,find(NcId))
            Gamma   = (alpha-alpha_in_file(Id))/L;
        else
            % for complex difference step
            %                                     aoa_EqC   = real(aoa) + imag(aoa);
            Id            = ceil(real(alpha)/L - alpha_in_file(1)/L); %ceil(aoa_EqC/L - alfa_in_file(1)/L);
            LinId        = Id + (IDSeg-1)*Nrow_cl;
            Gamma     = (alpha-alpha_in_file(Id))/L; %(real(aoa(NcId,ib))-alfa_in_file(Id))/L + (aoa(NcId,ib))/L;
        end

        cl  = (1-Gamma).*cl_in_file(LinId) + Gamma.*cl_in_file(LinId+1);
        cd = (1-Gamma).*cd_in_file(LinId) + Gamma.*cd_in_file(LinId+1);
        cm = (1-Gamma).*cm_in_file(LinId) + Gamma.*cm_in_file(LinId+1);  % positive moment = nose up convention
    end

end
