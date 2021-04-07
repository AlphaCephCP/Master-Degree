% Converts conserved quantities into primitive variables
% for the Euler equations.
%
% Parameters:
% - Q : vector of conserved quantities,
%       Q = [rho u v w E]
%
% Returns:
% - V : vector of primitive variables
%   
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/PDE.f90
%

function V = Cons2PrimEuler(Q)
    % import global variables
    global gamma; % specific heat ratio

    irho = 1/Q(1);
    p = (gamma-1)*(Q(5) - 0.5*sum(Q(2:4).^2)*irho); % fluid pressure

    % V = [rho u1 u2 u3 p]
    V(1)   = Q(1);
    V(2:4) = irho*Q(2:4);
    V(5)   = p;
end
