% Converts primitive variables into conservative form
% for the Euler equations.
%
% Parameters:
% - V : vector of primitive variables,
%       V = [rho u1 u2 u3 p]
%
% Returns:
% - Q : vector of conserved quantities    
%
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/PDE.f90
%

function Q = Prim2ConsEuler(V)
    % import global variables
    global gamma; % specific heat ratio
 
    rho = V(1);
 
    % Q = [rho u v w E]   
    Q(1)   = rho; 
    Q(2:4) = rho*V(2:4);
    Q(5)   = V(5)/(gamma-1) + 0.5*rho*(sum(V(2:4).^2));
end
