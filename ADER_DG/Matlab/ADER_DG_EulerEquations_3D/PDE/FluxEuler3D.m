% Evaluates the flux vector
% 
% readable version available at 
% http://www.theoretical-physics.net/dev/fluid-dynamics/euler.html
% section 'Conservative Form of the Euler Equations'
%
% Parameters:
% - Q     : vector of conserved quantities
%           Q = [rho u v w E] (= [density, momentum, total energy])
% 
% Returns:
% - F    : a collection of column vectors, one flux vector per dimension
% 
% This code closely follows Michael Dumbser's code given 
% in ADER_DG_EulerEquations_3D/PDE.f90.
%


function F=FluxEuler3D(Q)
    % import global variables
    global gamma; % specific heat ratio

    irho = 1/Q(1);
    p = (gamma-1)*(Q(5) - 0.5*sum(Q(2:4).^2)*irho); % fluid pressure

    % flux (row) vectors
    % x
    Fx(1)   = Q(2);
    Fx(2:5) = Q(2:5);
    Fx(5)   = Fx(5)+p;
    Fx(2:5) = irho*Q(2).*Fx(2:5);
    Fx(2)   = Fx(2)+p;

    % y
    Fy(1)   = Q(3);
    Fy(2:5) = Q(2:5);
    Fy(5)   = Fy(5)+p;
    Fy(2:5) = irho*Q(3).*Fy(2:5);
    Fy(3)   = Fy(3)+p;

    % z
    Fz(1)   = Q(4);
    Fz(2:5) = Q(2:5);
    Fz(5)   = Fz(5)+p;
    Fz(2:5) = irho*Q(4).*Fz(2:5);
    Fz(4)   = Fz(4)+p;

    F = [Fx', Fy', Fz'];
    
end


