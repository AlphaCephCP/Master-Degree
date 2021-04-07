% Computes the eigenvalues of the three-dimensional Euler equations
%
% Parameters:
% - Q       : vector of conserved quantities
% - nv      : normal vector
%
% Returns:
% - Lambda  : vector of eigenvalues
%
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/PDE.f90
%

function Lambda=EigenvaluesEuler(Q,nv)
    % import global variables
    global gamma; % specific heat ratio

    irho = 1/Q(1);
    u    = (Q(2:4)*nv')*irho;                          % normal velocity
    p    = (gamma-1)*(Q(5) - 0.5*sum(Q(2:4).^2)*irho); % fluid pressure
    c    = sqrt(gamma*p*irho);                         % sound speed

    % remark: the order of the eigenvalues does not matter; the 
    %         Riemann solver takes care of this
    Lambda = [u-c,u,u,u,u+c];
end
