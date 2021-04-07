% Assembles the stiffness matrix.
%
% Parameters:
% - N    : degree of the approximation polynomial
% - xGPN : Gaussian quadrature points, typically the Gauss-Legendre nodes
% - wGPN : Gaussian weights
%
% Returns:
% - K   : the stiffness matrix
%
% Extracted from Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/Init.f90
%

function K=stiffnessMatrix(N,xGPN,wGPN)
    K = zeros(N+1,N+1);

    for i=1:N+1
        xiGPN = xGPN(i); % i-th Gaussian quadrature point ...
        wiGPN = wGPN(i); % ... and the corresponding weight
        [phi,phi_xi] = BaseFunc1D(N,xGPN,xiGPN); % Lagrange basis
        
        % assemble stiffness matrix
        for k=1:N+1
            for j=1:N+1
                K(k,j)=K(k,j)+wiGPN*phi_xi(k)*phi(j);
            end
        end
    end
end
