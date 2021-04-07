% Assembles the mass matrix and its inverse.
%
% Parameters:
% - N    : degree of the approximation polynomial
% - xGPN : Gaussian quadrature points, typically the Gauss-Legendre nodes
% - wGPN : Gaussian weights
%
% Returns:
% - M    : the mass matrix of size (N+1)x(N+1)
% - invM : the inverse of the mass matrix 
%
% Extracted from Michael Dumbser's code given in ADERDG3D/Init.f90
%

function [M,iM]=massMatrix(N,xGPN,wGPN)
    M = zeros(N+1,N+1);

    for i=1:N+1
        xiGPN = xGPN(i); % i-th Gaussian quadrature point ...
        wiGPN = wGPN(i); % ... and the corresponding weight
        [phi,phi_xi] = BaseFunc1D(N,xGPN,xiGPN); % Lagrange basis
        
        % assemble mass matrix
        for k=1:N+1
            for j=1:N+1
                M(k,j)=M(k,j)+wiGPN*phi(k)*phi(j);
            end
        end
    end

    % for the Gauss-Legendre nodes this yields
    % M = diag(wGPN)

    % inverse of the mass matrix
    iM=inv(M);
end
