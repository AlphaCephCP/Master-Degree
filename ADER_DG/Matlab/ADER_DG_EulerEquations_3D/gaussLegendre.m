% Compute the Gauss-Legendre nodes and weights on the 
% interval [a,b]. N is the number of quadrature nodes.
%
% Parameters:
% - N   : number of Gauss-Legendre nodes
% - a   : left interval endpoint
% - b   : right interval endpoint
%
% Returns:
% - x   : Gauss-Legendre nodes in [a,b]
% - w   : weights
%
%
% The script is based on:
% [1] Gubner, Gaussian Quadrature and the Eigenvalue Problem, 2014,
%     available at http://gubner.ece.wisc.edu/gaussquad.pdf
%

function [x,w]=gaussLegendre(N,a,b)
    % compute three-term Legendre recurrence coefficients
    % see [1] example 15
    % phi_n+1(x) = (x-alpha_n) phi_n(x) - beta_n phi_{n-1}(x)
    % where (i)  alpha_n = 0
    %       (ii) beta_n  = 1/(4-n^(-2))
    alpha  = zeros(1,N);
    beta   = [arrayfun(@(n)1/(4-n^(-2)), 1:N-1)];

    % compute roots of Legendre polynomials in [-1,1] from
    % an eigenvalue decomposition of the symmetric, tridiagonal
    % Jacobi matrix, see [1] Theorem 12
    J = diag(alpha,0) + diag(sqrt(beta),-1) + diag(sqrt(beta),1);
    [V,D] = eig(J);
    x = diag(D);      % Gauss-Legendre nodes
    [x,i] = sort(x);
    w = 2*V(1,i).^2;  % Gauss-Legendre weights
     
    % map nodes and weights onto the interval [a,b]
    x = (b-a)/2 .* x + (a+b)/2;
    w = (b-a)/2 .* w;
end
