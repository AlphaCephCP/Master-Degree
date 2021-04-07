% Compute the Gauss-Lobatto nodes and weights on the
% interval [a,b]. N is the number of quadrature nodes. 
%
% In contrast to the Gauss-Legendre nodes, the Gauss-Lobatto
% nodes include the end points of the interval. The 
% remaining (N-2) nodes are chosen in such a way that the 
% quadrature accuracy is maximised.
%
% The script is based on:
% [1] http://www.scientificpython.net/pyblog/lobatto-quadrature (primarily)
% [2] Gubner, Gaussian Quadrature and the Eigenvalue Problem,
%     available at http://gubner.ece.wisc.edu/gaussquad.pdf, 2014
%
% Last updated: Nov 11 2015
%

function [x,w]=gaussLobatto(N,a,b)
    N=N-1;

    % n-the Cartesian basis vector
    eN = [zeros(1,N-1) 1]';

    % compute three-term Legendre recurrence coefficients
    % see [2] example 15
    % phi_n+1(x) = (x-alpha_n) phi_n(x) - beta_n phi_{n-1}(x)
    % where (i)  alpha_n = 0
    %       (ii) beta_n  = n^2/(4n^2-1)
    alpha  = zeros(1,N);
    beta   = [2 arrayfun(@(n)n^2/(4*n^2-1), 1:N-1)];


    % create tridiagonal Jacobi matrix, see [2] Theorem 12
    gamma = sqrt(beta);
    J = diag(alpha,0) + diag(gamma(2:N),1) + diag(gamma(2:N),-1);


    % See [1] on how to reduce the Gauss-Lobatto problem 
    % to the Gauss-Legendre problem. We need to solve 
    % 3 systems to determine (a) omega_1, ..., omega_{N+1} and
    %                        (b) gamma_1, ..., gamma_N
    % Nomenclature as in [1].

    % system 1
    J1    = J - a*eye(N);
    delta = J1\eN;

    % system 2
    J2 = J - b*eye(N);
    mu = J2\eN;

    % alpha_{N+1} and gamma_n are still missing
    % therefore solve system 3
    M   = [[ 1   -delta(N) ]
           [ 1   -mu(N)    ]];
    rhs = [a b]'; 
    soln = M\rhs; % soln = [ alpha_{N+1}, beta^2_N ]
    
    
    % now we are back to the Gauss-Legendre (GL) eigenvalue problem in [-1,1]
    % set up the tridiagonal Jacobi matrix...
    gammaGL = [gamma sqrt(soln(2))];
    alphaGL = [alpha soln(1)];
    T = diag(alphaGL,0) + diag(gammaGL(2:end),1) + diag(gammaGL(2:end),-1);

    % ... and solve the eigenvalue problem, see [2]
    [V,D] = eig(T);
    x = diag(D);      % Lobatto nodes
    [x,i] = sort(x);
    w = 2*V(1,i).^2;  % Lobatto weights


    % map nodes and weights onto interval [a,b]
    x = (b-a)/2 .* x + (a+b)/2;
    w = (b-a)/2 .* w;

end

