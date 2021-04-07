% Computes the Lagrange basis polynomial phi
% and its derivative phi_xi at xi.
%
% Parameters:
% - N      : degree of the interpolation polynomial
% - xin    : interpolation points (typically Gauss-Legendre nodes)
% - xi     : coordinate in [0,1] where to evaluate the basis 
%
% Returns:
% - phi    : Lagrange basis functions 
% - phi_xi : derivative of phi w.r.t. xi
%
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/DGBasis.f90
%

function [phi,phi_xi]=BaseFunc1D(N,xin,xi)
    phi     = ones(1,N+1);  % collection of Lagrange basis functions evaluated at xi
    phi_xi  = zeros(1,N+1); % collection of corresponding derivatives

    % Lagrange basis functions
    % l_m(x)  = \prod_{0\leq j \leq N+1 ,j\neq m} \frac{x-x_j}{x_m - x_j}
    % phi(xi) = [ l_0(xi) l_1(xi) ... l_{N+1}(xi) ]
    for m=1:N+1
        for j=1:N+1
            if j==m
                continue
            end
            phi(m) = phi(m) *(xi-xin(j))/(xin(m)-xin(j));
        end
    end
    

    % Derivatives of the Lagrange basis functions
    % l_m'(x) = d/dx (l_m(x)) 
    %         = \sum_{i \neq m} \frac{1}{x_m-x_i} \prod_{j \neq m, j \neq i} \frac{x-x_j}{x_m-x_j}
    % phi'(xi) = [ l_0'(xi)  l_1'(xi) ... l_{N+1}'(xi) ]
    for m=1:N+1
        for i=1:N+1
            if i==m
                continue
            end
            tmp = 1;
            for j=1:N+1
                if (j==i) || (j==m)
                    continue
                end
                tmp = tmp*(xi-xin(j))/(xin(m)-xin(j));
            end
            phi_xi(m) = phi_xi(m)+tmp/(xin(m)-xin(i));
        end
    end

end
