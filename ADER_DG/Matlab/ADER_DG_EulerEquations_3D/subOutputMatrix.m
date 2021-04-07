% Computes a mapping of the non-equispaced Gauss-Legendre
% nodes onto equispaced nodes.
%
% Parameters:
% - N    : degree of the interpolation polynomial
% - xGPN : interpolation points (typically Gauss-Legendre nodes)
% - nDim : number of space dimensions we simulate
% 
% Returns:
% - subOutputMatrix : matrix used for plotting on an equispaced subgrid
% 
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/Init.f90
%
function subOutputMatrix = subOutputMatrix(N,nDim,xGPN)
    % import global variables
    global nDOF;

    subOutputMatrix = zeros((N+1)^nDim,(N+1)^nDim);

    phi_i = zeros(N+1);
    phi_j = zeros(N+1);
    phi_k = zeros(N+1);

    subxi = linspace(0,1,N+1); % equispaced points on the subgrid
    % evaluate the basis functions at our plotting points
    cnt = 0;
    for k=1:N+1
        for j=1:N+1
            for i=1:N+1
                cnt = cnt+1;
                phi_i = BaseFunc1D(N,xGPN,subxi(i));
                phi_j = BaseFunc1D(N,xGPN,subxi(j));
                phi_k = BaseFunc1D(N,xGPN,subxi(k));
                count = 0;
                for kk=1:nDOF(3)
                    for jj=1:nDOF(2)
                        for ii=1:nDOF(1)
                            count = count+1;
                            aux = [phi_i(ii), phi_j(jj), phi_k(kk)];
                            val = prod(aux(1:nDim));
                            subOutputMatrix(count,cnt) = val;
                        end
                    end
                end
                %
            end
        end
    end
end
