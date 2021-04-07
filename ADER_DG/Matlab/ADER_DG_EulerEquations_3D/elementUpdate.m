% Updates the solution to the new time step.
%
% Parameters:
% - luh  : nonlinear flux tensor in each space-time DOF (l=element-local)
%          (nVar,nDOF(1),nDOF(2),nDOF(3))
% - lduh : spatial degrees of freedom                   (d=delta)
%          (nVar,nDOF(1),nDOF(2),nDOF(3))
% - dt   : (global) time step
%
% Returns:
% - luh  : the updated element-local coefficients of the DG polynomial
% 
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/CalcTimeStep.f90
%

function luh = elementUpdate(luh, lduh, dt)
    % import global variables
    global wGPN;
    global nVar;
    global nDOF;

    % project onto volume update lduh
    for k=1:nDOF(3)
        for j=1:nDOF(2)
            % no loop over i!
            constWeights = diag(wGPN(k)*wGPN(j));            
            weights = diag(wGPN)*constWeights;

            % multiply with the inverse of the mass matrix. For Gauss-Legendre
            % nodes, we simply need to divide by the Gaussian weights
            invWeights = inv(weights); 
            
            lduh(:,:,j,k) = lduh(:,:,j,k)*invWeights;
        end
    end

    % sum the contribution to the spatial DOFs
    luh = luh + dt*lduh;
end
