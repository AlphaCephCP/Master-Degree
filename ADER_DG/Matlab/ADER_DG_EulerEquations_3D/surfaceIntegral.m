% Computes the surface integral for one element.
%
% Parameters:
% - lFBoundary  : time-averaged nonlinear flux tensor in each space-time DOF
% - lduh        : current status of the update DOFs as emerged from volume integral
% 
% Returns:
% - lduh        : update DOF for one element
%
% Based in Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/SurfaceIntegral.f90.
%

function lduh = surfaceIntegral(lduh,lFBoundary)
    % import global variables
    global nVar;
    global nDOF;
    global wGPN;
    global dx;      % spacing in x,y,z direction
    global FCoeff;  % [FLcoeff, FRcoeff], extrapolation coefficients to the left
                    % and right boundary

    FLcoeff = FCoeff(:,1)';
    FRcoeff = FCoeff(:,2)';

    % compute weighted average on both sides of the face (2d)
    % outer product gives all combinations of weights
    weights = wGPN'*wGPN;

    % Now multiply the numerical fluxes on the surfaces with the test functions 
    % and compute the surface integrals 
    
    % (1) x-faces
    for k=1:nDOF(3)
        for j=1:nDOF(2)
            for iVar=1:nVar
                lduh(iVar,:,j,k) = lduh(iVar,:,j,k) - weights(j,k)/dx(1)*(lFBoundary(iVar,2,j,k)*FRcoeff - lFBoundary(iVar,1,j,k)*FLcoeff);
            end
        end
    end     

    % (2) y-faces
    for k=1:nDOF(3)
        for i=1:nDOF(1)
            for iVar=1:nVar
                lduh(iVar,i,:,k) = lduh(iVar,i,:,k) - reshape(weights(i,k)/dx(2)*(lFBoundary(iVar,4,i,k)*FRcoeff - lFBoundary(iVar,3,i,k)*FLcoeff),1,1,4);
            end
        end
    end

    % (3) z-faces
    for j=1:nDOF(2)
        for i=1:nDOF(1)
            for iVar=1:nVar
                lduh(iVar,i,j,:) = lduh(iVar,i,j,:) - reshape(weights(i,j)/dx(3)*(lFBoundary(iVar,6,i,j)*FRcoeff - lFBoundary(iVar,5,i,j)*FLcoeff),1,1,1,4);
            end
        end
    end
end







