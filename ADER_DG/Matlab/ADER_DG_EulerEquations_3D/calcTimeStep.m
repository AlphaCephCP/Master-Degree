% Computes the maximum admissable time step according to the CFL
% condition for a generic nonliner PDE.
%
% Parameters:
% - curTime     : current time
% - tEnd        : end of simulation time
% - CFLmax      : maximum possible CFL number for the DG scheme
%
% Returns:
% - dt          : the time step
%
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/CalcTimeStep.f90
%

function dt=calcTimeStep(curTime,tEnd,CFLmax)
    % import global variables 
    global nDOF;  
    global nElem;
    global dx;      % mesh spacing along each dimension
    global nDim;
    global uh;      % the coefficients of the DG polynomial

    normalVectors = eye(nDim);  % 3x3

    dt = 1e20; 
    for iElem=1:nElem
        for k=1:nDOF(3)         % z
            for j=1:nDOF(2)     % y
                for i=1:nDOF(1) % x
                    denom = 0.0;
                    for iDim=1:nDim
                       % A method is stable only if the largest eigenvalue
                       % is < 1 for all reduced wavenumbers
                       Lambda = EigenvaluesEuler(uh(:,i,j,k,iElem)',normalVectors(:,iDim)');
                       lambdaMax = max(abs(Lambda));
                       denom = denom + lambdaMax/dx(iDim);
                    end
                    dt = min(dt, CFLmax/denom);
                end
            end
        end
    end
   
    if (curTime+dt > tEnd)
        dt = tEnd-curTime;
    end
end



