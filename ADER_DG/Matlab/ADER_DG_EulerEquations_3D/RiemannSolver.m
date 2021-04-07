% Solves the Riemann Problem for one face.
%
% Parameters:
% - iFace   : The face that is to be updated through solving the Riemann problem
%
% Returns:
% - iFace   : updated face struct
%
% Based on Michael Dumbser's code given in
% ADER_DG_EulerEquations_3D/RiemannSolver.f90
%

function iFace=RiemannSolver(iFace)
    % import global variables
    global nVar;  % number of PDE variables
    global wGPN;  % Gauss-Legendre weights 
    global nDOF;  % degrees of freedom
    global gamma; % specific heat ratio


    % Compute the average states from the left and the right, which 
    % we need to compute the numerical dissipation 
    qAvgL = zeros(nVar,1);
    qAvgR = zeros(nVar,1);


    % compute weighted average on both sides of the face (2d)
    % outer product gives all combinations of weights
    weights = wGPN'*wGPN;

    % consider the 3d matrix as nDOFz (cardinality) matrices of size (nDOFy x nVar)
    % Each row of the 2d matrix sized (nDOFy x nVar) is weighted; then 
    % the rows are summed up. This gives the average value in z direction
    % on one Gauss-Legendre point. Repeated for all nodes in z direction.
    for k=1:nDOF(3)
        for j=1:nDOF(2)
            qAvgL = qAvgL + weights(j,k)*iFace.qL(:,j,k);
            qAvgR = qAvgR + weights(j,k)*iFace.qR(:,j,k);       
        end
    end  
 

    % Here, we implement a very simple Rusanov scheme with scalar 
    % dissipation (smax*Id). We can change this into a more 
    % sophisticated Osher or HLLEM Riemann solver whenever needed!        
    eigValVecL = EigenvaluesEuler(qAvgL',iFace.nv);
    eigValVecR = EigenvaluesEuler(qAvgR',iFace.nv);
    smax = max(abs([eigValVecL eigValVecR])); % maximum wave speed


    % We now compute the numerical flux. Note that the scheme is at the moment written 
    % in CONSERVATION FORM => no fluctuations, but real fluxes. 
    % Later, this will be converted into the left and right fluctuations. 
    % readable version in Toro (2009): Eqn (10.55)
    iFace.FL = 0.5 * (iFace.FR + iFace.FL) - 0.5*smax*(iFace.qR - iFace.qL);
    iFace.FR = iFace.FL;

    lFbndL = iFace.FL;
    lFbndR = iFace.FR;
end
