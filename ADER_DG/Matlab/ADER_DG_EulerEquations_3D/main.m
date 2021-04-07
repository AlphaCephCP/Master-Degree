% ADER-DG solver for the three-dimensional Euler Equations.
%
% Based on Michael Dumbser's code.
%

clear all;
close all;


addpath('PDE/');

format long;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Definition of variables,
%                                     initialisation,
%                                     precomputations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -----------------       This part can be modified by the user     ------------------- 

global N;               % polynomial degree of our approximation in space and time, PNPM
N               = 3;
CFL             = 0.9;  % Courant-Friedrichs-Lewy number < 1 
global nVar;            % number of variables of the PDE system
nVar            = 5;


% -----------------           Do NOT change the stuff below         ------------------- 
% The following variables contain important information about the numerical method. Do NOT change.  

global nDim;            % number of space dimensions. For readability we drop
nDim  = 3;
                        % the support for multiple dimensions in one single code.
                        % Note: We do no longer have a variable d as present 
                        % in the Fortran code. 

global xGPN;            % Gauss-Legendre nodes used for our basis 
global wGPN;            % Gauss-Legendre weights used for our basis
[xGPN,wGPN]=gaussLegendre(N+1,0,1); % precompute nodes and weights

nDOFx = N+1;            % number of degrees of freedom in x direction
nDOFy = N+1;            % number of degrees of freedom in y direction
nDOFz = N+1;            % number of degrees of freedom in z direction
nDOFt = N+1;            % number of degrees of freedom in time
global nDOF;            % number of degrees of freedom in space i=1:3 and time i=4
        nDOF  = [nDOFx, % Note: The nDOF array in the Fortran code uses i=0 for time 
                nDOFy,  
                nDOFz, 
                nDOFt];   

% PDE variables
global gamma;           % specific heat ratio, gamma = 1+R/c_v
gamma       = 1.4;
global nVar;            % number of variables of the PDE system
nVar        = 5;

% init everything related to timestepping
timestep    = 0;        % current time step [1 ... NMAX]
NMAX        = 100000;   % max. number of time steps
curTime     = 0.0;      % current time
tEnd        = 0.25;     % final time

% Maximum possible CFL numbers for the DG scheme, in terms of PNPM it is N==M.
PNPMTable   = [1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015];
CFL_DG      = @(N)PNPMTable(N+1); % extract the stability limit for a 
                                  % DG polynomial given its degree N
CFLmax      = CFL_DG(N)*CFL;      % stability limit




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Definition of the domain,
%                                      mesh setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -----------------       This part can be modified by the user     ------------------- 

xL    = [-0.5, -0.5, -0.5];                     % lower-left corner of the domain. 
xR    = [+0.5, +0.5, +0.5];                     % upper-right corner of the domain. 

IMAX = 2;                                       % number of elements in x,y,z direction
JMAX = 2;
KMAX = 2;

% -----------------           Do NOT change the stuff below         ------------------- 

nVertexPerElem = 2^nDim;                        % number of vertices per element
nEdgePerElem   = nDim*2^(nDim-1);               % number of edges per element
nFacePerElem   = 2*nDim;                        % number of faces per element

global dx;
global nElem;
global nNode;
global nFace;
global x;

dx      = (xR-xL)./[IMAX JMAX KMAX];     % regular mesh spacing
nElem   = IMAX*JMAX*KMAX;                % total number of elements
nNode   = (IMAX+1)*(JMAX+1)*(KMAX+1);    % total number of nodes
nFace   = nElem*3 + KMAX*IMAX ...        % total number of faces
                       + JMAX*IMAX ...
                       + JMAX*KMAX;             

x       = zeros(nDim,nNode);             % physical position of the nodes in the grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Matrices and operators
%                            (nomenclature as in Fortran code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MM,iMM]   = deal(zeros(N+1,N+1));              % Element mass matrix and its inverse
global Kxi;
Kxi = zeros(N+1,N+1);                    % Element stiffness matrix
dudx       = zeros(N+1,N+1);                    % Discrete derivative operator
[FLm, FLp] = deal(zeros(N+1,N+1));              % Left flux matrices, values attained on the left (m, -) 
                                                % and the right (p, +) side of the interface
[FRm, FRp] = deal(zeros(N+1,N+1));              % Right flux matrices, values attained on the left (m, -) 
                                                % and the right (p, +) side of the interface
[FLcoeff,...                                    % extrapolation coefficient to the left...
  FRcoeff] = deal(zeros(N+1));                  % ... and right boundary                                           
global F0;
global F1;
F0  = zeros(N+1);                        % time flux matrix (left)
F1  = zeros(N+1,N+1);                    % time flux matrix (right)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Main data structures of the 
%                                      ADER-DG scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global uh;
uh = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nElem);  % the coefficients of the DG polynomial
duh = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nElem); % the update coefficients of the DG polynomial
qBoundary = zeros(nVar,6,nDOF(2),nDOF(3),nElem); % time-averaged boundary-extrapolated data for the state vector Q in the element
FBoundary = zeros(nVar,6,nDOF(2),nDOF(3),nElem); % time-averaged boundary-extrapolated values for the normal flux F * n
% 6 faces in a 3D setup
qhi = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nElem);
Fhi_x = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nElem);
Fhi_y = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nElem);
Fhi_z = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nElem);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Set Connectivity Information
%                                        between 
%                                elements, nodes and faces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define connectivity between vertices settings for one element
% define node coordinates and the node numbers and connect elements and nodes
% idx2NodeNumber = zeros(IMAX+1,   % position of the nodes in terms of 
%                        JMAX+1,   % row, column and page indices (i,j,k)
%                        KMAX+1);  % in Fortran code called: idxn
idx2NodeNumber = zeros(IMAX+1, JMAX+1, KMAX+1);  % in Fortran code called: idxn

% idx2ElemNumber = zeros(IMAX,     % the lower left corner of each element
%                        JMAX,     % gives its position in terms of row, 
%                        KMAX);    % column and page indices (i,j,k)
idx2ElemNumber = zeros(IMAX,JMAX,KMAX);    % column and page indices (i,j,k)

% connectivity from the element to the nodes.
% An elements "owns" its vertices
global verticesOfElem;
verticesOfElem = zeros(nVertexPerElem,nElem);

% (1) define the connectivity between index coordinates and node numbers.
nodeNumber = 1;
for k=1:KMAX+1
    for j=1:JMAX+1
        for i=1:IMAX+1
            x(:,nodeNumber) = xL(:) + [i-1;j-1;k-1].*dx(:);    % physical coordinates
            idx2NodeNumber(i,j,k) = nodeNumber; % index coordinates
            nodeNumber = nodeNumber+1;
        end
    end
end


% (2) define the connectivity between the elements and the nodes.
elemNumber = 1;
for k=1:KMAX
    for j=1:JMAX
        for i=1:IMAX
            idx2ElemNumber(i,j,k) = elemNumber;
            verticesOfElem(:,elemNumber) = [ 
                idx2NodeNumber(i,j,k),           
                idx2NodeNumber(i+1,j,k),         
                idx2NodeNumber(i,j+1,k),
                idx2NodeNumber(i+1,j+1,k),
                idx2NodeNumber(i,j,k+1),
                idx2NodeNumber(i+1,j,k+1),
                idx2NodeNumber(i,j+1,k+1),
                idx2NodeNumber(i+1,j+1,k+1)      
            ];
            % now the element knows the IDs of its nodes

            elemNumber = elemNumber + 1;
        end
    end
end

% (3) define the connectivity between the faces and the elements.
%
% for the time being, follow the Fortran code and use an array of structs
% we could eliminate the branches easily but for readability we keep them.

% x faces
faceNumber = 0;
global Face; 
for k=1:KMAX
    for j=1:JMAX
        for i=1:IMAX+1 % !
            faceNumber = faceNumber+1;
            if i==1 
                Face(faceNumber).Left  = 0; % n/a
                Face(faceNumber).Right = idx2ElemNumber(i,j,k);
            elseif i==IMAX+1
                Face(faceNumber).Left  = idx2ElemNumber(i-1,j,k);
                Face(faceNumber).Right = 0; % n/a
            else % non-boundary
                Face(faceNumber).Left = idx2ElemNumber(i-1,j,k);
                Face(faceNumber).Right = idx2ElemNumber(i,j,k);
            end
            Face(faceNumber).nv = [1.0, 0.0, 0.0];  % set face normal vector
            Face(faceNumber).orientation = 1;       % for convenience
        end
    end
end
% y faces
for k=1:KMAX
    for j=1:JMAX+1 % !
        for i=1:IMAX
            faceNumber = faceNumber+1;
            if j==1
                Face(faceNumber).Left  = 0;  % n/a
                Face(faceNumber).Right = idx2ElemNumber(i,j,k);
            elseif j==JMAX+1
                Face(faceNumber).Left  = idx2ElemNumber(i,j-1,k);
                Face(faceNumber).Right = 0;  % n/a
            else
                Face(faceNumber).Left  = idx2ElemNumber(i,j-1,k);
                Face(faceNumber).Right = idx2ElemNumber(i,j,k);
            end
            Face(faceNumber).nv = [0.0, 1.0, 0.0];  % set face normal vector
            Face(faceNumber).orientation = 2;       % for convenience
        end
    end
end
% z faces
for k=1:KMAX+1 % !
    for j=1:JMAX
        for i=1:IMAX
            faceNumber = faceNumber+1;
            if k==1
                Face(faceNumber).Left  = 0;  % n/a
                Face(faceNumber).Right = idx2ElemNumber(i,j,k);
            elseif k==KMAX+1
                Face(faceNumber).Left  = idx2ElemNumber(i,j,k-1);
                Face(faceNumber).Right = 0;  % n/a
            else
                Face(faceNumber).Left  = idx2ElemNumber(i,j,k-1);
                Face(faceNumber).Right = idx2ElemNumber(i,j,k);
            end
            Face(faceNumber).nv = [0.0, 0.0, 1.0];  % set face normal vector
            Face(faceNumber).orientation = 3;       % for convenience
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Output-related 
%                                Settings 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the DOFs are on the non-equispaced Gauss-Legendre nodes; for 
% plotting, however, we want equispaced points including the 
% boundary values
SubOutputMatrix = subOutputMatrix(N,nDim,xGPN);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Compute some of the important 
%                                      matrices 
%                                in the ADER-DG method 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[MM,iMM] = massMatrix(N,xGPN,wGPN);
Kxi = stiffnessMatrix(N,xGPN,wGPN);
dudx = iMM * Kxi'; % discrete derivative operator 
[FLm,FLp,FRm,FRp,FLcoeff,FRcoeff] = fluxMatrices(N,xGPN);

% The time flux matrices for the ADER-DG predictor method are given by the 
% principle of upwinding in time (causality principle) 
F0 = FLcoeff; % upwinding in time = information comes from smaller times 
F1 = FRm;     % upwinding in time = information comes from smaller times 
K1 = F1 - Kxi;% stiffness matrix in reference element
global iK1;
global FCoeff;
iK1 = inv(K1); % matrix inverse
FCoeff = [FLcoeff' FRcoeff'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    Compute the
%                                  initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uh = initialField();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    Let's Start!
%                           Main loop of the ADER-DG scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time measurement
[tTmp,tPredictor,tVolumeIntegral,tRiemannSolver,tSurfaceIntegral] = deal(0.0);
tStart = tic;

curTime = 0.0;
for timestep=1:NMAX 
    if curTime >= tEnd
        break;
    end


    % (1) Compute the time step size according to the CFL condition 
    dt = calcTimeStep(curTime,tEnd,CFLmax);

    tTmp = toc(tStart);
    % (2) ADER predictor step 
    for iElem=1:nElem
        [qhi(:,:,:,:,iElem),       ...
         Fhi_x(:,:,:,:,iElem),     ...
         Fhi_y(:,:,:,:,iElem),     ...
         Fhi_z(:,:,:,:,iElem),     ...
         qBoundary(:,:,:,:,iElem), ...
         FBoundary(:,:,:,:,iElem)] = spaceTimePredictor(uh(:,:,:,:,iElem),dt);
    end
    tPredictor = tPredictor + (toc(tStart)-tTmp); % overall time spent on predictor
                                                  % across the entire simulation

    tTmp = toc(tStart);
    % (3) Compute the element volume integral
    for iElem=1:nElem
        duh(:,:,:,:,iElem) = volumeIntegral(qhi(:,:,:,:,iElem),   ...
                                            Fhi_x(:,:,:,:,iElem), ...
                                            Fhi_y(:,:,:,:,iElem), ...
                                            Fhi_z(:,:,:,:,iElem));
    end
    tVolumeIntegral = tVolumeIntegral + (toc(tStart)-tTmp);

    % (4) Set the boundary conditions (could depend on space and time) 
    boundaryConditions()


    % (5) Set the face data 
    %
    % If the face is on the boundary of our domain, copy the values
    % from the inside to the outside. 
    %
    % Check out the Fortran code for its definition of the faces. As we don't 
    % have pointer assignments at our disposal, we fix face-adjacent values manually.
    %
    shape = [nVar,nDOF(2),nDOF(3)]; % used to omit tensor dimensions of 1
    for iFace=1:nFace
        if Face(iFace).Left==0 % there is no left neighbour element
            switch Face(iFace).orientation
                case 1 % x faces
                    Face(iFace).qR = reshape(qBoundary(:,1,:,:,Face(iFace).Right),shape);
                    Face(iFace).FR = reshape(FBoundary(:,1,:,:,Face(iFace).Right),shape);
                    % qL and FL are set in boundaryConditions()
                case 2 % y faces
                    Face(iFace).qR = reshape(qBoundary(:,3,:,:,Face(iFace).Right),shape);
                    Face(iFace).FR = reshape(FBoundary(:,3,:,:,Face(iFace).Right),shape);
                case 3 % z face
                    Face(iFace).qR = reshape(qBoundary(:,5,:,:,Face(iFace).Right),shape);
                    Face(iFace).FR = reshape(FBoundary(:,5,:,:,Face(iFace).Right),shape);
                otherwise
                    disp('Face(#no).orientation not well-defined');
                    return;
            end % switch
        elseif Face(iFace).Right==0 % there is no right neighbour element
            switch Face(iFace).orientation
                case 1 % x faces
                    Face(iFace).qL = reshape(qBoundary(:,2,:,:,Face(iFace).Left),shape);
                    Face(iFace).FL = reshape(FBoundary(:,2,:,:,Face(iFace).Left),shape);
                case 2 % y faces
                    Face(iFace).qL = reshape(qBoundary(:,4,:,:,Face(iFace).Left),shape);
                    Face(iFace).FL = reshape(FBoundary(:,4,:,:,Face(iFace).Left),shape);
                case 3 % z faces
                    Face(iFace).qL = reshape(qBoundary(:,6,:,:,Face(iFace).Left),shape);
                    Face(iFace).FL = reshape(FBoundary(:,6,:,:,Face(iFace).Left),shape);
                otherwise
                    disp('Face(#no).orientation not well-defined');
                    return;
            end % switch
        else % we do have a left and a right neighbour element
            switch Face(iFace).orientation
                case 1 % x faces
                    Face(iFace).qR = reshape(qBoundary(:,1,:,:,Face(iFace).Right),shape);
                    Face(iFace).qL = reshape(qBoundary(:,2,:,:,Face(iFace).Left ),shape);
                    Face(iFace).FR = reshape(FBoundary(:,1,:,:,Face(iFace).Right),shape);
                    Face(iFace).FL = reshape(FBoundary(:,2,:,:,Face(iFace).Left ),shape);
                case 2 % y faces
                    Face(iFace).qR = reshape(qBoundary(:,3,:,:,Face(iFace).Right),shape);
                    Face(iFace).qL = reshape(qBoundary(:,4,:,:,Face(iFace).Left ),shape);
                    Face(iFace).FR = reshape(FBoundary(:,3,:,:,Face(iFace).Right),shape);
                    Face(iFace).FL = reshape(FBoundary(:,4,:,:,Face(iFace).Left ),shape);
                case 3 % z faces
                    Face(iFace).qR = reshape(qBoundary(:,5,:,:,Face(iFace).Right),shape);
                    Face(iFace).qL = reshape(qBoundary(:,6,:,:,Face(iFace).Left ),shape);
                    Face(iFace).FR = reshape(FBoundary(:,5,:,:,Face(iFace).Right),shape);
                    Face(iFace).FL = reshape(FBoundary(:,6,:,:,Face(iFace).Left ),shape);
                otherwise
                    disp('Face(#no).orientation not well-defined');
                    return;
            end % switch
        end % if
    end % for

    % (6) Solve the Riemann problems
    for iFace=1:nFace

        tTmp = toc(tStart);
        Face(iFace) = RiemannSolver(Face(iFace)); % update face data
        tRiemannSolver = tRiemannSolver + (toc(tStart) - tTmp);


        % update flux for the surface integral
        % If we had pointers, we could skip the data update below. The
        % Riemann solver should not know about the orientation of the 
        % faces; so we update our data structure manually.
        switch Face(iFace).orientation
            case 1 % x faces
                if Face(iFace).Left==0
                    FBoundary(:,1,:,:,Face(iFace).Right) = Face(iFace).FR;
                elseif Face(iFace).Right==0
                    FBoundary(:,2,:,:,Face(iFace).Left ) = Face(iFace).FL;
                else
                    FBoundary(:,1,:,:,Face(iFace).Right) = Face(iFace).FR;
                    FBoundary(:,2,:,:,Face(iFace).Left ) = Face(iFace).FL;
                end
            case 2 % y faces
                if Face(iFace).Left==0
                    FBoundary(:,3,:,:,Face(iFace).Right) = Face(iFace).FR;
                elseif Face(iFace).Right==0
                    FBoundary(:,4,:,:,Face(iFace).Left ) = Face(iFace).FL;
                else
                    FBoundary(:,3,:,:,Face(iFace).Right) = Face(iFace).FR;
                    FBoundary(:,4,:,:,Face(iFace).Left ) = Face(iFace).FL;
                end
            case 3 % z faces
                if Face(iFace).Left==0
                    FBoundary(:,5,:,:,Face(iFace).Right) = Face(iFace).FR;
                elseif Face(iFace).Right==0
                    FBoundary(:,6,:,:,Face(iFace).Left ) = Face(iFace).FL;
                else
                    FBoundary(:,5,:,:,Face(iFace).Right) = Face(iFace).FR;
                    FBoundary(:,6,:,:,Face(iFace).Left ) = Face(iFace).FL;
                end
        end
    end

    % (7) Compute the surface integrals of the test function multiplied 
    %     with the numerical flux 
    tTmp = toc(tStart);
    for iElem=1:nElem
        duh(:,:,:,:,iElem) = surfaceIntegral(duh(:,:,:,:,iElem),FBoundary(:,:,:,:,iElem));
    end
    tSurfaceIntegral = tSurfaceIntegral + (toc(tStart)-tTmp);

    % (8) Do the element update
    for iElem=1:nElem
        uh(:,:,:,:,iElem) = elementUpdate(uh(:,:,:,:,iElem),duh(:,:,:,:,iElem),dt);
    end

    % (9) Update the timestep
    curTime = curTime + dt;
    

end % loop over timesteps

disp('Total number of timesteps executed: ');
timestep
disp('Total runtime: ');
tElapsed = toc(tStart)
disp('Time spent on the kernels: ')
tPredictor
tVolumeIntegral
tRiemannSolver
tSurfaceIntegral

% export final solution (linearised data)
%fileID = fopen('uh.dat','w');
%fprintf(fileID,'%.12f \n',uh(:));
%fclose(fileID);

