% Sets the boundary conditions for the 3D Euler equations.
% 
% This code closely follows Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/BoundaryConditions.f90. The loops have 
% been replaced with matrix manipulations. The global face
% data structure is directly updated
%
% Remark: The iterations do not take the actual orientation of the face 
%         into account 
%

function boundaryConditions()
    % import global variables
    global Face;
    global nFace;
    global nDOF;
    global nVar;

    % fix boundary data
    primBC = [1, 0, 0, 0, 1];         % boundary conditions in primitive variables ...
    consBC = Prim2ConsEuler(primBC)'; % ... converted into conservation form

    % fluxes [fx, fy, fz]
    Fluxes = FluxEuler3D(consBC);

    % shape of the qL, qR, FL, FR matrices
    shape = [nVar, nDOF(2), nDOF(3)]; % note the remark above

    % set boundary conditions
    % 
    % For the moment, we use simple extrapolation (copy from inside the domain) 
    %
    % enumeration of all stuff (elements, nodes, faces) always starts at 1. When
    % we encounter a 0 (zero), we know this face is on the boundary
    for iFace=1:nFace
        if Face(iFace).Left==0
            fluxThroughFace = Fluxes * Face(iFace).nv';

            % boundary conditions are the same across the face
            % create linear data (repmat) and then interpret it as a
            % (nVar x nDOF(2) x nDOF(3)) matrix (reshape)
            Face(iFace).qL = reshape(repmat(consBC,          nDOF(2)*nDOF(3),1),shape); 
            Face(iFace).FL = reshape(repmat(fluxThroughFace, nDOF(2)*nDOF(3),1),shape);
        end % if
        if Face(iFace).Right==0
            fluxThroughFace = Fluxes * Face(iFace).nv';

            Face(iFace).qR = reshape(repmat(consBC,          nDOF(2)*nDOF(3),1),shape);
            Face(iFace).FR = reshape(repmat(fluxThroughFace, nDOF(2)*nDOF(3),1),shape);
        end % if
    end % for


end
