% Computes the initial conditions on all nodes in conservative form. Here, we assume 
% a nodal basis. Otherwise, we would have to do formal L2 projection,  i.e. integration 
% of the initial condition and multiplication with the inverse mass matrix. 
%
% Returns:
% - u0 : initial data for all elements along all space DOFs
%
% Based on Michael Dumbser's code given in ADER_DG_EulerEquations_3D/Init.f90. 
%
% Note that the loops over the DOF have been moved.
%

function u0=initialField()
    % import global variables
    global xGPN;                            % Gauss-Legendre nodes used for our basis 
    global x;                               % physical position of the nodes in the grid
    global verticesOfElem;                  % connectivity from the element to the nodes
    global nDOF;                            % degrees of freedom
    global dx;                              % mesh spacing
    global nVar;                            % number of variables of the PDE system
    global nElem;                           % total number of elements

    % local variables
    sigma = [0.05; 0.05; 0.05];             % half-width
    VBase = [1.0,0.0,0.0,0.0,1.0];          % base state, primitive variables
    ampl  = [0.0,0.0,0.0,0.0,1.0e-3];       % pertubation amplitude vector


    for iElem=1:nElem
        x0 = x(:,verticesOfElem(1,iElem));  % physical coordinate of lower left vertex
                                            % in element
        % compute inital condition
        for k=1:nDOF(3)
            for j=1:nDOF(2)
                for i=1:nDOF(1)
                    xNode = x0 + ([xGPN(i), xGPN(j), xGPN(k)].*dx)';
                    V0 = VBase(:) + ampl(:).*exp(-0.5*sum(xNode.^2./sigma.^2));
                    u0(:,i,j,k,iElem) = Prim2ConsEuler(V0);
                end
            end
        end

    end % loop over elements
                         
end
