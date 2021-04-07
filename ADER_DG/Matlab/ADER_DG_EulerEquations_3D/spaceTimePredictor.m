% Computes the space-time predictor.
%
% \hat{q}^{n+1} = inv(K1) * ( F0 * \hat{u} - Kxi \hat{F}^n )
%
% Parameters
% - luh : spatial degrees of freedom (nVar,nDOF(1),nDOF(2),nDOF(3))
% - dt  : timestep
%
% Returns:
% - lqhi    : time-averaged space-time degrees of freedom
% - lFhi_x  : time-averaged nonlinear flux tensor in each x direction
% - lFhi_y  : time-averaged nonlinear flux tensor in each y direction
% - lFhi_z  : time-averaged nonlinear flux tensor in each z direction
% - lqbnd   : time-averaged space-time degrees of freedom
% - lFbnd   : time-averaged nonlinear flux tensor in each space-time DOF
%
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/SpaceTimePredictor.f90. 
%
% Modifications:
% - all possible combinations of weights are precomputed (we can do this because
%   the code is 3D only)
% - the flux tensor lFh is split into its 3 components, which are stored in separate
%   tensors lFh_x, lFh_y, lFh_z for each direction
% - transpose(iK1) is computed beforehand (not in the loop)
% - joined computations in loop constructs have been separated as far as possible
% - general deloopification
%
% Note: Some trial code snippets demand for further fiddling around. Ignore these
% for the time being :-)
%


function [lqhi,lFhi_x,lFhi_y,lFhi_z,lQbnd,lFbnd] = spaceTimePredictor(luh,dt)
    % import global variables
    global nVar;
    global N;
    global nDOF; % 1..3 -> x,y,z; 4 -> time
    global nDim;
    global dx;
    global F0;
    global F1;
    global iK1;
    global Kxi;
    global wGPN;
    global FCoeff; % [FLcoeff, FRcoeff], extrapolation coefficients to the left
                   % and right boundary
    FLCoeff = FCoeff(:,1);
    FRCoeff = FCoeff(:,2);

    % output variables
    lqhi = zeros(nVar,nDOF(1),nDOF(2),nDOF(3));
    %lFhi = zeros(nVar,nDim,nDOF(1),nDOF(2),nDOF(3));
    lFhi_x = zeros(nVar,nDOF(1),nDOF(2),nDOF(3));
    lFhi_y = zeros(nVar,nDOF(1),nDOF(2),nDOF(3)); 
    lFhi_z = zeros(nVar,nDOF(1),nDOF(2),nDOF(3)); 

    lQbnd = zeros(nVar,6,nDOF(2),nDOF(3)); % 6 = 6 faces
    lFbnd = zeros(nVar,6,nDOF(2),nDOF(3));

    % local variables
    lqh = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(4)); % space-time degrees of freedom
    rhs0 = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(4));
    rhs = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(4));
    %lFh = zeros(nVar,nDim,nDOF(1),nDOF(2),nDOF(3),nDOF(4)); % nonlinear flux tensor in each space time DOF -> replaced with the three vectors
    lFh_x = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(4));
    lFh_y = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(4));
    lFh_z = zeros(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(4));

    tol = 1e-7; % tolerance


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preliminary stuff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % we always need the transpose of inv(K1)
    iK1_T = transpose(iK1);

    % W_ijk = w_i * w_j * w_k, generalised outer product
    W = bsxfun(@times,permute(wGPN'*wGPN,[1 3 2]),wGPN);
    

    % Trivial initial guess. Take over solution from previous timestep
    % (later on we will use time-extrapolated polynomials)
    for i=1:nDOF(4)
      lqh(:,:,:,:,i) = luh;
    end
    
    % Compute rhs0 = F0 * \hat{u}
    for k=1:nDOF(3)
        for j=1:nDOF(2)
            for i=1:nDOF(1)
                rhs0(:,i,j,k,:) = (W(i,j,k)*luh(:,i,j,k))*F0;
            end
        end
    end
     

    % for the C++ version
    % rhs = rhs0;
    
    
    % Say, W(:)' = [w1 w2 w3] and vecLength = 2. Then
    % weightVector = [w1 w1 w2 w2 w3 w3] implements a "replication by sections"
    % total length of weightVector: (nVar * (N+1))  * (N+1)^3
    % it's massive overhead, but we need the very same vector several times 
    vecLength = nVar * (N+1);
    weightVector = cell2mat(arrayfun(@(x) ones(1,vecLength)*x, W(:)', 'UniformOutput',false));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Let's start
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Discrete Picard iterations. This set of nested loops should (theoretically) 
    % be a dream for vectorization, since they are rather independent... 
    for iter=1:N+1
        for l=1:nDOF(4) % time DOF
            % (1) Compute the fluxes (once these fluxes are available, the 
            % subsequent operations are independent from each other) 
            for k=1:nDOF(3)
                for j=1:nDOF(2)
                    for i=1:nDOF(1)
                        %lFh(:,:,i,j,k,l) = FluxEuler3D(lqh(:,i,j,k,l));
                        flux = FluxEuler3D(lqh(:,i,j,k,l)); % 5 x 3
                        lFh_x(:,i,j,k,l) = flux(:,1);
                        lFh_y(:,i,j,k,l) = flux(:,2);
                        lFh_z(:,i,j,k,l) = flux(:,3);
                    end
                end
            end          

            % (2) Compute the "derivatives" (contributions of the stiffness matrix)
            % (2a) x direction
            for k=1:nDOF(3)
                for j=1:nDOF(2)
                    rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) - W(l,j,k)*dt/dx(1)* lFh_x(:,:,j,k,l)*Kxi;
                end
            end

            % (2b) y direction
            % Matlab needs the permutation stuff (can't handle matrices containing strides)
            matrix_y = permute(lFh_y, [1 3 2 4 5]);
            rhs_y = permute(rhs,[1 3 2 4 5]);
            for k=1:nDOF(3)
                for i=1:nDOF(1)
                    %rhs(:,i,:,k,l) = rhs(:,i,:,k,l) - W(l,i,k)*dt*dx(2)*lFh_y(:,i,:,k,l)*Kxi;
                    rhs_y(:,:,i,k,l) = rhs_y(:,:,i,k,l) - W(l,i,k)*dt/dx(2)*matrix_y(:,:,i,k,l)*Kxi;
                end
            end
            rhs = permute(rhs_y,[1 3 2 4 5]);

            % (2c) z direction
            matrix_z = permute(lFh_z,[1 4 2 3 5]);
            rhs_z = permute(rhs,[1 4 2 3 5]);
            for j=1:nDOF(2)
                for i=1:nDOF(1)
            %        rhs(:,i,j,:,l) = rhs(:,i,j,:,l) - W(l,i,j)*dt/dx(3)*lFh_z(:,i,j,:,l)*Kxi;
                    rhs_z(:,:,i,j,l) = rhs_z(:,:,i,j,l) - W(l,i,j)*dt/dx(3)*matrix_z(:,:,i,j,l)*Kxi;

                end
            end
            rhs = permute(rhs_z,[1 3 4 2 5]);

        end % time DOF
        
        % (3) Multiply with (K1)^(-1) to get the discrete time integral of the 
        % discrete Picard iteration 
        rhs_t = permute(rhs,[1 5 2 3 4]); % let matrix operation range over time DOFs
        lqh_t = permute(lqh,[1 5 2 3 4]); % Q - time - x - y - z
        for k=1:nDOF(3)
            for j=1:nDOF(2)
                for i=1:nDOF(1)     
                    % lqh(:,i,j,k,:) = 1./W(i,j,k) * rhs(:,i,j,k,:) * transpose(iK1);
                    lqh_t(:,:,i,j,k) = 1./W(i,j,k) * (rhs_t(:,:,i,j,k) * iK1_T);
                end
            end
        end

        % TODO: Check if we can make the permutation of lqh permanent. Below the
        % matrix multiplications seem to always range over the time DOFs
        lqh = permute(lqh_t,[1 3 4 5 2]); % Q - x - y - z - time


        % No termination criterion. We always process all Picard iterations.
       
    end % end of Picard iteration


    % (4) Immediately compute the time-averaged space-time polynomials 
    
    % for the time being, redo the permutation and interpret as multiple matrices
    % that have been glued together and exec the matrix-vector multiplication
    % 
    % Column major format prohibits writing 
    % lqhi = lqh_t * wGPN';
    %
    % C++ gives us the conversion for free.
    % 
    % +------------------------------------+   +-----------+
    % |         lqh_t(:,:,1,1,1)           |   |  wGPN(1)  |
    % +------------------------------------+   +-----------+
    % |         lqh_t(:,:,2,1,1)           | * |    ...    |
    % +------------------------------------+   +-----------+
    % |                ...                 |   | wGPN(N+1) |
    % +------------------------------------+   +-----------+
    % | lqh_t(:,:,nDOF(1),nDOF(2),nDOF(3)) |
    % +------------------------------------+
    %
    % 

    lqh_t   = permute(lqh,  [1 5 2 3 4]);
    lFh_x_t = permute(lFh_x,[1 5 2 3 4]);
    lFh_y_t = permute(lFh_y,[1 5 2 3 4]);
    lFh_z_t = permute(lFh_z,[1 5 2 3 4]);
    for k=1:nDOF(3)
        for j=1:nDOF(2)
            for i=1:nDOF(1)
                % (nVar x 1) = (nVar x nDOF(4)) * (nDOF(4) x 1)
                lqhi(:,i,j,k) = reshape(lqh_t(:,:,i,j,k),nVar,nDOF(4))*wGPN';

                % for each direction: multiply lFh_<direction> with weights
                % (nVar x 1) = (nVar x nDOF(4)) * (nDOF(4) x 1)
                lFhi_x(:,i,j,k) = reshape(lFh_x_t(:,:,i,j,k),nVar,nDOF(4))*wGPN';
                lFhi_y(:,i,j,k) = reshape(lFh_y_t(:,:,i,j,k),nVar,nDOF(4))*wGPN';
                lFhi_z(:,i,j,k) = reshape(lFh_z_t(:,:,i,j,k),nVar,nDOF(4))*wGPN';
            end
        end
    end

    % (5) Compute the bounday-extrapolated values for Q and F*n
    %
    % (5a) x-direction: face 1 and face 2
    %
    for k=1:nDOF(3)
        for j=1:nDOF(2)
            lQbnd(:,1,j,k) = lqhi(:,:,j,k)*FLCoeff;
            lQbnd(:,2,j,k) = lqhi(:,:,j,k)*FRCoeff;
            lFbnd(:,1,j,k) = lFhi_x(:,:,j,k)*FLCoeff;
            lFbnd(:,2,j,k) = lFhi_x(:,:,j,k)*FRCoeff;
        end
    end
    %
    % (5b) y-direction: face 3 and face 4
    %
    for k=1:nDOF(3)
        for i=1:nDOF(1)
            lQbnd(:,3:4,i,k) = reshape(lqhi(:,i,:,k),nVar,nDOF(2))*FCoeff;
            lFbnd(:,3:4,i,k) = reshape(lFhi_y(:,i,:,k),nVar,nDOF(2))*FCoeff;
        end
    end
    %
    % (5c) z-direction: face 5 and face 6
    %
    for j=1:nDOF(2)
        for i=1:nDOF(1)
            lQbnd(:,5:6,i,j) = reshape(lqhi(:,i,j,:),nVar,nDOF(3))*FCoeff;
            lFbnd(:,5:6,i,j) = reshape(lFhi_z(:,i,j,:),nVar,nDOF(3))*FCoeff;
        end
    end

end
