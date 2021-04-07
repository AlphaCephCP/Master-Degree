% Computes the volume integral for one element.
%
% Parameters:
% - lqhi    : space-time degrees of freedom
% - lFhi_x  : x direction of nonlinear flux tensor in each space-time DOF
%             (nVar, nDOF(1),nDOF(2),nDOF(3))
% - lFhi_y  : y direction
%             (nVar,nDOF(1),nDOF(2),nDOF(3))  
% - lFhi_z  : z direction
%             (nVar,nDOF(1),nDOF(2),nDOF(3))  
%
% Returns:
% - lduh    : update DOF for one element
%
% Based in Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/VolumeIntegral.f90.
%
% Modifications:
% - the flux tensor has been split up into its three dimensions
%
% TODO: why not change the order of the DOFs in flux tensors?
%

function lduh = volumeIntegral(lqhi,lFhi_x,lFhi_y,lFhi_z)
    % import global variables
    global nVar;
    global nDOF;
    global wGPN;
    global Kxi;
    global dx;      % spacing in x,y,z direction

    % We always need the transpose of Kxi. Let's precompute is once.
    Kxi_T = transpose(Kxi);

    % all combinations of weights (~ 3D outer product)
    %weights = bsxfun(@times,permute(wGPN'*wGPN,[1 3 2]),wGPN);

    % Initialise the update DOF (spatial DOF)
    lduh = zeros(nVar,nDOF(1),nDOF(2),nDOF(3));

    %
    % (1) x-direction
    %
    weight = (wGPN'*wGPN)./dx(1);
    for k=1:nDOF(3)
        for j=1:nDOF(2)
            lduh(:,:,j,k) = lduh(:,:,j,k) + (lFhi_x(:,:,j,k)*Kxi_T) .* weight(j,k);
        end
    end

    % 
    % (2) y-direction
    %
    % Matlab can't handle matrices containing strides. Hence we change the data layout 
    % in order to have contiguous data in memory
    lduh_y = permute(lduh,[1 3 2 4]);

    weight = (wGPN'*wGPN)./dx(2);
    for k=1:nDOF(3)
        for i=1:nDOF(1)
            lduh_y(:,:,i,k) = lduh_y(:,:,i,k) + (reshape(lFhi_y(:,i,:,k),nVar,nDOF(2))*Kxi_T) .* weight(i,k);
            % ... (lFhi_y(:,:,i,k)*Kxi_T) .* weight(i,k); % after reordering? (see TODO) 
        end
    end
    lduh = permute(lduh_y,[1 3 2 4]); % undo permutation

    % 
    % (3) z-direction
    %
    lduh_z = permute(lduh,[1 4 2 3]);

    weight = (wGPN'*wGPN)./dx(3);
    for j=1:nDOF(2)
        for i=1:nDOF(1)
            lduh_z(:,:,i,j) = lduh_z(:,:,i,j) + (reshape(lFhi_z(:,i,j,:),nVar,nDOF(3))*Kxi_T) .* weight(i,j);
        end
    end
    lduh = permute(lduh_z,[1 3 4 2]);
end



