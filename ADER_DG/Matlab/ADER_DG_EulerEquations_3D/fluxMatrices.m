% Computes the components which contribute to the flux matrix on the reference element.
%
% Parameters:
% - N       : degree of the approximation polynomial
% - xGPN    : quadrature points, typicall the Gauss-Legendre nodes
%
% Returns:
% - FLm     : left contribution to the left flux matrix,  (m = left  of the interface)
% - FLp     : right contribution to the left flux matrix  (p = right of the interface)
% - FRm     : left contribution to the right flux matrix  (m = left  of the interface)
% - FRp     : right contribution to the right flux matrix (p = right of the interface)
% - FLcoeff : extrapolated data onto the left boundary
% - FRcoeff : extrapolated data onto the right boundary
%
%
% Based on Michael Dumbser's code given in 
% ADER_DG_EulerEquations_3D/Init.f90. Note that the final 
% flux matrix is not assembled here.
%

function [FLm,FLp,FRm,FRp,FLcoeff,FRcoeff]=fluxMatrices(N,xGPN)
    phi0 = BaseFunc1D(N,xGPN,0); % interpolation polynomial 
                                 % running through the quadrature points
                                 % evaluated at 0.
    phi1 = BaseFunc1D(N,xGPN,1); % ... evaluated at 1.

    FLm = phi0'*phi1;
    FLp = phi0'*phi0;
    FRm = phi1'*phi1;
    FRp = phi1'*phi0;

    FLcoeff = phi0;
    FRcoeff = phi1;
end
