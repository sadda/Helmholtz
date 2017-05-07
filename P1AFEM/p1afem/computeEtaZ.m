function etaZ = computeEtaZ(x,coordinates,elements)

%computeEtaZ: computes ZZ-type error estimator for P1-finite element solution
%
%Usage:
%
%etaZ = computeEtaZ(x,coordinates,elements)
%
%Comments:
%
%    The column vector X contains the nodal values of a P1 finite element
%    function uh. The corresponding finite element mesh is given in terms of
%    coordinates and elements. The function computes the averaged gradient
%    of uh by use of the Clement operator Ah. The output provides the L2-norm
%    || grad(uh) - Ah(grad(uh)) ||_{L2(Omega)}
%
%    The function returns the column vector etaZ where etaZ(J) is the
%    squared error indicator associated with the j-th element, i.e.,
%    etaZ(J) = || grad(uh) - Ah(grad(uh)) ||_{L2(T_j)}^2.
%    These values may be used to mark triangles for refinement. In particular, the 
%    value of the ZZ-type error estimator is given by sqrt(sum(etaZ).
%    
%Remark:
%
%    This program is a supplement to the paper 
%    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
%    by S. Funken, D. Praetorius, and P. Wissgott. The reader should 
%    consult that paper for more information.   
%
%Authors:
% 
%    S. Funken, D. Praetorius, P. Wissgott  10-07-08

nE = size(elements,1);
nC = size(coordinates,1);
%*** First vertex of elements and corresponding edge vectors 
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;
%*** Vector of element volumes 2*|T|
area2 = d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
%*** Elementwise integrated gradient of uh --> 2*int(T,grad uh)
u21 = x(elements(:,2))-x(elements(:,1));
u31 = x(elements(:,3))-x(elements(:,1));
dudx = d31(:,2).*u21 - d21(:,2).*u31; 
dudy = d21(:,1).*u31 - d31(:,1).*u21;
%*** Compute coefficients for Clement interpolant Jh(grad uh)
zArea2 = accumarray(elements(:),[area2;area2;area2],[nC 1]);
qz = [accumarray(elements(:),[dudx;dudx;dudx],[nC 1])./zArea2, ...
      accumarray(elements(:),[dudy;dudy;dudy],[nC 1])./zArea2];
%*** Compute ZZ-refinement indicators
dudx = dudx./area2;
dudy = dudy./area2;
sz = [dudx dudx dudx dudy dudy dudy] - reshape(qz(elements,:), nE,6);
etaZ = (sum(sz.^2,2) + sum(sz.*sz(:,[2 3 1 5 6 4]),2)).*area2/12;
