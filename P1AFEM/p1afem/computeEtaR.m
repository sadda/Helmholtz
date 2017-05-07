function etaR = computeEtaR(x,coordinates,elements,dirichlet,neumann,f,g)

%computeEtaR: computes residual-based error estimator for finite element
%             solution of Laplace problem with mixed Dirichlet-Neumann
%             boundary conditions.
%
%Usage:
%
%etaR = computeEtaR(x,coordinates,elements,dirichlet,neumann,f,g)
%
%Comments:
%
%    The column vector X contains the nodal values of the P1 finite element
%    solution. The corresponding finite element mesh is given in terms of
%    coordinates, elements, dirichlet and neumann. The volume force <f.m> and
%    the Neumann data <g.m> are given as M-files. Either file is assumed
%    to take N evaluation points as (N x 2) matrix and to return an (N x 1) 
%    column vector.
%
%    The function returns the column vector etaR where etaR(J) is the
%    squared error indicator associated with the j-th element. These values
%    may be used to mark triangles for refinement. In particular, the 
%    value of the residual error estimator is given by sqrt(sum(etaR).
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

[edge2nodes,element2edges,dirichlet2edges,neumann2edges] ...
     = provideGeometricData(elements,dirichlet,neumann);
%*** First vertex of elements and corresponding edge vectors
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;
%*** Vector of element volumes 2*|T|
area2 = d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
%*** Compute curl(uh) = (-duh/dy, duh/dx)
u21 = repmat(x(elements(:,2))-x(elements(:,1)), 1,2);
u31 = repmat(x(elements(:,3))-x(elements(:,1)), 1,2);
curl = (d31.*u21 - d21.*u31)./repmat(area2,1,2);
%*** Compute edge terms hE*(duh/dn) for uh
dudn21 = sum(d21.*curl,2);
dudn13 = -sum(d31.*curl,2);
dudn32 = -(dudn13+dudn21);
etaR = accumarray(element2edges(:),[dudn21;dudn32;dudn13],[size(edge2nodes,1) 1]);
%*** Incorporate Neumann data
if ~isempty(neumann)
  cn1 = coordinates(neumann(:,1),:);
  cn2 = coordinates(neumann(:,2),:);
  gmE = feval(g,(cn1+cn2)/2);
  etaR(neumann2edges) = etaR(neumann2edges) - sqrt(sum((cn2-cn1).^2,2)).*gmE;
end
%*** Incorporate Dirichlet data
etaR(dirichlet2edges) = 0;
%*** Assemble edge contributions of indicators
etaR = sum(etaR(element2edges).^2,2);
%*** Add volume residual to indicators
fsT = feval(f,(c1+(d21+d31)/3));
etaR = etaR + (0.5*area2.*fsT).^2;
