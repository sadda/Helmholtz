function etaH = computeEtaH(x,coordinates,elements,dirichlet,neumann,f,g)

%computeEtaH: computes hierarchical error estimator for finite element
%             solution of Laplace problem with mixed Dirichlet-Neumann
%             boundary conditions.
%
%Usage:
%
%etaH = computeEtaH(x,coordinates,elements,dirichlet,neumann,f,g)
%
%Comments:
%
%    The column vector x contains the nodal values of the P1 finite element
%    solution. The corresponding finite element mesh is given in terms of
%    coordinates, elements, dirichlet and neumann. The volume force <f.m> and
%    the Neumann data <g.m> are given as M-files. Either file is assumed
%    to take N evaluation points as (N x 2) matrix and to return an (N x 1) 
%    column vector.
%
%    The function returns the column vector etaH where etaH(J) is the
%    squared error indicator associated with the j-th element. These values
%    may be used to mark triangles for refinement. In particular, the 
%    value of the hierarchical error estimator is given by sqrt(sum(etaH).
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
[edge2nodes,element2edges,dirichlet2edges,neumann2edges] ...
    = provideGeometricData(elements,dirichlet,neumann);
%*** First vertex of elements and corresponding edge vectors 
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;
%*** Vector of element volumes 2*|T|
area2 = d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
%*** Elementwise gradient of uh
u21 = repmat( (x(elements(:,2))-x(elements(:,1)))./area2, 1,2);
u31 = repmat( (x(elements(:,3))-x(elements(:,1)))./area2, 1,2);
grad = (d31.*u21 - d21.*u31)*[0 -1 ; 1 0];
%*** Elementwise integrated gradients of hat functions --> 2*int(T,grad Vj)
grad1 = [d21(:,2)-d31(:,2) d31(:,1)-d21(:,1)];
grad2 = [d31(:,2) -d31(:,1)];
grad3 = [-d21(:,2) d21(:,1)];
%*** Compute volume contribution of rT (contribution of element bubble)
fsT = feval(f,c1+(d21+d31)/3);
rT = area2.*fsT/120;
%*** Compute volume contributions of rE edgewise (contribution of edge bubbles)
rE = repmat(area2.*fsT/24,1,3) - [sum((grad1+grad2).*grad,2) ...
                                  sum((grad2+grad3).*grad,2) ...
                                  sum((grad3+grad1).*grad,2)]/6;
rE = accumarray(element2edges(:),rE(:),[size(edge2nodes,1) 1]);
%*** Incorporate Neumann contributions to rE
if ~isempty(neumann)
    cn1 = coordinates(edge2nodes(neumann2edges,1),:);
    cn2 = coordinates(edge2nodes(neumann2edges,2),:);
    gmE = feval(g,(cn1+cn2)/2);
    rE(neumann2edges) = rE(neumann2edges) + sqrt(sum((cn2-cn1).^2,2)).*gmE/6;
end
%*** Incorporate Dirichlet data to rE
rE(dirichlet2edges) = 0;
%*** Compute error indicators
etaH = rT.^2 + accumarray(repmat(1:nE,1,3)',rE(element2edges(:)).^2,[nE 1]);
