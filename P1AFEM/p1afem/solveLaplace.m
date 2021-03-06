function [x,energy] = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,uD)

%solveLaplace: computes P1-finite element solution for the two dimensional
%              Laplace equation with mixed Dirichlet-Neumann boundary 
%              condition
%
%    solveLaplace solves Laplace equation 
%      - div(grad(u)) = f                   in Omega
%                   u = u_D                 on the Dirichlet boundary
%              d/dn u = g                   on the Neumann boundary
%    on a geometry described by triangles. 
%
%Usage:
%
%[x,energy] = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,ud)
%
%Comments:
%
%    solveLaplace expects as input a finite element mesh described by the 
%    fields coordinates, elements, dirichlet and neumann. The volume
%    force f, the Neumann data g, and the (inhomogeneous) Dirichlet data
%    uD are given as M-files <f.m>, <g.m>, and <uD.m>. Either of these 
%    M-files is assumed to take n evaluation points as (n x 2) matrix and to
%    return an (n x 1) column vector. 
%
%    solveLaplace assembles the Galerkin data and solves the resulting 
%    linear system of equations to obtain the P1 finite element solution 
%    of the Laplace problem. The function returns a column vector X which
%    contains the nodal values of the FEM solution. Additionally, solveLaplace 
%    provides the energy of the discrete solution uh, i.e. 
%    energy = || grad(uh) ||_{L2(Omega)}^2.
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
x = zeros(nC,1);
%*** First vertex of elements and corresponding edge vectors 
c1 = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;
%*** Vector of element areas 4*|T|
area4 = 2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));
%*** Assembly of stiffness matrix
I = reshape(elements(:,[1 2 3 1 2 3 1 2 3])',9*nE,1);
J = reshape(elements(:,[1 1 1 2 2 2 3 3 3])',9*nE,1);
a = (sum(d21.*d31,2)./area4)';
b = (sum(d31.*d31,2)./area4)';
c = (sum(d21.*d21,2)./area4)';
A = [-2*a+b+c;a-b;a-c;a-b;b;-a;a-c;-a;c];
A = sparse(I,J,A(:));
%*** Prescribe values at Dirichlet nodes
dirichlet = unique(dirichlet);
x(dirichlet) = feval(uD,coordinates(dirichlet,:));
%*** Assembly of right-hand side
fsT = feval(f,c1+(d21+d31)/3);
b = accumarray(elements(:),repmat(12\area4.*fsT,3,1),[nC 1]) - A*x;
if ~isempty(neumann)
    cn1 = coordinates(neumann(:,1),:);
    cn2 = coordinates(neumann(:,2),:);
    gmE = feval(g,(cn1+cn2)/2);
    b = b + accumarray(neumann(:),...
                      repmat(2\sqrt(sum((cn2-cn1).^2,2)).*gmE,2,1),[nC 1]);
end
%*** Computation of P1-FEM approximation
freenodes = setdiff(1:nC, dirichlet);
x(freenodes) = A(freenodes,freenodes)\b(freenodes);
%*** Compute energy || grad(uh) ||^2 of discrete solution
energy = x'*A*x;
