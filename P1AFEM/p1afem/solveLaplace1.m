function [x,energy] = solveLaplace1(coordinates,elements,dirichlet,neumann,f,g,uD)

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
%[X,ENERGY] = SOLVELAPLACE1(COORDINATES,ELEMENTS,DIRICHLET,NEUMANN,F,G,UD)
%
%Comments:
%
%    solveLaplace expects as input a finite element mesh described by the 
%    fields COORDINATES, ELEMENTS, DIRICHLET, and NEUMANN. The volume
%    force F, the Neumann data G, and the (inhomogeneous) Dirichlet data
%    UD are given as M-files <f.m>, <g.m>, and <uD.m>. Either of these 
%    M-files is assumed to take n evaluation points as (n x 2) matrix and to
%    return an (n x 1) column vector. 
%
%    solveLaplace assembles the Galerkin data and solves the resulting 
%    linear system of equations to obtain the P1 finite element solution 
%    of the Laplace problem. The function returns a column vector X which
%    contains the nodal values of the FEM solution. Additionally, solveLaplace 
%    provides the energy of the discrete solution uh, i.e. 
%    ENERGY = || grad(uh) ||_{L2(Omega)}^2.
%
%Remark:
%
%    This program is a supplement to the paper 
%    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
%    by S. Funken, D. Praetorius, and P. Wissgott. This is the "medium" 
%    version of solveLaplace from Listing 2. The reader should consult the
%    paper for more information.   
%
%Authors:
% 
%    S. Funken, D. Praetorius, P. Wissgott  13-07-09

%*** Assembly of stiffness matrix in linear complexity
nE = size(elements,1);
nC = size(coordinates,1);
x = zeros(nC,1);
I = zeros(9*nE,1);
J = zeros(9*nE,1);
A = zeros(9*nE,1);
for i = 1:nE
   nodes = elements(i,:);
   B = [1 1 1 ; coordinates(nodes,:)'];
   grad = B \ [0 0 ; 1 0 ; 0 1];
   idx = 9*(i-1)+1:9*i;
   tmp = [1;1;1]*nodes;
   I(idx) = reshape(tmp',9,1);
   J(idx) = reshape(tmp,9,1);
   A(idx) = det(B)/2*reshape(grad*grad',9,1);
end
A = sparse(I,J,A,nC,nC);
%*** Prescribe values at Dirichlet nodes
dirichlet = unique(dirichlet);
x(dirichlet) = feval(uD,coordinates(dirichlet,:));
%*** Assembly of right-hand side
b = -A*x;
for i = 1:size(elements,1)
    nodes = elements(i,:);
    sT = [1 1 1]*coordinates(nodes,:)/3;
    b(nodes) = b(nodes) + det([1 1 1 ; coordinates(nodes,:)'])*feval(f,sT)/6;
end
for i = 1:size(neumann,1)
    nodes = neumann(i,:);
    mE = [1 1]*coordinates(nodes,:)/2;
    b(nodes) = b(nodes) + norm([1 -1]*coordinates(nodes,:))*feval(g,mE)/2;
end
%*** Computation of P1-FEM approximation
freenodes = setdiff(1:nC, dirichlet);
x(freenodes) = A(freenodes,freenodes)\b(freenodes);
%*** Compute energy || grad(uh) ||^2 of discrete solution
energy = x'*A*x;