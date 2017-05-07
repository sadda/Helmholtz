function [x,coordinates,elements,indicators] ...
    = adaptiveAlgorithm(coordinates,elements,dirichlet,neumann,f,g,uD,nEmax,theta)

%adaptiveAlgorithm  adaptive finite element algorithm for two-dimensional 
%                   Laplace equation
%
%    adaptiveAlgorithm computes the P1-FEM solution of the Laplace 
%    equation 
%      - div(grad(u)) = f                   in Omega
%                   u = u_D                 on the Dirichlet boundary
%              d/dn u = g                   on the Neumann boundary
%    on a geometry described by triangles. 
%
%    adaptiveAlgorithm is the implementation of the following iterative
%    process for a given initial finite element mesh:
%
%          1) compute the discrete solution (via solveLaplace.m)
%          2) compute the elementwise error indicators (via computeEtaR.m)
%          3) mark the elements for refinement (via the Doerfler criterion)
%          4) refine the finite element mesh (via refineNVB.m)
%          5) return to 1)
%
%Usage:
%
%[x,coordinates,elements,indicators] ...
%    = adaptiveAlgorithm(coordinates,elements,dirichlet,neumann,f,g,uD,nemax,theta)
%
%Comments:
%
%    adaptiveAlgorithm expects as input an initial finite element mesh 
%    described by the fields coordinates, elements, dirichlet and neumann.
%
%    Volume force and boundary data are given as M-files <f.m>, <g.m>, and 
%    <uD.m>. Either of these M-files is assumed to take N evaluation
%    points as (N x 2) matrix and to return a (N x 1) column vector.
%
%    The stopping criterion is realized via the maximal number of elements 
%    NEMAX of the finite element mesh, and the computation is stopped a
%    soon as the current number of elements is larger than NEMAX.
%
%    The parameter THETA in (0,1) corresponds to the marking of elements by
%    use of the Doerfler marking with respect to the residual-based error
%    estimator. THETA = 1 means that (essentially) all elements are marked
%    for refinement, whereas small THETA leads to highly adapted meshes.
%
%    The function returns the adaptively generated mesh in terms of 
%    the matrices coordinates and elements, the vector x of the nodal
%    values of the P1-FEM solution as well the corresponding
%    refinement indicators, i.e., the value of the error estimator is
%    given by sqrt(sum(indicators)).
%
%    Please see <example1.m> for an example on the use of the function.
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

while 1
    %*** Compute discrete solution
    x = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,uD);
    %*** Compute refinement indicators
    indicators = computeEtaR(x,coordinates,elements,dirichlet,neumann,f,g);
    %*** Stopping criterion
    if size(elements,1) >= nEmax
        break
    end
    %*** Mark elements for refinement
    [indicators,idx] = sort(indicators,'descend');
    sumeta = cumsum(indicators);
    ell = find(sumeta>=sumeta(end)*theta,1);
    marked = idx(1:ell);
    %*** Refine mesh
    [coordinates,elements,dirichlet,neumann] = ...
       refineNVB(coordinates,elements,dirichlet,neumann,marked);
end
