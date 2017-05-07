%example1 computes and visualizes the P1-FEM solution for some Laplace
%equation via an adaptive algorithm
%
%    example1 solves Laplace equation 
%      - div(grad(u)) = 1                 in Omega
%                   u = 0                 on the Dirichlet boundary
%              d/dn u = 0/1               on the Neumann boundary
%    on a geometry described by triangles. 
% 
%    First, the finite element mesh is loaded. Therefore, the program
%    expects the files <coordinates.dat>, <elements.dat>, <dirichlet.dat>,
%    and <neumann.dat> to be located at the files' directory.
%
%    Second, the initial mesh is uniformly refined by 4 successive calls of
%    refineRGB.m to provide a uniform finite element mesh with N=3072 elements.
%
%    Third, we set the parameters and call adaptiveAlgorithm.m, which 
%    returns the finite element solution for the problem as well as the
%    final mesh.
%
%    Finally, we visualize the discrete solution via the Matlab function
%    trisurf.m.
%
%Remark:
%
%    This program is a supplement to the paper 
%    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
%    by S. Funken, D. Praetorius, and P. Wissgott. The reader should 
%    consult that paper for more information. In Section 6.1 of the
%    paper, this example is discussed in more detail.
%
%Authors:
% 
%    S. Funken, D. Praetorius, P. Wissgott  10-07-08

warning off
addpath(strcat('..',filesep))
warning on

%*** Parameters
nEmax = 2e4;
theta = 0.25;

%*** Initialization
load elements.dat
load coordinates.dat
load dirichlet.dat
load neumann.dat

%*** Initial refinement yields a uniform mesh with N=3072 elements
for i=1:4
    [coordinates,elements,dirichlet,neumann] ...
        = refineRGB(coordinates,elements,dirichlet,neumann,1:size(elements,1));
end

%*** Call adaptive algorithm
[x,coordinates,elements,indicators] ...
    = adaptiveAlgorithm(coordinates,elements,dirichlet,neumann,@f,@g,@uD,nEmax,theta);

%*** Visualization
trisurf(elements(:,1:3),coordinates(:,1),coordinates(:,2),x','facecolor','interp')
view(15,22)
