function [coordinates,newElements,varargout] ...
           = refineMRGB(coordinates,elements,varargin)

%refineRGB: local refinement of finite element mesh by red-green-blue
%           refinement, where marked elements are red-refined.
%
%Usage:
%
%[coordinates,elements,dirichlet,neumann] ...
%    = refineMRGB(coordinates,elements,dirichlet,neumann,marked)
% 
%    refineMRGB expects as input a finite element mesh described by the 
%    fields coordinates, elements, dirichlet and neumann. The vector 
%    marked contains the indices of elements which are refined by
%    red-refinement, i.e. which are split into 4 similar triangles.
%    refineMRGB realizes an algorithm from [Carstensen, Constr. Approx.
%    20 (2004)] which guarantees H1-stability of the L2-projection,
%    where classical RGB-refinement and NVB-refinement are merged.
%    Note that elements(j,3) provides the index of the newest vertex of 
%    the j-th element, and the edge spanned by the other two vertices
%    is the refinement edge used for mesh-closure by NVB.
% 
%    The function returns the refined mesh in terms of the same data as
%    for the input.
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
%    S. Funken, D. Praetorius, P. Wissgott  13-10-11

% function [coordinates,newElements,varargout] ...
%            = refineMNVB(coordinates,elements,varargin)

% [COORDINATES,ELEMENTS,DIRICHLET,NEUMANN] ...
%   = refineMNVB(COORDINATES,ELEMENTS,DIRICHLET,NEUMANN,MARKED)
%
% If an element is marked, we mark all of its edges for refinement.
% Whereas usual newest vertex bisection would use bisec(3) for marked
% elements, we use red-refinement. The closure of the refined 
% triangulation is done by certain newest-vertex bisection
%
% MNVB = modified new vertex bisection, where bisec(3) is replaced
% by red-refinement if all edges of an element are marked for refinement.

markedElements = varargin{end};
nE = size(elements,1);
%*** Obtain geometric information on edges
[edge2nodes,element2edges,boundary2edges{1:nargin-3}] ...
    = provideGeometricData(elements,varargin{1:end-1});
%*** Mark edges for refinement
edge2newNode = zeros(max(max(element2edges)),1);
edge2newNode(element2edges(markedElements,:)) = 1;
swap = 1;
while ~isempty(swap)
    markedEdge = edge2newNode(element2edges);
    swap = find( ~markedEdge(:,1) & (markedEdge(:,2) | markedEdge(:,3)) );
    edge2newNode(element2edges(swap,1)) = 1;
end
%*** Generate new nodes
edge2newNode(edge2newNode~=0) = size(coordinates,1) + (1:nnz(edge2newNode));
idx = find(edge2newNode);
coordinates(edge2newNode(idx),:) ...
    = (coordinates(edge2nodes(idx,1),:)+coordinates(edge2nodes(idx,2),:))/2;
%*** Refine boundary conditions
for j = 1:nargout-2
    boundary = varargin{j};
    if ~isempty(boundary)
        newNodes = edge2newNode(boundary2edges{j});
        markedEdges = find(newNodes);
        if ~isempty(markedEdges)
            boundary = [boundary(~newNodes,:); ...
                        boundary(markedEdges,1),newNodes(markedEdges); ...
                        newNodes(markedEdges),boundary(markedEdges,2)];
        end
    end
    varargout{j} = boundary;
end
%*** Provide new nodes for refinement of elements
newNodes = edge2newNode(element2edges);
%*** Determine type of refinement for each element
markedEdges = (newNodes~=0);
none = ~markedEdges(:,1);
bisec1   = ( markedEdges(:,1) & ~markedEdges(:,2) & ~markedEdges(:,3) );
bisec12  = ( markedEdges(:,1) &  markedEdges(:,2) & ~markedEdges(:,3) );
bisec13  = ( markedEdges(:,1) & ~markedEdges(:,2) &  markedEdges(:,3) );
red      = ( markedEdges(:,1) &  markedEdges(:,2) &  markedEdges(:,3) );
%*** Generate element numbering for refined mesh
idx = ones(nE,1);
idx(bisec1)  = 2; %*** bisec(1): newest vertex bisection of 1st edge
idx(bisec12) = 3; %*** bisec(2): newest vertex bisection of 1st and 2nd edge
idx(bisec13) = 3; %*** bisec(2): newest vertex bisection of 1st and 3rd edge
idx(red)     = 4; %*** red: element is refined into 4 similar elements

idx = [1;1+cumsum(idx)];
%*** Generate new elements
newElements = zeros(idx(end)-1,3);
newElements(idx(none),:) = elements(none,:);
newElements([idx(bisec1),1+idx(bisec1)],:) ...
    = [elements(bisec1,3),elements(bisec1,1),newNodes(bisec1,1); ...
       elements(bisec1,2),elements(bisec1,3),newNodes(bisec1,1)];
newElements([idx(bisec12),1+idx(bisec12),2+idx(bisec12)],:) ...
    = [elements(bisec12,3),elements(bisec12,1),newNodes(bisec12,1); ...
       newNodes(bisec12,1),elements(bisec12,2),newNodes(bisec12,2); ...
       elements(bisec12,3),newNodes(bisec12,1),newNodes(bisec12,2)]; 
newElements([idx(bisec13),1+idx(bisec13),2+idx(bisec13)],:) ...
    = [newNodes(bisec13,1),elements(bisec13,3),newNodes(bisec13,3); ...
       elements(bisec13,1),newNodes(bisec13,1),newNodes(bisec13,3); ...
       elements(bisec13,2),elements(bisec13,3),newNodes(bisec13,1)];
%newElements([idx(bisec123),1+idx(bisec123),2+idx(bisec123),3+idx(bisec123)],:) ...
%    = [newNodes(bisec123,1),elements(bisec123,3),newNodes(bisec123,3); ...
%       elements(bisec123,1),newNodes(bisec123,1),newNodes(bisec123,3); ...
%       newNodes(bisec123,1),elements(bisec123,2),newNodes(bisec123,2); ...
%       elements(bisec123,3),newNodes(bisec123,1),newNodes(bisec123,2)];
newElements([idx(red),1+idx(red),2+idx(red),3+idx(red)],:) ...
    = [elements(red,1),newNodes(red,1),newNodes(red,3); ...
       newNodes(red,1),elements(red,2),newNodes(red,2); ...
       newNodes(red,3),newNodes(red,2),elements(red,3); ...
       newNodes(red,2),newNodes(red,3),newNodes(red,1)];
