function [coordinates,newElements,varargout] ...
           = refineNVB5(coordinates,elements,varargin)

%refineNVB5: local refinement of finite element mesh by newest vertex
%           bisection, where marked elements are refined by five bisections
%
%Usage:
%
%[coordinates,elements,dirichlet,neumann] ...
%    = refineNVB5(coordinates,elements,dirichlet,neumann,marked)
%
%Comments:
%
%    refineNVB5 expects as input a finite element mesh described by the 
%    fields coordinates, elements, dirichlet and neumann. The vector 
%    marked contains the indices of elements which are refined by five bisections.
%    Further elements will be refined by newest vertex bisection to obtain
%    a regular triangulation.  Note that elements(j,3) provides the index
%    of the newest vertex of the j-th element.
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
%    S. Funken, D. Praetorius, P. Wissgott  10-07-08

nE = size(elements,1);
markedElements = zeros(nE,1); 
markedElements(varargin{end}) = 1; 
%*** Obtain geometric information on edges
[edge2nodes,element2edges,boundary2edges{1:nargin-3}] ...
    = provideGeometricData(elements,varargin{1:end-1});
%*** Mark edges for refinement
edge2newNode = zeros(max(max(element2edges)),1);
%markedElements
edge2newNode(element2edges(markedElements==1,:)) = 1;
swap = 1;
while ~isempty(swap)
    markedEdge = edge2newNode(element2edges);
    swap = find( ~markedEdge(:,1) & (markedEdge(:,2) | markedEdge(:,3)) );
    edge2newNode(element2edges(swap,1)) = 1;
end
%*** Generate new nodes (interior point of marked elements) 
markedElements(markedElements~=0) = size(coordinates,1) + (1:nnz(markedElements));
idx = find(markedElements);
coordinates(markedElements(idx),:) ...
    = (coordinates(elements(idx,1),:) + coordinates(elements(idx,2),:) + 2*coordinates(elements(idx,3),:))/4;
%*** Generate new nodes (midpoint of marked edges)
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
bisec123 = ( markedEdges(:,1) &  markedEdges(:,2) &  markedEdges(:,3) & ~markedElements);
bisec5   = ( markedEdges(:,1) &  markedEdges(:,2) &  markedEdges(:,3) & markedElements);
%*** Generate element numbering for refined mesh
idx = ones(nE,1);
idx(bisec1)   = 2; %*** bisec(1): newest vertex bisection of 1st edge
idx(bisec12)  = 3; %*** bisec(2): newest vertex bisection of 1st and 2nd edge
idx(bisec13)  = 3; %*** bisec(2): newest vertex bisection of 1st and 3rd edge
idx(bisec123) = 4; %*** bisec(3): newest vertex bisection of all edges
idx(bisec5)   = 6; %*** bisec(5): newest vertex bisection of all edges + interior node property
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
newElements([idx(bisec123),1+idx(bisec123),2+idx(bisec123),3+idx(bisec123)],:) ...
    = [newNodes(bisec123,1),elements(bisec123,3),newNodes(bisec123,3); ...
       elements(bisec123,1),newNodes(bisec123,1),newNodes(bisec123,3); ...
       newNodes(bisec123,1),elements(bisec123,2),newNodes(bisec123,2); ...
       elements(bisec123,3),newNodes(bisec123,1),newNodes(bisec123,2)];
newElements([idx(bisec5),1+idx(bisec5),2+idx(bisec5),3+idx(bisec5),4+idx(bisec5),5+idx(bisec5)],:) ...
    = [newNodes(bisec5,3),newNodes(bisec5,1),markedElements(bisec5); ...
       elements(bisec5,3),newNodes(bisec5,3),markedElements(bisec5); ...
       elements(bisec5,1),newNodes(bisec5,1),newNodes(bisec5,3); ...
       newNodes(bisec5,1),elements(bisec5,2),newNodes(bisec5,2); ...
       newNodes(bisec5,2),elements(bisec5,3),markedElements(bisec5); ...
       newNodes(bisec5,1),newNodes(bisec5,2),markedElements(bisec5)];