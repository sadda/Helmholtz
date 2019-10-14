function [coordinates,newElements,varargout] ...
           = refineRGB(coordinates,elements,varargin)

%refineRGB: local refinement of finite element mesh by red-green-blue
%           refinement, where marked elements are red-refined.
%
%Usage:
%
%[coordinates,elements,dirichlet,neumann] ...
%    = refineRGB(coordinates,elements,dirichlet,neumann,marked)
% 
%    refineRGB expects as input a finite element mesh described by the 
%    fields coordinates, elements, dirichlet and neumann. The vector 
%    marked contains the indices of elements which are refined by
%    red-refinement, i.e. which are split into 4 similar triangles.
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

markedElements = varargin{end};
nE = size(elements,1);
%*** Sort elements such that first edge is longest
dx = coordinates(elements(:,[2,3,1]),1)-coordinates(elements,1);
dy = coordinates(elements(:,[2,3,1]),2)-coordinates(elements,2);
[hT,idxMax] = max(reshape(dx.^2+dy.^2,nE,3),[],2);
idx = ( idxMax==2 );
elements(idx,:) = elements(idx,[2,3,1]); 
idx = ( idxMax==3 );
elements(idx,:) = elements(idx,[3,1,2]); 
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
edge2newNode(find(edge2newNode)) = size(coordinates,1) + [1:nnz(edge2newNode)];
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
bisec1  = ( markedEdges(:,1) & ~markedEdges(:,2) & ~markedEdges(:,3) );
bisec12 = ( markedEdges(:,1) &  markedEdges(:,2) & ~markedEdges(:,3) );
bisec13 = ( markedEdges(:,1) & ~markedEdges(:,2) &  markedEdges(:,3) );
red     = ( markedEdges(:,1) &  markedEdges(:,2) &  markedEdges(:,3) );
%*** Generate element numbering for refined mesh
idx = ones(nE,1);
idx(bisec1)  = 2; %*** green = newest vertex bisection of 1st edge
idx(bisec12) = 3; %*** blue (right) = newest vertex bisection of 1st and 2nd edge
idx(bisec13) = 3; %*** blue (left) = newest vertex bisection of 1st and 3rd edge
idx(red)     = 4; %*** red refinement
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
newElements([idx(red),1+idx(red),2+idx(red),3+idx(red)],:) ...
    = [elements(red,1),newNodes(red,1),newNodes(red,3); ...
       newNodes(red,1),elements(red,2),newNodes(red,2); ...
       newNodes(red,3),newNodes(red,2),elements(red,3); ...
       newNodes(red,2),newNodes(red,3),newNodes(red,1)];
