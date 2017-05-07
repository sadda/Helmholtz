function out = g(x)
%Neumann boundary data: d/dn u = g

out = zeros(size(x,1),1);
out(x(:,1)==-1 & x(:,2)>=0) = 1;
out(x(:,2)==1 & x(:,1)>=0) = 1;
