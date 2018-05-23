function [p] = Newton_interp(x,y,z)
% Usage: p = Newton(x,y,z)
%
% This routine is based off of formula 8.16 from 
% Numerical Mathematics by Quarteroni.  It uses the 
% Newton Divided Difference Method to compute the 
% Newton Interpolating polynomial of a function at a point z recursively.  
% %
% Inputs:   x  nodal locations [array of length n] 
%           y  data values [array of length n]
%           z  evaluation point
% Outputs:  p  p(z)
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

% check inputs
if (length(x) ~= length(y)) 
   error('Newton error: (x,y) have different sizes');
end

p=0; %initialize p
F(:,1)=y'; %initialize divided difference matrix

n=max(size(y));

% Newton divided differences based off of program 66 
% from Numerical Mathematics by Quarteroni
for j=2:n
    for i=j:n
        F(i,j)=(F(i-1,j-1)-F(i,j-1))./(x(i-j+1)-x(i));
    end
end

% iteration to add k=1:n sums
for k=1:n
   % initialize w (the kth nodal polynomial)
   w = 1;
   % iterate over data to construct wk(z)
   for j=1:k
      % exclude the k-th data point
      if (j ~= k)
          w = w*(z-x(j));
      end
   end
   % construct the interpolating polynomial by equation 8.16
   p=p+F(k,k)*w;
   % F(k,k) because the relevant coefficient is on the diagonal
end






