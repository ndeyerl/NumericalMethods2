function [p] = Hermite_interp(x,y,dy,z)
% Usage: p = Hermite(x,y,dy,z)
%
% This routine is based off of example 8.6 from
% Numerical Mathematics by Quarteroni, since in this case we only have 1
% derivative of f (so mi=1).  It evaluates the Hermite interpolating
% polynomial of a function if the first derivative of the function is
% supplied (also known as the osculatory interpolation) at an evaluation point z.
% %
% Inputs:   x  nodal locations [array of length nu] 
%           y  data values [array of length nu]
%           dy data derivative values [array of length nu]
%           z  evaluation point
% Outputs:  p  p(z)
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

% check inputs
if (length(x) ~= length(y)) 
   error('Hermite error: (x,y) have different sizes');
end

if (length(y) ~= length(dy)) 
   error('Hermite error: (y,dy) have different sizes');
end

n = max(size(x));

%initialize polynomial p
p=0;

% iteration to add i=1:n sums
for i=1:n
    
    % lagrange basis function and derivative of lagrange basis function,
    % used in Ai and Bi
    l = 1;
    lprime = 0;
    % iterate over data to construct l(z)
    for j=1:n
        % exclude the i-th data point
        if (j ~= i)
            l = l.*(z-x(j))./(x(i)-x(j)); % based off Reynolds' lagrange.m
            lprime=lprime+1./(x(i)-x(j)); % lprime from page 350
        end
    end

    %compute A(i) and B(i)
    Ai=(1-((2.*(z-x(i))).*lprime)).*(l.^2);
    Bi=(z-x(i)).*(l.^2);
    
    % construct the interpolating polynomial by example 8.6
    p=p+(y(i).*Ai)+(dy(i).*Bi);

end