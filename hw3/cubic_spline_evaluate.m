function s = cubic_spline_evaluate(x, f, M, z)
% Usage: s = cubic_spline_evaluate(x, f, M, z)
%
% This routine evaluates the cubic spline of f at node x at an evaluation 
% point z by determining the nodal interval z is in, calling the appr-
% -opriate coefficients Mk which should have already been computed
% at the nodes and function values (x,f) using cubic_spline_coefficients
% then evaluating s at the evaluation point z. It is based off of 
% Quarteroni's Numerical Mathematics 8.7.1 and Dan Reynolds' 8.7.1 lecture notes. 
%
% %
% Inputs:   x  nodal locations 
%           f  data values 
%           M  cubic spline coefficients
%           z  evaluation point
% Outputs:  s  s(z) cubic spline evaluate at the evaluation point z
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

%initialize s
s=0;

% using sj formula
for j=1:length(x)-1 %1:n-1 because we're checking between xj and xj+1
    % check what interval z is in
    if le(x(j),z) && le(z,x(j+1))
        h=x(j+1)-x(j);
        s = (1/(6*h))*M(j)*((x(j+1)-z).^3)+((1/(6*h))*M(j+1)*(z-x(j)).^3)+...
        ((z-x(j))*((1/h)*(f(j+1)-f(j))-(h/6)*(M(j+1)-M(j))))+f(j)-((h.^2)/6)*M(j);
        return
    else
    end
    
end

end
       
