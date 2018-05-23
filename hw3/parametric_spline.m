function v = parametric_spline(x, m)
% Usage: v = parametric_spline(x, m)
%
% This routine is meant to calculate a parametric spline for a 2 column 
% vector of (x,y) locations in R^2.  It uses the functions
% cubic_spline_coefficients and cubic_spline_evaluate to build an interpol-
% -ating parametric curve at m points.  The program is based off of hints 
% from Dan Reynolds MATH 6316 homework 3.pdf.
% %
% Inputs:   x  nodal locations in R^2
%           m  integer equal to at least 1, specifies number of eval points
% Outputs:  v  output coordinates for an interpolating curve
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

n=max(length(x));

xcoord=x(1,:); % take the first column of input x as the xk's
ycoord=x(2,:); % take the second column of input x as the yk's
l=[];
l(1)=0; % initialize first line segment length as 0

% loop over line segment lengths as in homework3.pdf
for j=2:n
    l(j)=sqrt( (xcoord(j)-xcoord(j-1)).^2+(ycoord(j)-ycoord(j-1)).^2 );
end

% initialize tk
t=zeros(n,1);

% loop over line segment length sums as in homework3.pdf
for k=1:n
    for j=1:k
        t(k)=t(k)+l(j);
    end
end

% call the spline coefficients for the separate coordinates against 
% the new knots tk
Mx = cubic_spline_coefficients(t,xcoord);
My = cubic_spline_coefficients(t,ycoord);

% initialize m data times
z = linspace(0,t(n),m);

% initialize vectors to hold the splines
sy=[];
sx=[];

% loop over the m data times to assemble the cubic splines for each 
% coordinate vector
for k=1:length(z)
    sx(k) = cubic_spline_evaluate(t, xcoord, Mx, z(k));
    sy(k) = cubic_spline_evaluate(t, ycoord, My, z(k));
end
% return matrix v
v=[sx;sy];

end

