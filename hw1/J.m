function Jac = J(u)
% Usage: Jac = J(u)
%
% This routine forms the Jacobian matrix for the nonlinear residual
% function associated with the semiconductor device simulation
% problem from section 7.4.1.
%
% Daniel R. Reynolds
% Math6316 @ SMU
% Spring 2016

% access globally-stored variables
global K A;
Jac = A + 2*K*diag(cosh(u));

% end of function
