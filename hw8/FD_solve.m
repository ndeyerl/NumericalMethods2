function [u,x] = FD_solve(alpha, gamma, f, lambda, mu, eta, theta, ga, gb, a, b, n)
% usage: [u,x] = FD_solve(alpha, gamma, f, lambda, mu, eta, theta, ga, gb, a, b, n)
% 
% This is a function which solves a reaction diffusion equation in 1 dimension
% with Robin mixed boundary conditions.  It solves the left endpoint using
% a ghost node, and the right endpoint using a one sided finite difference 
% method.  It is based off of Dan Reynolds 1d heat equation script.
% 
% Inputs:  alpha = pde coefficient function alpha(x)
%          gamma = pde coefficient function gamma(x)
%          f     = pde right hand side f(x)
%          lambda= robin condition coefficient (left endpoint)
%          mu    = robin condition coefficient (left endpoint)
%          eta   = robin condition coefficient (right endpoint)
%          theta = robin condition coefficient (right endpoint)
%          ga    = robin condition rhs (left endpoint)
%          gb    = robin condition rhs (right endpoint)
%          n     = number of meshpoints
% Outputs: u     = solution to the pde
%          x     = nodes
% Nicole Deyerl
% Math 6316, SMU
% Spring 2016


% set spatial mesh size
h = (b-a)/n;

% set up nodes and dual nodes
x=linspace(a,b,n+1);
%xd=linspace(a+h/2,b+h/2,n+1);

% create matrix storage
A = sparse(n+1,n+1);

% interior rows and RHS
for i=2:n
	A(i,i-1) = alpha(x(i)-h/2);
	A(i,i) = -(alpha(x(i)-h/2)+alpha(x(i)+h/2)+(h.^2)*gamma(x(i)));
	A(i,i+1) = alpha(x(i)+h/2);
    F(i)= -(h.^2)*feval(f,x(i));
end

% boundary rows
A(1,1) = (2/(h.^2))*alpha(x(1)+h/2)-((2*lambda)/(h*mu))+gamma(x(1));
A(1,2) = -(2/(h.^2))*(alpha(x(1)+h/2));
A(n+1,n-1) = theta*alpha(x(n+1))/(2*h);
A(n+1,n) = -(4*theta*alpha(x(n+1))/(2*h));
A(n+1,n+1) = eta+(3*theta*alpha(x(n+1))/(2*h));

%boundary RHS
F(1)=feval(f,x(1))-(2/(h*mu))*ga;
F(n+1)=gb;
F=F';

u=A\F;

end