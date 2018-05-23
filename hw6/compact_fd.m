function [ddf,x] = compact_fd(f,a,b,n)
% Usage: [ddf,x] = compact_fd(f, a, b, n)
% 
% This is a routine used to approximate the second derivative of a function
% with periodic boundary conditions based on the layout
% h^2[alphm1*u(i-1)+u(i)+alph1*u(i+1)]=betm1*f(i-1)+b0*f(i)+bet1*f(i+1)
% It is an order 4 method.
% 
% Inputs:   f    function handle
%           a    left endpoint of interval
%           b    right endpoint of interval
%           n    number of unique nodes to use (note x0=xn)
% Outputs:  ddf  2nd derivative approximation of f
%           x    the unique nodes at which f was evaluated for its deriv.
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016


x=linspace(a,b,n+1);
%remove one node for being nonunique
x=x(2:n+1);

% form matrix A
A=diag((1/10)*ones(n-1,1),-1)+diag((1/10)*ones(n-1,1),1)+diag(ones(n,1));
A(1,n)=1/10;
A(n,1)=1/10;


% form right hand side b
rhs=zeros(n,1);
beta=6/5;
gamma=12/5;
h=(b-a)/n;
for i=2:n-1
    rhs(i)=(1/(h.^2))*(beta*feval(f,x(i-1))-gamma*feval(f,x(i))+beta*feval(f,x(i+1)));
end
% code the first and last equations separately because they
% contain the periodic nodes/end and beg of period
rhs(1)=(1/(h.^2))*(beta*feval(f,x(n))-gamma*feval(f,x(1))+beta*feval(f,x(2)));
rhs(n)=(1/(h.^2))*(beta*feval(f,x(n-1))-gamma*feval(f,x(n))+beta*feval(f,x(1)));

% linsolve for the solution vector
ddf=A\rhs;

end

