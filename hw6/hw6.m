% hw6.m
%
% This is a script to test compact_fd.m to approximate the second
% derivative of a periodic function f on the interval [0,4pi] using
% different numbers of nodes.  It will also compare the error between the
% real and approximate second derivative for different numbers of nodes.
% 
% Nicole Deyerl
% Math6316
% Spring 2016

clear;

% define function f(x)=e^cos(x) and true second derivative
f=@(x) exp(cos(x));
ddftrue=@(x) ((sin(x).^2)-cos(x)).*exp(cos(x));

% boundaries
a=0;
b=4*pi;

% node amounts
n=zeros(1,9);
for j=0:8
    n(j+1)=5*2.^j;
end

%initialize error and convergence
error=zeros(1,9);
convergence=zeros(1,8);

%loop over n values to calculate second deriv. of f and compare with true
%second deriv. of f
for j=1:9
    if j==1
        fprintf('\n n = %i\n',n(j));
        [ddf,x] = compact_fd(f,a,b,n(j));
        ddft=ddftrue(x)'; % calculate true second deriv @ nodes
        error(1,j)=max(abs(ddf-ddft)); % calculate error between true and 
        fprintf('  error = %9.3e', error(j)) % approx soln
    else
        fprintf('\n n = %i\n',n(j));
        [ddf,x] = compact_fd(f,a,b,n(j));
        ddft=ddftrue(x)'; % calculate true second deriv @ nodes
        error(1,j)=max(abs(ddf-ddft));
        % calculate convergence, based off of homework 5 convergence
        % formula
        convergence=(log(error(j))-log(error(j-1)))/(log(n(j))-log(n(j-1)));
        fprintf('  error = %9.3e,  convergence = %9.3e', error(j), convergence)
    end
end

