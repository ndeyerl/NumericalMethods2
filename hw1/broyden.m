function [x,its]=broyden(Ffun,Jfun,x0,maxit,atol,rtol,output)
% Usage: [x,its]=broyden(Ffun,Jfun,x0,maxit,atol,rtol,output)
%
% This routine uses the Broyden method to approximate a root of
% the nonlinear system of equations F(x)=0.  The iteration
% ceases if the maximum number of iterations has been met, or if 
% the norm on s, the update value, is less than 
% atol+rtol*||x||.
%
% inputs:   Ffun     nonlinear function name/handle
%           Jfun     Jacobian function name/handle
%           x0        initial guess at solution
%           maxit    maximum allowed number of iterations
%           atol     absolute solution tolerance
%           rtol     relative solution tolerance
%           output   flag (true/false) to output iteration history/plot
% outputs:  x        approximate solution
%           its      number of iterations taken
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016


% check input arguments
if (floor(maxit) < 1) 
   fprintf('broyden: maxit = %i < 1. Resetting to 10\n',...
           floor(maxit)); 
   maxit = 10;
end
if (rtol < 10*eps)
   fprintf('broyden: rtol = %g < %g. Resetting to %g\n',...
           rtol, 10*eps, 10*eps)
   rtol = 10*eps;
end
if (atol < 0)
   fprintf('broyden: atol = %g < 0. Resetting to %g\n',...
           atol, 10*eps)
   atol = 10*eps;
end

%check size of function matrix, exit if not square
[n,m]=size(Ffun);
if n~=m
    error('Only square systems');
end

% initialize error history arrays
abs_err_hist = zeros(maxit,1);
rel_err_hist = zeros(maxit,1);

its=0; %initialize number of iterations at 0
err=1; %set error to 1 so that the while loop can initialize
fk=zeros(n,1); 
fk1=fk; 
x=x0; %initial guess

J=feval(Jfun,x); %evaluate jacobian of F at x0
fk=feval(Ffun,x); %evaluate F at x0

Q=J; %set Q to the jacobian of F at x0

xnorm = norm(x); %dummy norm of x to initialize the while loop

%broyden loop
while its<maxit && err>atol+rtol*xnorm %exit if iterations too large or error small enough
    s=-Q\fk; %h from newton.m, solve by gaussian elim
    x=s+x; %update
    err=norm(s); %check error, same as hnorm
    xnorm=norm(x);
    
    if (output)
        fprintf('   iter %3i, \t||h|| = %g, \t||h||/||x|| = %g\n',...
              its+1, err, err/xnorm); %I added one to its because otherwise...
        abs_err_hist(its+1) = err;    %the number of the iterations is off by 1
        rel_err_hist(its+1) = err/xnorm;
    end
	if err > atol+rtol*xnorm
        fk1=feval(Ffun,x); %evaluate f
        Q=Q+1/(s'*s)*fk1*s'; %update Q
	end
its=its+1; %update number of iterations
fk=fk1; %replace old (kth) f with new (k+1th) f
end

% plot convergence history (if requested)
if(output)
   str1=sprintf('Broyden convergence history -- absolute error with rtol=%g and atol=%g', rtol, atol);
   str2=sprintf('Broyden convergence history -- relative error with rtol=%g and atol=%g', rtol, atol);
   figure()
   subplot(2,1,1)
   semilogy(1:its, abs_err_hist(1:its), 'b-', ...
            1:its, atol*ones(its,1), 'k:')
   xlabel('iteration'), title(str1)
   subplot(2,1,2)
   semilogy(1:its, rel_err_hist(1:its), 'b-', ...
            1:its, rtol*ones(its,1), 'k:')
   xlabel('iteration'), title(str2)
end

end