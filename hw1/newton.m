function [x,its] = newton(Ffun, Jfun, x, maxit, rtol, atol, output)
% Usage: [x,its] = newton(Ffun, Jfun, x, maxit, rtol, atol, output)
%
% This routine uses the Newton method to approximate a root of
% the nonlinear system of equations F(x)=0.  The iteration ceases 
% when the following condition is met:
%
%    ||xnew - xold|| < atol + rtol*||xnew||
%
% inputs:   Ffun     nonlinear function name/handle
%           Jfun     Jacobian function name/handle
%           x        initial guess at solution
%           maxit    maximum allowed number of iterations
%           rtol     relative solution tolerance
%           atol     absolute solution tolerance
%           output   flag (true/false) to output iteration history/plot
% outputs:  x        approximate solution
%           its      number of iterations taken
%
% D.R. Reynolds
% Math6316 @ SMU
% Spring 2016

% check input arguments
if (floor(maxit) < 1) 
   fprintf('newton: maxit = %i < 1. Resetting to 10\n',...
           floor(maxit)); 
   maxit = 10;
end
if (rtol < 10*eps)
   fprintf('newton: rtol = %g < %g. Resetting to %g\n',...
           rtol, 10*eps, 10*eps)
   rtol = 10*eps;
end
if (atol < 0)
   fprintf('newton: atol = %g < 0. Resetting to %g\n',...
           atol, 10*eps)
   atol = 10*eps;
end

% initialize error history arrays
abs_err_hist = zeros(maxit,1);
rel_err_hist = zeros(maxit,1);

% begin iteration
for its=1:maxit

   % evaluate function and derivative
   F = feval(Ffun,x);
   J = feval(Jfun,x);

   % compute Newton update, new guess at solution
   h = -J\F;
   x = x + h;

   % check for convergence and output diagnostics
   hnorm = norm(h);
   xnorm = norm(x);
   if (output)
      fprintf('   iter %3i, \t||h|| = %g, \t||h||/||x|| = %g\n',...
              its, hnorm, hnorm/xnorm);
      abs_err_hist(its) = hnorm;
      rel_err_hist(its) = hnorm/xnorm;
   end
   if (hnorm < atol + rtol*xnorm)
      break;
   end

end

% plot convergence history (if requested)
if(output)
   str1=sprintf('Newton convergence history -- absolute error with rtol=%g and atol=%g', rtol, atol);
   str2=sprintf('Newton convergence history -- relative error with rtol=%g and atol=%g', rtol, atol);
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

% end of function
