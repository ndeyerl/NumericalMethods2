% Script to test various solvers on the semiconductor device simulation
% problem from section 7.4.1.
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

% remove all existing variables
clear

% general parameters for all tests
maxit  = 20;
atol   = [1e-6,1e-8,1e-10,1e-12];
rtol   = [1e-5,1e-7,1e-9,1e-11];


% Set up problem
global K A;      % declare some global variables for access by Jacobian function
n = 2000;        % discretization size -- should be even (for how b is set up below)
K = 6.77e-6;
h = 1/(n+1);
phi = @(u) 2*K*sinh(u);
b = [-ones(n/2,1); ones(n/2,1)];
u0 = zeros(n,1);                        % initial guess vector


% set up storage for overall run statistics
% ntests = 2+2*length(alphas)+2*length(pvals);
% stats = zeros(ntests,2);
tests = 0;

%---- run with lambda value .049 for homework 1  ----%
lambda = 0.049;
A = (lambda/h)^2*(2*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1));
F = @(u) A*u + phi(u) - b;              % nonlinear residual function
disp('Tests with lambda = 0.049:')

% run Newton solver 4 times for different rtol, atol vals
disp('Newton tests')
for j=1:4    
    str3=sprintf('rtol=%g and atol=%g',rtol(j),atol(j));
    disp(str3)
    st = tic;
    [u,its] = newton(F, 'J', u0, maxit, rtol(j), atol(j), true);
    ttot = toc(st);
    fprintf('   solution time = %g\n', ttot);
    statsn(j,:) = [its, ttot];
end

% run Broyden solver 4 times for different rtol, atol vals
disp('Broyden tests')
for j=1:4
    str3=sprintf('rtol=%g and atol=%g',rtol(j),atol(j));
    disp(str3)
    st = tic;
    [u,its] = broyden(F, 'J', u0, maxit, rtol(j), atol(j), true);
    ttot = toc(st);
    fprintf('   solution time = %g\n', ttot);
    statsb(j,:) = [its, ttot];
end

%---- output final stats summarizing run time and iteration counts----%
fprintf(' ======================================================\n');
fprintf('  solver                                             | its | runtime\n');
fprintf(' ------------------------------------------------------\n');
for i=1:4
fprintf('  Newton:  lambda=0.049, rtol=%g, atol=%g      | %3i | %g\n',rtol(i),atol(i),statsn(i,:));
fprintf('  Broyden: lambda=0.049, rtol=%g, atol=%g,     | %3i | %g\n',rtol(i),atol(i),statsb(i,:));
end
fprintf(' ======================================================\n');


% end of script