% % hw8.m
% 
% This is a script used to run FD_solve, which is a solver for a reaction 
% diffusion equation in 1 dimension with robin boundary conditions.  This 
% script tests the accuracy of FD_solve using different mesh sizes.
% 
% Nicole Deyerl
% Math 6316, SMU
% Spring 2016

clear

% problem specification
k=[1:8];
n = 10*2.^k;   % mesh sizes
a = 0;   % [a,b] domain
b = 2*pi;   
f = @(x) (4*x*(sin(2*x)+2*cos(2*x))+6*(1+x.^2)*(cos(2*x)-2*sin(2*x))); % pde rhs 
alpha = @(x) (1+x.^2); %variable pde coefficient
gamma = @(x) (2+2*x.^2); %variable pde coefficient

%robin condition coefficients
lambda=2; %left endpoint
mu=1; %left endpoint
eta=1; %right endpoint
theta=1/(-1-4*pi.^2); %right endpoint

%robin condition rhs's
gb=5; % right rhs
ga=-2; %left rhs

% true solution
U = @(x) cos(2*x)-2*sin(2*x);

% loop over mesh sizes
for j=1:length(n)
    
    %evaluate the approximated solution
    [u,x] = FD_solve(alpha, gamma, f, lambda, mu, eta, theta, ga, gb, a, b, n(j));
    
    %evaluate the true solution at the nodes
    utrue=feval(U,x)';
    
    % plot true solution and approximation
    figure()
    plot(u)
    hold on
    plot(utrue,'r')
    % label axes, legend and title
    xlabel('x'), ylabel('u')
    legend('u(x)','utrue(x)')
    title(['Reaction Diffusion Equation Solutions for mesh size n = ' num2str(n(j))])

    % compute errors
    err = utrue - u;
    err_u(j) = max(abs(err));

end

h=(b-a)./n;

% output results
fprintf('Results for FD_solve: \n')
%print results for i=1 since we can't compute convergence for i=1
fprintf('   h = %10g,  err = %.2e \n',...
           h(1),err_u(1))    
       
%compute and print the convergence rates and error
for i=2:length(h)
    fprintf('   h = %10g,  err = %.2e,  rate = %g\n',...
           h(i),err_u(i),log(err_u(i)/err_u(i-1))/log(h(i)/h(i-1)))

end
