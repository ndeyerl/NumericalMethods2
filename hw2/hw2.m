% Script to test newton and hermite birkoff polynomial interpolation 
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

% remove all existing variables
clear

% declare the function and derivative of the function we want to interpolate
f=@(x) atan(2*x.^2); 
fprime=@(x) 4*x./(4*(x.^4)+1);
L=3; % absolute value of bound value on x and z
m=201; % number of evaluation points
z=linspace(-L,L,m); % evaluation points for polynomial interpolation
fvals=f(z); % f evaluated at z, needed for plots
nvec=[5, 11, 21, 41]; % n values (polynomial degrees)

for n=nvec
    xn=linspace(-L,L,n+1); % newton nodes
    yn=f(xn); % value of f at each of the nodes

    xh=linspace(-L,L,(1/2)*(n+1)); % hermite nodes
    yh=f(xh); % value of f at each of the nodes
    dy=fprime(xh); % value of dy at each of the nodes
    
    for i=1:length(z) % using newton, determine p or PIn of f at each z value
        pn(i) = Newton_interp(xn, yn, z(i));
    end
    
    for i=1:length(z) % using hermite, determine p or PIn of f at each z value
        ph(i)= Hermite_interp(xh,yh,dy,z(i));
    end
    
    % polynomial plot
    figure('units','normalized','position',[.1 .1 .6 .4])
    plot(z, fvals, 'k--', z, pn, 'r--', z, ph, 'b-')
    legend('f(xi)','Newton interpolant','Hermite interpolant')
    xlabel('x'), ylabel('y')
    title(sprintf('Newton and Hermitian interpolant: f vs p_{%i}f, H_{%i}',n,n))
    hold on
    
    % error plot
    figure('units','normalized','position',[.1 .1 .6 .4])
    semilogy(z, abs(fvals-pn), 'r--', z, abs(fvals-ph), 'b-')
    legend('Newton interpolant error','Hermite interpolant error')
    xlabel('x'), ylabel('y')
    title(sprintf('Newton and Hermitian error: E_{%i}f',n))

end
  