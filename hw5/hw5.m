% hw5.m
%
% This is a script to test composite_quad.m using the 3-point Gaussian
% quadrature which forces 1 endpoint, and compare it with the files Gauss3.m
% and Simpson.m by computing absolute error and outputting convergence rates
% for increasing numbers of intervals.  It is based off of Dan Reynolds' 
% composite_test.m file. 
% 
% Nicole Deyerl
% Math6316
% Spring 2016

clear;
% set the integration interval, integrand function and interval amounts
a = 0; 
b = pi;
f = @(x) exp(2.*x).*cos(3.*x);
m=zeros(1,7);
for j=0:6
    m(j+1)=4*2^j;
end

% set the true integral
Itrue = 1/13*(exp(2*b)*(3*sin(3*b)+2*cos(3*b)))-1/13*(exp(2*a)*(3*sin(3*a)+2*cos(3*a)));
fprintf('\n I = %16.10e\n',Itrue)

%initialize error vectors to use for convergence 
errorforced=zeros(1,6);
errorgauss=zeros(1,6);
errorsimp=zeros(1,6);

% iterate over a variety of m values
for i=1:7
   if i==1
       
       fprintf('\n m = %i\n',m(i));

       % 3-point Gaussian Quad with 1 forced end
       fprintf('Gauss3 Forced: ')
       [Iapprox,nf] = composite_quad(f, a, b, m(i));
       errorforced(i) = abs(Itrue-Iapprox)/abs(Itrue); %absolute error
       fprintf(' Imn = %16.10e,  error = %9.3e,  #f = %i\n', Iapprox, errorforced(i), nf)

       % Gauss3 Test
       fprintf('       Gauss3:  ')
       [Iapprox,nf] = composite_quadrature(f, a, b, 'Gauss3', m(i));
       errorgauss(i) = abs(Itrue-Iapprox)/abs(Itrue); %absolute error
       fprintf('Imn = %16.10e,  error = %9.3e,  #f = %i\n', Iapprox, errorgauss(i), nf)

       % Simpson test
       fprintf('      Simpson:  ')
       [Iapprox,nf] = composite_quadrature(f, a, b, 'Simpson', m(i));
       errorsimp(i) = abs(Itrue-Iapprox)/abs(Itrue); %absolute error
       fprintf('Imn = %16.10e,  error = %9.3e,  #f = %i\n', Iapprox, errorsimp(i), nf)

   else

       fprintf('\n m = %i\n',m(i));

       % 3-point Gaussian Quad with 1 forced end
       fprintf('Gauss3 Forced: ')
       [Iapprox,nf] = composite_quad(f, a, b, m(i));
       errorforced(i) = abs(Itrue-Iapprox)/abs(Itrue); %absolute error
       convergence=(log(errorforced(i))-log(errorforced(i-1)))/(log(m(i))-log(m(i-1)));
       fprintf(' Imn = %16.10e,  error = %9.3e,  convergence = %9.3e, #f = %i\n', Iapprox, errorforced(i), convergence, nf)

       % Gauss3 Test
       fprintf('       Gauss3:  ')
       [Iapprox,nf] = composite_quadrature(f, a, b, 'Gauss3', m(i));
       errorgauss(i) = abs(Itrue-Iapprox)/abs(Itrue); %absolute error
       convergence=(log(errorgauss(i))-log(errorgauss(i-1)))/(log(m(i))-log(m(i-1)));
       fprintf('Imn = %16.10e,  error = %9.3e,  convergence = %9.3e, #f = %i\n', Iapprox, errorgauss(i), convergence, nf)

       % Simpson test
       fprintf('      Simpson:  ')
       [Iapprox,nf] = composite_quadrature(f, a, b, 'Simpson', m(i));
       errorsimp(i) = abs(Itrue-Iapprox)/abs(Itrue); %absolute error
       convergence=(log(errorsimp(i))-log(errorsimp(i-1)))/(log(m(i))-log(m(i-1)));
       fprintf('Imn = %16.10e,  error = %9.3e,  convergence = %9.3e, #f = %i\n', Iapprox, errorsimp(i), convergence, nf)
   end

end

% plot the integrand over this interval
xvals = linspace(a, b, 2000);
fvals = f(xvals);
plot(xvals, fvals)
xlabel('x'), ylabel('y'), title('Integrand')

% end of script
