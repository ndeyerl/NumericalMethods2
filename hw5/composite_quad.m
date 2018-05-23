function [Imn, nf] = composite_quad(f, a, b, m)
% Usage: [Imn, nf] = composite_quad(f, a, b, m)
%
% This code approximates the integral of input f(x) on the interval from 
% a to b, using a composite Gaussian quadrature formula with m uniformly 
% sized subintervals.  This quadrature formula is a 3-point formula which 
% forces one of the nodes to equal 1.  This code was motivated by homework
% 5 and based off of Dan Reynolds' composite_quadrature.m and Gauss3.m files.
%
% Inputs:  f        integrand function
%          a        left end point of integration
%          b        right end point of integration
%          m        number of subintervals
%
% Outputs: Imn      value of the composite quadrature result
%          nf       total number of calls to f
%
% Nicole Deyerl
% Math6316
% Spring 2016

% check inputs
if (b < a) 
   error('composite_quadrature error: b < a!');
end
if (m < 1) 
   error('composite_quadrature error: m < 1!');
end

% initialize results
Imn = 0;
nf = 0;

%define weights
w = [(8/9)-(1/(3*sqrt(6))), (1/18)*(16+sqrt(6)), (2/9)];

% iterate over subintervals
for i=1:m
    
   % create new endpoints based on number of subintervals
   leftend=a+(i-1)*(b-a)*(1/m);
   rightend=a+i*(b-a)*(1/m);
   
   % determine 3 nodes within the subinterval based on xk's we calculated
   % using legendre polynomials
   x = (leftend+rightend)/2 + (rightend-leftend)/2*[-(1/5)-(sqrt(6)/5), -(1/5)+(sqrt(6)/5), 1];
   
   % call quadrature formula on this subinterval
   In = (rightend-leftend)/2*(sum(w.*feval(f,x)));
   nflocal=3; % 3 calls to f since we're doing 3-point Gaussian Quadrature

   % increment outputs
   Imn = Imn + In;
   nf  = nf  + nflocal;

end
