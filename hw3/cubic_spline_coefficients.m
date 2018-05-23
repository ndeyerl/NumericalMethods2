function M = cubic_spline_coefficients(x, f)
% Usage: M = cubic_spline_coefficients(x, f)
%
% This routine uses a matrix to solve for the M cubic spline coefficients
% representing the second derivatives of the cubic spline, at the nodes x
% for the values f.  After using this function one may input the coefficients
% M into the function cubic_spline_evaluate.m to use the coefficients to 
% interpolate f at the nodes x. It is based off of problem 3 of homework 3,
% Quarteroni's Numerical Mathematics 8.7.1 and Dan Reynolds' 8.7.1 lecture notes.  
% %
% Inputs:   x  nodal locations 
%           f  data values 
% Outputs:  M  cubic spline coefficients
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

% check inputs
if (length(x) ~= length(f)) 
   error('Cubic spline error: (x,f) have different sizes');
end

if (f(1) ~= f(end)) 
   error('Periodic boundary condition error: f0 and fn are not equal');
end

np1=max(size(x));
n=np1-1; % x is size 'n+1'

% initialize coefficients running from 0:n-1 for n total 
mu=zeros(n,1);
lam=zeros(n,1);
d=zeros(n,1);

% initialize matrix A of size 'n+1'
A=zeros(np1,np1);

% loop to build coefficients based on lecture notes 8.7.1
for j=2:n
    hj=x(j)-x(j-1);
    hjpo=x(j+1)-x(j);
    
    mu(j-1)=hj/(hj+hjpo);
    lam(j-1)=hjpo/(hj+hjpo);
    d(j-1)=(6/(hj+hjpo))*( ((f(j+1)-f(j))/hjpo) + ((f(j-1)-f(j))/hj));
end

% build nth coefficients since the index was shifted backward to start
% the coefficients saving from 0

mu(n)=hj/(hj+hjpo); % can reuse same hj, hjpo since the last hj, hjpo saved
lam(n)=hjpo/(hj+hjpo);% were hj(n), hjpo(n) 
d(n+1)=(6/(hj+hjpo))*( ((f(n+1)-f(n))/hjpo) + ((f(n-1)-f(n))/hj));


% build tridiagonal matrix A from coefficients
for j=1:n-1
    A(j,j)=mu(j);
    A(j,j+1)=2;
    A(j,j+2)=lam(j);
end

% first boundary equation coming from s''(a)=s''(b)
A(np1-1,1)=1;
A(np1-1,np1)=-1;
d(np1-1)=0;

hone=x(2)-x(1);
hn=x(np1)-x(np1-1);

% second boundary equation coming from s'(a)=s'(b)
A(np1,1)=(hone/3);
A(np1,2)=(hone/6);
A(np1,np1-1)=(hn/6);
A(np1,np1)=(hn/3);
d(np1)=((1/hone)*(f(2)-f(1)))+((1/hn)*(f(n-1)-f(n)));

% solve the matrix using backwards substitution for the 
% cubic spline coefficients Mn
M=A\d;
end
