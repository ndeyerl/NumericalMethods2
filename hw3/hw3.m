% Script to test parametric cubic spline interpolation by making 3 
% parametric plots.
%
% Nicole Deyerl
% Math6316 @ SMU
% Spring 2016

% remove all existing variables
clear 

x=[1,2,3,4,5,6,4,3,2,1]; %10 points
y=[1,4,7,4,1,-2,-4,-7,-4,1];
x=[x;y];
m=2000; % 2000 evaluation points

plot(x(1,:),x(2,:),'ro','MarkerFaceColor','r','MarkerSize',10)
hold on
v = parametric_spline(x, m);
plot(v(1,:),v(2,:),'b')
 
x=[1,2,3,4,5,6,5,4,5,6,7,8,1]; %13 points
y=[0,1,3,6,10,11,7,5,4,3,2,1,0];
x=[x;y];
m=2000; % 2000 evaluation points

figure %new plot
plot(x(1,:),x(2,:),'ro','MarkerFaceColor','r','MarkerSize',10)
hold on
v = parametric_spline(x, m);
plot(v(1,:),v(2,:),'b')
 
x=[1,2,3,4,5,6,7,6,5,4,3,2,1]; %13 points
y=[0,1,3,6,10,11,7,10,12,13,5,3,0];
x=[x;y];
m=2000; % 2000 evaluation points

figure %new plot
plot(x(1,:),x(2,:),'ro','MarkerFaceColor','r','MarkerSize',10)
hold on
v = parametric_spline(x, m);
plot(v(1,:),v(2,:),'b')
