function unew = AB3_step(f, t, u, u1, u2, h)
% usage: unew = AB3_step(f, t, u, u1, u2, h)
%
% Adams Bashforth method of order 3 solver for one step of the scalar-valued 
% ODE problem, 
%     y' = f(t,y), t in tspan,
%     y(t0) = y0.
% Based off Dan Reynolds' RK2_step.
%
% Inputs:  f = function name for ODE right-hand side, f(t,y)
%          t = current time
%          y = current solution
%          h = time step size
% Outputs: unew = updated solution
%
% Nicole Deyerl
% Math 6316, SMU
% Spring 2016

b0=23/12;
b1=-4/3;
b2=5/12;
tnm2=t-2*h;
tnm1=t-h;
fn=feval(f,t,u);
fnm1=feval(f,tnm1,u1);
fnm2=feval(f,tnm2,u2);

% get un+1
unew=u+h*(b0*fn+b1*fnm1+b2*fnm2);

end