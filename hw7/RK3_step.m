function unew = RK3_step(f, t, u, h)
% usage: unew = RK3_step(f, t, u, h)
%
% Runge-Kutta of order 3 solver for one step of the scalar-valued ODE
% problem, 
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

% get K1, K2, K3
K1=feval(f,t,u);
K2=feval(f, t+h/2, u+(1/2)*h*K1);
K3=feval(f, t+h, u-h*K1+2*h*K2);

% get un+1
unew=u+h*((1/6)*K1+(2/3)*K2+(1/6)*K3);

end