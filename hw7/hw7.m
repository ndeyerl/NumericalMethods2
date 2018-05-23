% Driver to test various ODE solvers for the test problem
%    y' = -y(t)+2cos(t),   t in [0,10]
%    y(0) = 1.
% Based off of Dan Reynolds' driver1.m
% Nicole Deyerl
% Math 6316, SMU
% Spring 2016
clear

% ODE RHS , initial condition, etc.
f   = @(t,y) -y+2*cos(t);
y0  = 1;
t0  = 0;
Tf  = 10;
k=1:9;
n=10*2.^k;
h=10./n;

    % iterate over our h values
    for j=1:length(h)

           % set up tests for this h value
           tvals = t0:h(j):Tf;
           ytrue = sin(tvals)+cos(tvals); %via wolframalpha

           % initialize outputs
           RK3AB3(1)    = y0;

           % try out our methods
           for i=1:length(tvals)-1
              hstep = tvals(i+1)-tvals(i);
              t = tvals(i);
              if i<=2 %use RK3 for first two steps
                  RK3AB3(i+1)    = RK3_step(f, t, RK3AB3(i), hstep);
              else %then switch to AB3
                  RK3AB3(i+1)    = AB3_step(f, t, RK3AB3(i),...
                                     RK3AB3(i-1),RK3AB3(i-2),hstep);
              end 
           end

        % compute errors
        err = ytrue - RK3AB3;
        err_RK3AB3(j)  = max(abs(err));
    end

% output results
disp('Results for RK3AB3:')
err = err_RK3AB3;
fprintf('   h = %10g,  err = %.2e\n',h(1),err(1))

%compute and print the convergence rates and error
if j==1 %can't compute convergence for j=1
    for i=2:length(h)
         fprintf('   h = %10g,  err = %.2e,  rate = %g\n',...
           h(i),err(i))    
    end
else
    for i=2:length(h)
        fprintf('   h = %10g,  err = %.2e,  rate = %g\n',...
               h(i),err(i),log(err(i)/err(i-1))/log(h(i)/h(i-1)))
    end
end
