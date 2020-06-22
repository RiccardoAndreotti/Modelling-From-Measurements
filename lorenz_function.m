function [t,y] = lorenz_function(x0,r)

dt=0.01; T=8; t=0:dt:T; 
b=8/3; sig=10; 

Lorenz = @(t,x)([ sig * (x(2) - x(1))       ; ...
                  r * x(1)-x(1) * x(3) - x(2) ; ...
                  x(1) * x(2) - b*x(3)         ]);              

ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);

[t,y] = ode45(Lorenz,t,x0); 

end

