
function [u,v] = doublegyre(x,y,t)
A = 0.25;
epsilon = 0.25;
omega = 2.0 * pi / 10.0;

u = calcU(x,y,t, A, epsilon, omega);
v = calcV(x,y,t, A, epsilon, omega);
end

