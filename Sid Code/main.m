clc;
clear all;
close all;

constants;
global beta;
for i = 1:1:20
    [T,X,J] = dynamics(u,x0);
    
    [~,P] = costate(T,X,u);
    u_star = min_Hamilt(P);
    theta = optimality(X,P,u,u_star);
    lambda = armijo(theta,J,u,u_star);
    u = (1-lambda).*u + (lambda).*u_star;
    J_stor = [J_stor J];
    THETA = [THETA theta];
%     K = [K k];
end

plot(J_stor);
title('Cost vs iterations');
figure(2)
plot(T,u);
title('u');
figure(3)
plot(T,u_star);
title('u*');
figure(3);
plot(T,X);
title('Water levels');

    