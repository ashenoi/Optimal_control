function  lambda = armijo(theta,J,u,u_star)
global dt;
global Tf;
global beta;
global alpha;
J_a = 0;
X_a = [];
i = 1;
x_a = [2;2];
lambda = 0.001;
for k = 1:1:40
    u_new = (1-beta^k).*u + beta^k.*u_star;
    [T,X_a,J_a] = dynamics(u_new,x_a);
    if (J_a - J) <= alpha*beta^k*theta;
        lambda = beta^k;
        break;
    end
end

end

