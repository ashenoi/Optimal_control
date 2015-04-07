function theta = optimality(X,P,u,u_star)
global dt;
theta = 0;
for i = 1:1:length(X)
    f = @(i,u)[u(i,1)-sqrt(X(1,i)) ; sqrt(X(1,i))-sqrt(X(2,i))];
    
    H1 = P(:,i)'*f(i,u) + 2*(X(2,i) - 3)^2;
    H2 = P(:,i)'*f(i,u_star) + 2*(X(2,i) - 3)^2;
    
    theta = theta + (H2-H1)*dt;
end

end

   