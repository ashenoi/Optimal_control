function [T,X,J] = dynamics(u,x0)
global Tf;
global dt;
i = 0;
J = 0;
X = [];
T = [];
x1 = x0(1);
x2 = x0(end);
% X = [x1;x2];
for t = 0:dt:Tf
    i = i+1;
    x1_dot = u(i,1) - sqrt(x1);
    x2_dot = sqrt(x1) - sqrt(x2);
    x1 = x1 + x1_dot*dt;
    x2 = x2 + x2_dot*dt;
    J = J + (x2-3).^2*dt;
%     if x1<=0, x1 = 0;  end;
%     if x2<=0, x2 = 0; end;
    X = [X [x1;x2]];
    T = [T t];
  
end
J = 2*J;
end
