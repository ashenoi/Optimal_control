x0 = [2;2];

global Tf;
Tf = 10;
global dt;
dt = 0.01;
u = ones(Tf/dt + 1,1);
u_star = zeros(Tf/dt,1);
global beta;
beta = 0.5;
global alpha;
alpha = 0.5;
J_stor = [];
THETA = [];
K = [];