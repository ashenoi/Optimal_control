function robot_path_planning()
close all;
global tf;
global x0;
global pf;
global dt;
global alpha;
global beta;
global rho;
global x_goal;
global x_obs;

alpha = 2;
beta = 0.1;
rho = 0.01;
x_goal = [ 0 ; 2.5];
x_obs = [-0.5 ; 1.5];

format long g
tf=3;
x0=[0.05;0];
pf=[0;0];
dt=0.001;
T=0:dt:tf;

max_iter = 10;
iter_num = 0;
N=length(T);
U=zeros(2,N);
a = 0.6;
b = 0.5;
J_array = [];
last_j=2;
while(iter_num<max_iter)
    j = last_j-2;
    X = compute_x(U);
    P = compute_p(X,U);
    J = compute_cost(X,U)
    H = compute_h(P,U,X);
    U_star = compute_u_star(P,X);
    H_star = compute_h(P,U_star,X);
    theta = compute_theta(H_star,H);
    J_array = [J_array,J];
    U1 = ((1-b^j)*U + b^j*U_star);
    X1 = compute_x(U1);
    J1 = compute_cost(X1,U1);
    while((J1-J)> a*b^j*theta)
        j = j+1;
        U1 = ((1-b^j)*U + b^j*U_star);
        X1 = compute_x(U1);
        J1 = compute_cost(X1,U1);  
    end
    last_j=j;
    U=U1;
    iter_num=iter_num +1;
end


figure;
plot(T,U);
figure;
plot(X(1,:),X(2,:));
figure;
plot(J_array);

disp('Done');
end

function fg = fg(x)
global x_goal;
fg = (x_goal - x)/norm(x_goal - x);
end

function fcl = fcl(x)
global x_obs;
fcl = [0 -1;1 0]*((x_obs - x)/norm(x_obs - x));
end

function fcc = fcc(x)
global x_obs;
fcc = [0 1;-1 0]*((x_obs - x)/norm(x_obs - x));
end

function L1 = L1(x1,u1)
global alpha;
global beta;
global rho;
global x_goal;
global x_obs;

L1 = rho*norm(x_goal-x1)^2 + alpha*exp(-(norm(x_obs-x1)^2)/beta);
end

function fg_prime = fg_prime(x1)
global x_goal;
fg_prime = ((-(norm(x_goal-x1))*eye(2)) +((x_goal - x1)*(x_goal-x1)')/(norm(x_goal-x1)))/(norm(x_goal-x1))^2;
end

function fcl_prime = fcl_prime(x1)
global x_obs;
fcl_prime = [0 -1;1 0]*((-(norm(x_obs-x1))*eye(2)) +((x_obs - x1)*(x_obs-x1)')/(norm(x_obs-x1)))/(norm(x_obs-x1))^2;
end

function fcc_prime = fcc_prime(x1)
global x_obs;
fcc_prime = [0 1;-1 0]*((-(norm(x_obs-x1))*eye(2)) +((x_obs - x1)*(x_obs-x1)')/(norm(x_obs-x1)))/(norm(x_obs-x1))^2;
end


function L1_prime = L1_prime(x1,u1)
global x_goal;
global x_obs;
global alpha;
global beta;
global rho;

L1_prime =-(2*rho*(x_goal-x1)' + 2*(alpha/beta)*exp(-(norm(x_obs-x1)^2)/beta)*(x_obs-x1)') ;
end

function X1=compute_x(U1)
global dt;
global x0;

X1=[];
x= x0;
N =length(U1);
for i=1:N
    X1=[X1,x];
    x=x+dt*( (U1(1,i))*fg(x) + U1(2,i)*fcl(x) +(1-U1(1,i)-U1(2,i))*fcc(x));    
end
end

function P1=compute_p(X1,U1)

global dt;
global pf;
global r;
N=length(U1);
P1=[];
p=pf;

for i=1:N
    P1=[p,P1];
    p=p-dt*(-((U1(1,i))*fg_prime(X1(:,i)) + U1(2,i)*fcl_prime(X1(:,i)) + (1-U1(1,i)-U1(2,i))*fcc_prime(X1(:,i)))'*p -(L1_prime(X1(:,i),U1(:,i)))');
end
end

function J1 = compute_cost(X1,U1)
global dt;
global r;
J1= 0;
N=length(U1);
for i =1:N
    J1=J1+dt*(L1(X1(:,i),U1(:,i)));
end
end

function H1 = compute_h(P1,U1,X1)
global r;
N=length(U1);
H1=[];
for i=1:N
    h=P1(:,i)'*((U1(1,i))*fg(X1(:,i)) + U1(2,i)*fcl(X1(:,i)) + (1-U1(1,i)-U1(2,i))*fcc(X1(:,i)))+L1(X1(:,i),U1(:,i));
    H1=[H1,h];
end
end

function U1_star = compute_u_star(P1,X1)
N=length(X1);
U1_star =[];
u= [0;0];
for i=1:N
    if P1(:,i)'*fg(X1(:,i)) < P1(:,i)'*fcl(X1(:,i)) && P1(:,i)'*fg(X1(:,i)) < P1(:,i)'*fcc(X1(:,i))
        u(1)=1;
        u(2)=0;
    elseif P1(:,i)'*fcl(X1(:,i)) < P1(:,i)'*fg(X1(:,i)) && P1(:,i)'*fcl(X1(:,i)) < P1(:,i)'*fcc(X1(:,i))
        u(1)=0;
        u(2)=1;       
    elseif P1(:,i)'*fcc(X1(:,i)) < P1(:,i)'*fg(X1(:,i)) && P1(:,i)'*fcc(X1(:,i)) < P1(:,i)'*fcl(X1(:,i))
        u(1)=0;
        u(2)=0;       
    end

    U1_star=[U1_star,u];
end

end

function theta1 = compute_theta(H1_star,H1)
global dt;
N=length(H1);
theta1=sum((H1_star-H1))*dt;
end
