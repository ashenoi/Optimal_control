function [X,U,J_array,U_star,T] = double_tank_problem_proj()
global tf;
global x0;
global pf;
global dt;
global r;
format long g
tf=10;
x0=[2.000;2.0000];
pf=[0;0];
dt=0.01;
T=0:dt:tf;
r = 2.5+T*0.1;
%r = 3*ones(1,length(T));
max_iter = 500;
iter_num = 0;
N=length(T);
U=zeros(1,N);
alpha = 0.6;
beta = 0.5;
J_array = [];
last_j=2;
while(iter_num<max_iter)
    if last_j -2 < 0
        j = 0;
    else
        j = last_j-2;
    end
    X = compute_x(U);
    P = compute_p(X,U);
    J = compute_cost(X,U)
    H = compute_h(P,U,X);
    U_star = compute_u_star(P,X);
    H_star = compute_h(P,U_star,X);
    theta = compute_theta(H_star,H);
    J_array = [J_array,J];
    U1 = ((1-beta^j)*U + beta^j*U_star);
    X1 = compute_x(U1);
    J1 = compute_cost(X1,U1);
    while((J1-J)> alpha*beta^j*theta)
        j = j+1;
        U1 = ((1-beta^j)*U + beta^j*U_star);
        X1 = compute_x(U1);
        J1 = compute_cost(X1,U1);  
    end
    last_j=j;
    U=U1;
    iter_num=iter_num +1;
end



end

function f1 = f1(x)
f1 = [ 1-x(1)^0.5 ; x(1)^0.5-x(2)^0.5];
end

function f2 = f2(x)
f2 = [ 2-x(1)^0.5 ; x(1)^0.5-x(2)^0.5];
end

function L1 = L1(x1,u1,r1)
L1=2*(x1(2)-r1)^2;
end

function f1_prime = f1_prime(x1)
f1_prime = [ -0.5*(x1(1))^(-0.5) 0; 0.5*x1(1)^(-0.5) -0.5*x1(2)^(-0.5)];
end

function f2_prime = f2_prime(x1)
f2_prime = [ -0.5*(x1(1))^(-0.5) 0; 0.5*x1(1)^(-0.5) -0.5*x1(2)^(-0.5)];
end

function L1_prime = L1_prime(x1,u1,r1)
L1_prime = [0, 4*(x1(2) -r1)];
end

function X1=compute_x(U1)
global dt;
global x0;

X1=[];
x= x0;
N =length(U1);
for i=1:N
    X1=[X1,x];
    x=x+dt*( (1-U1(i))*f1(x) + U1(i)*f2(x));    
end
end

function P1=compute_p(X1,U1)

global dt;
global pf;
global r;
N=length(U1);
P1=[];
p=pf;

for i=0:N-1
    P1=[p,P1];
    p=p-dt*(-((1-U1(end-i))*f1_prime(X1(:,end-i)) + U1(end-i)*f2_prime(X1(:,end-i)))'*p-(L1_prime(X1(:,end-i),U1(end-i),r(end-i)))');
end
end

function J1 = compute_cost(X1,U1)
global dt;
global r;
J1= 0;
N=length(U1);
for i =1:N
    J1=J1+dt*(L1(X1(:,i),U1(i),r(i)));
end
end

function H1 = compute_h(P1,U1,X1)
global r;
N=length(U1);
H1=[];
for i=1:N
    h=P1(:,i)'*((1-U1(i))*f1(X1(:,i)) + U1(i)*f2(X1(:,i)))+L1(X1(:,i),U1(i),r(i));
    H1=[H1,h];
end
end

function U1_star = compute_u_star(P1,X1)
N=length(X1);
U1_star =[];

% for i=1:N
%     if P1(:,i)'*f1(X1(:,i)) < P1(:,i)'*f2(X1(:,i))
%         u=0;
%     else
%         u=1;
%     end
%     U1_star=[U1_star,u];
% end
U1_star = (P1'*f1(X1) >= P1'*f2(X1))';

end

function theta1 = compute_theta(H1_star,H1)
global dt;
N=length(H1);
theta1=sum((H1_star-H1))*dt;
end