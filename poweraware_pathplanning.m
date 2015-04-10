function [X,U,J_array,U_star,T] = poweraware_pathplanning(c)
global tf;
global X0;
global x0;
global x7;
global pf;
global dt;
global C;

format long g
tf=15;
C = c;
x0=[0;0];
x1=[2;1];
x2=[4;-4];
x3=[7;-5];
x4=[12;7];
x5=[15;0];
x6=[19;3];
x7=[20;0];
X0=[x1;x2;x3;x4;x5;x6];
pf=zeros(12,1);
dt=0.01;
T=0:dt:tf;

max_iter = 500;
iter_num = 0;
N=length(T);
U=zeros(12,N);
a = 0.6;
b = 0.5;
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


end


function L1 = L1(x1,u1)
global x0;
global x7;
global C;
L1 = 0;
for i=1:7
    if i == 1
        L1 = L1 + norm(x1(2*i-1:2*i)-x0)^2 + C*norm(u1(2*i-1:2*i));
    elseif i == 7
        L1 = L1 + norm(x7 - x1(2*(i-1)-1:2*(i-1)))^2; 
    else
        L1 = L1 + norm(x1(2*i-1:2*i)-x1(2*(i-1)-1:2*(i-1)))^2 + C*norm(u1(2*i-1:2*i));
    end
end
end

function X1=compute_x(U1)
global dt;
global X0;

X1=[];
x= X0;
N =length(U1);
for i=1:N
    X1=[X1,x];
    x=x+dt*(U1(:,i)) ;
end
end

function P1=compute_p(X1,U1)

global dt;
global pf;
global x0;
global x7;

N=length(U1);
P1=[];
p=pf;

for i=0:N-1
    P1=[p,P1];
    x = X1(:,end-i);
    x1 = x(1*2-1:1*2);
    x2 = x(2*2-1:2*2);
    x3 = x(3*2-1:3*2);
    x4 = x(4*2-1:4*2);
    x5 = x(5*2-1:5*2);
    x6 = x(6*2-1:6*2);

    p1dot = 2*( x2 - 2*x1 + x0 );
    p2dot = 2*( x3 - 2*x2 + x1 );
    p3dot = 2*( x4 - 2*x3 + x2 );
    p4dot = 2*( x5 - 2*x4 + x3 );
    p5dot = 2*( x6 - 2*x5 + x4 );
    p6dot = 2*( x7 - 2*x6 + x5 );
    pdot = [p1dot;p2dot;p3dot;p4dot;p5dot;p6dot];
    p=p-dt*pdot;
end
end

function J1 = compute_cost(X1,U1)
global dt;
J1= 0;
N=length(U1);
for i =1:N
    J1=J1+dt*(L1(X1(:,i),U1(:,i)));
end
end

function H1 = compute_h(P1,U1,X1)
N=length(U1);
H1=[];
for i=1:N
    h= P1(:,i)'*U1(:,i)+L1(X1(:,i),U1(:,i));
    H1=[H1,h];
end
end

function U1_star = compute_u_star(P1,X1)
global C;
N=length(X1);
U1_star =[];

for i=1:N
    u = [];
    for j=1:6
        p =  P1(2*j-1:2*j,i);
        pnorm = norm(P1(2*j-1:2*j,i));	
        if C - pnorm < 0
            u = [ u ; -p/(pnorm)];
        else
            u = [ u ; 0; 0];
        end
    end
    U1_star=[U1_star,u];
end

end

function theta1 = compute_theta(H1_star,H1)
global dt;
theta1=sum((H1_star-H1))*dt;
end
