clear all
close all

format long



for p=1:4
N=(p)*10; %Number of elements
a=0;
b=1;
h(p)=(b-a)/N;  %Step Size
x=a:h(p):b;

Ke= [1/h(p), -1/h(p) ; -1/h(p), 1/h(p)];  % Element stiffness matrix
Nq=3; %Number of quadrature points
qx=[-sqrt(3/5),0,sqrt(3/5)]; %Quadrature points
qw=[5/9,8/9,5/9]; %Quadrature weigh(p)ts

K=zeros(N+1,N+1);
F=zeros(N+1,1);

for i=1:N
    %Find load vector for each(p) element
   Fe=zeros(2,1);  %Element force
   for j =1:2
       for q=1:Nq
          basis=[(1-qx(q))/2,(1+qx(q))/2];
          pt=(x(i+1)-x(i))/2*qx(q)+(x(i+1)+x(i))/2;
          Fe(j)=Fe(j)+h(p)/2*qw(q)*f(pt)*basis(j);
       end
    
   end  
   %Assembly
   K([i,i+1],[i,i+1])= K([i,i+1],[i,i+1]) +Ke;    
   F([i,i+1],1)=F([i,i+1],1)+Fe;
    
end    

%boundary conditions
boundary_nodes=[1,N+1];
boundary_val=[0,0];
F=F-boundary_val(1)*(K(:,1))-boundary_val(2)*(K(:,end));
%F=F-boundary_val.*K(:,boundary_nodes);

full_nodes=1:N+1;
inner_nodes=setdiff(full_nodes,boundary_nodes);


%Solve for th(p)e system

u=zeros(N+1,1);
u(inner_nodes)=K(inner_nodes,inner_nodes)\F(inner_nodes);

u(1)=boundary_val(1);
u(N+1)=boundary_val(2);

%Visualization
plot(x,u);
hold on


%%%%Compare with exact solution

exact =x.*(1-x).*exp(x);
error(p)=norm(exact'-u,2);
end

%%%Find order of convergence
for j=1:3
order(j)=log(error(j)/error(j+1))/log(h(j)/h(j+1));
end