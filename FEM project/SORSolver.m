%Code for SOR
close all;
clear all;
omega  = 1.8;
Err = 1;
tol = 10^-8;
N = 10;
A = rand(N);
B = zeros(N,1);
for i = 1:N
    B(i) = i*0.1;
end
L = zeros(N);
U = zeros(N);
D = zeros(N);
for i = 1:N
    for j = 1:N
        if i > j
          L(i,j)=  A(i,j)  ;
        elseif j > i
            U(i,j) = A(i,j);
        else 
            D(i,j) = A(i,j) ;
        end
    end
end
E = (L+U+D) - A;
u = zeros(N,1);
u1 = u;

for i = 1:N
    u1(i) = -inv((D + omega*L))*((omega*U + (1 - omega)*D)*u(i)) + inv((D + omega*L))*(omega*B);
end
