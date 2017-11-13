
close all
global nx ny N hx hy nnodes A1 B1 F1 Nc C Afinal D Sorted vel;

A = Afinal;%xlsread('Amatrix.xlsx');
BT = B1;%xlsread('BTmatrix.xlsx');
%D = xlsread('Dmatrix.xlsx');
F = F1;%xlsread('Ffmatrix.xlsx');


ep = 0.01;
B = BT';
K = A + (B*BT)/ep;
L = F' + (B*D)/ep;
%Replace by SOR function
X = matrix_solver(K,L);
X = X/abs(vel);
Theta = (D - BT*X)/ep;

 iter = 1;
for i = 1:nnodes
    k = 0;
    for j = 1:length(Sorted) 
        if Sorted(j) == i;
            k = 1;
            break
        else
        end
    end
    if  k ==0;
        P(iter,1) = Nc(i,1);
        P(iter,2) = Nc(i,2);
        P(iter,3) = X(2*iter -1);
        P(iter,4) = X(2*iter);
        iter = iter + 1;
    end
end

for i = 1:N
    Xp(i) = rem(i,nx);
    if Xp(i) == 0;
        Xp(i) = nx;
    else
    end
    Yp(i) = fix((i-1)/nx);
end
Pfinal = zeros(nx,ny);
iter =1;
for i=1:nx
    for j=1:ny
        
        Pfinal(i,j) = Theta(iter,1);
        iter=iter+1;
    end
end
xplot = linspace(0,1,nx);
yplot = linspace(0,1,ny);

[x,y] = meshgrid(xplot,yplot);
%Pressure plot
figure(1);
surf(x,y,Pfinal);
figure(2);
%Velocity plot
quiver(P(:,1),P(:,2),P(:,3),P(:,4));
