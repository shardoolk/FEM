clear all

global a b c d nx ny N hx hy nnodes A1 B1 F1 Nc C A BT D F Sorted;


% global nx ny C integration N nnodes Nc localnumber1 localnumber2

%a,b,c,d are the coordinate combinations
a = 0;
b = 1;
c = 0;              
d = 1;
%nx - number of elements in x direction
nx = 10;
%ny - number of elements in y direction
ny = 10;
hx = (b-a)/nx;
hy = (d-c)/ny;
N = nx*ny;
%C - connectivity matrix for velocity elements(9-nodes)
C = zeros(N,9);

%Assembling the connectivity matrix
for i = 1:N
    n = rem(i,nx);
    m = fix((i-1)/nx) + 1;
    if n == 0
        n = nx;
    else
    end
        C(i,1) = (2*n-1) + (2*m-2)*(2*nx + 1);
        C(i,2) = C(i,1) + 1;
        C(i,3) = C(i,2) + 1;
        C(i,4) = C(i,1) + (2*nx + 1);
        C(i,5) = C(i,4) + 1;
        C(i,6) = C(i,5) + 1;
        C(i,7) = C(i,4) + (2*nx +1);
        C(i,8) = C(i,7) + 1;
        C(i,9) = C(i,8) + 1;
    
end
%Node to Coordinate matrix
nnodes = (2*nx +1)*(2*ny+1);%Number of nodes in the velocity space
Nc = zeros(nnodes,2);
for i = 1:nnodes
    Nc(i,1) = (rem2((i),nx)-1)*hx/2;
    Nc(i,2) = fix((i-1)/(2*nx+1))*hy/2;
end
% Elements of Intergration - structures for storing the elements of
% integration
integration  = zeros(nnodes,nnodes,4);
localnumber1 = zeros(nnodes,nnodes,4);
localnumber2 = zeros(nnodes,nnodes,4);
for i = 1:nnodes
    for j = 1:nnodes
        l = 1;
        for k = 1:(nx*ny)
            p = 0;
            q = 0;
            for m = 1:9
                if C(k,m) == i;
                    p = m;
                else
                end
                if C(k,m) == j;
                    q = m;
                else
                end
            end
            if (p ~= 0 && q ~= 0)
                integration(i,j,l) = k;
                localnumber1(i,j,l) = p;
                localnumber2(i,j,l) = q;
                l = l + 1;
            end
        end
    end
end


%Gaussian quadrature thrid order weights
w1 = 0.555555;
w2 = 0.888888;
w3 = 0.555555;
%Gaussian quad points
z1 = 0.774596669;
z2 = 0;
z3 = -z1;

%Gaussian integration for each of Aij matrix entry
for i = 1:nnodes
    for j = 1:nnodes
        Aint = 0;
        AintTemp = 0;
        for p = 1:4
            element = integration(i,j,p);
            if element == 0
                break
            end
            a1 = Nc(C(element,1),1);
            a2 = Nc(C(element,1),2);
            b1 = Nc(C(element,9),1);
            b2 = Nc(C(element,9),2);
            localn1 = localnumber1(i,j,p);
            localn2 = localnumber2(i,j,p);
            xj = Nc(C(element,localn1),1);
            yj = Nc(C(element,localn1),2);
            xk = Nc(C(element,localn2),1);
            yk = Nc(C(element,localn2),2);
            alpha1 = Nc(C(element,rem1(localn1 + 1)),1);
            beta1 = Nc(C(element,rem1(localn1 + 2)),1);
            gamma1 = Nc(C(element,rem1(localn1 + 3)),2);
            delta1 = Nc(C(element,rem1(localn1 + 6)),2);
            
            alpha2 = Nc(C(element,rem1(localn2 + 1)),1);
            beta2 = Nc(C(element,rem1(localn2 + 2)),1);
            gamma2 = Nc(C(element,rem1(localn2 + 3)),2);
            delta2 = Nc(C(element,rem1(localn2 + 6)),2);
            
            Aint1 = w1*acal1(z1,xj,yj,alpha1,beta1,alpha2,beta2,gamma1,delta1,a1,b1) + w2*acal1(z2,xj,yj,alpha1,beta1,alpha2,beta2,gamma1,delta1,a1,b1) + w3*acal1(z3,xj,yj,alpha1,beta1,alpha2,beta2,gamma1,delta1,a1,b1);
            
            Aint2 = w1*acal2(z1,xk,yk,alpha2,beta2,gamma1,gamma2,delta1,delta2,a2,b2) + w2*acal2(z2,xk,yk,alpha2,beta2,gamma1,gamma2,delta1,delta2,a2,b2) + w3*acal2(z3,xk,yk,alpha2,beta2,gamma1,gamma2,delta1,delta2,a2,b2);
            
            Aint3 = w1*acal3(z1,xj,yj,alpha1,beta1,gamma1,gamma2,delta1,delta2,a1,a2,b1,b2) + w2*acal3(z2,xj,yj,alpha1,beta1,gamma1,gamma2,delta1,delta2,a1,a2,b1,b2) + w3*acal3(z3,xj,yj,alpha1,beta1,gamma1,gamma2,delta1,delta2,a1,a2,b1,b2);
            
            Aint4 = w1*acal4(z1,xk,yk,alpha1,alpha2,beta1,beta2,gamma2,delta2,a1,a2,b1,b2) + w2*acal4(z2,xk,yk,alpha1,alpha2,beta1,beta2,gamma2,delta2,a1,a2,b1,b2) + w3*acal4(z3,xk,yk,alpha1,alpha2,beta1,beta2,gamma2,delta2,a1,a2,b1,b2);
            
            AintTemp = Aint1*Aint2 + Aint3*Aint4;
            
            Aint = Aint + AintTemp;
        end
        A1(i,j) = Aint;
    end
end


%xlswrite('A1matrix.xlsx',A1);          


%Data structures for the B matrix
integrationB  = zeros(nnodes,4);
localnumberB = zeros(nnodes,4);

%Data structure assembly
for i = 1:nnodes
   
        l = 1;
        for k = 1:(nx*ny)
            p = 0;
            q = 0;
            for m = 1:9
                if C(k,m) == i;
                    p = m;
                else
                end
                
            end
            if (p ~= 0 )
                integrationB(i,l) = k;
                localnumberB(i,l) = p;
                
                l = l + 1;
            end
        end
  
end

%Bij matrix entry gaussian quadrature
for j = 1:N
    for i = 1:nnodes
        AintB1=0;
        AintB2=0;
        for p = 1:4
            element = integrationB(i,p);
            if element == 0
                break
            end
            if element == j
            a1 = Nc(C(element,1),1);
            a2 = Nc(C(element,1),2);
            b1 = Nc(C(element,9),1);
            b2 = Nc(C(element,9),2);
            localn1 = localnumberB(i,p);
            xj = Nc(C(element,localn1),1);
            yj = Nc(C(element,localn1),2);
            alpha = Nc(C(element,rem1(localn1 + 1)),1);
            beta = Nc(C(element,rem1(localn1 + 2)),1);
            gamma = Nc(C(element,rem1(localn1 + 3)),2);
            delta = Nc(C(element,rem1(localn1 + 6)),2);
            
            Aint5 = w1*acal5(z1,xj,yj,alpha,beta,gamma,delta,a1,b1) + w2*acal5(z2,xj,yj,alpha,beta,gamma,delta,a1,b1) + w3*acal5(z3,xj,yj,alpha,beta,gamma,delta,a1,b1);
            
            Aint6 = w1*acal6(z1,a2,b2,gamma,delta) + w2*acal6(z2,a2,b2,gamma,delta) + w3*acal6(z3,a2,b2,gamma,delta);
            
            Aint7= w1*acal7(z1,b1,a1,xj,yj,alpha,beta,gamma,delta) + w2*acal7(z2,b1,a1,xj,yj,alpha,beta,gamma,delta) + w3*acal7(z3,b1,a1,xj,yj,alpha,beta,gamma,delta);
            
            Aint8= w1*acal8(z1,a2,b2,gamma,delta) + w2*acal8(z2,a2,b2,gamma,delta) + w3*acal8(z3,a2,b2,gamma,delta);
            AintB1 = Aint5*Aint6;
            AintB2 = Aint7*Aint8;
            
            break
            end
      end
        
      B1(j,2*i-1) = AintB1;
      B1(j,2*i) = AintB2;
    end
end


%xlswrite('B1matrix.xlsx',B1);


%Forcing function calculation
f1=0;
f2=0; %Negative of Rho*g

    for i = 1:nnodes
        AintF1=0;
        AintF1temp=0;
        for p = 1:4
            element = integrationB(i,p);
            if element == 0
                break
            end
   
            a1 = Nc(C(element,1),1);
            a2 = Nc(C(element,1),2);
            b1 = Nc(C(element,9),1);
            b2 = Nc(C(element,9),2);
            localn1 = localnumberB(i,p);
            xj = Nc(C(element,localn1),1);
            yj = Nc(C(element,localn1),2);
            alpha = Nc(C(element,rem1(localn1 + 1)),1);
            beta = Nc(C(element,rem1(localn1 + 2)),1);
            gamma = Nc(C(element,rem1(localn1 + 3)),2);
            delta = Nc(C(element,rem1(localn1 + 6)),2);
            
            Aint9 = w1*acal9(z1,a1,b1,alpha,beta,gamma,delta,xj,yj,f2) + w2*acal9(z2,a1,b1,alpha,beta,gamma,delta,xj,yj,f2) + w3*acal9(z3,a1,b1,alpha,beta,gamma,delta,xj,yj,f2);
            
            Aint10 = w1*acal10(z1,b2,a2,gamma,delta) + w2*acal10(z2,b2,a2,gamma,delta) + w3*acal10(z3,b2,a2,gamma,delta);
            
         
            AintF1temp = Aint9*Aint10;
            
            AintF1=AintF1+AintF1temp;
       
            
           
      end
        
      F1(2*i-1) = 0;
      F1(2*i) = AintF1;
    end


    
%xlswrite('F1matrix.xlsx',F1);
%xlswrite('Ncmatrix.xlsx',Nc);               
%xlswrite('ConnectivityMatrix.xlsx',C);