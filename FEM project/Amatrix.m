w1 = 0.555555;
w2 = 0.888888;
w3 = 0.555555;
z1 = 0.774596669;
z2 = 0;
z3 = -z1;
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
            
            Aint1 = w1*acal1(z1,xj,yj,alpha1,beta1,alpha2,beta2,gamma1,delta1,a1,b1) + w2*acal1(z2,xj,yj,alpha1,beta1,alpha2,beta2,gamma1,delta1,a1,b1) + w3*acal1(z3,xj,yj,aplha1,beta1,aplha2,beta2,gamma1,delta1,a1,b1);
            
            Aint2 = w1*acal2(z1,xk,yk,alpha2,beta2,gamma1,gamma2,delta1,delta2,a2,b2) + w2*acal2(z2,xk,yk,alpha2,beta2,gamma1,gamma2,delta1,delta2,a2,b2) + w3*acal2(z3,xk,yk,alpha2,beta2,gamma1,gamma2,delta1,delta2,a2,b2);
            
            Aint3 = w1*acal3(z1,xj,yj,alpha1,beta1,gamma1,gamma2,delta1,delta2,a1,a2,b1,b2) + w2*acal3(z2,xj,yj,alpha1,beta1,gamma1,gamma2,delta1,delta2,a1,a2,b1,b2) + w3*acal3(z3,xj,yj,alpha1,beta1,gamma1,gamma2,delta1,delta2,a1,a2,b1,b2);
            
            Aint4 = w1*acal4(z1,xk,yk,alpha1,alpha2,beta1,beta2,gamma2,delta2,a1,a2,b1,b2) + w2*acal4(z2,xk,yk,alpha1,alpha2,beta1,beta2,gamma2,delta2,a1,a2,b1,b2) + w3*acal4(z3,xk,yk,alpha1,alpha2,beta1,beta2,gamma2,delta2,a1,a2,b1,b2);
            
            AintTemp = Aint1*Aint2 + Aint3*Aint4;
            
            Aint = Aint + AintTemp;
        end
        A1(i,j) = Aint;
    end
end

            