function v4 = acal4(x,xk,yk,alpha1,alpha2,beta1,beta2,gamma2,delta2,a1,a2,b1,b2)
v4 = (1/((xk- alpha2)*(xk - beta2)*(yk - gamma2)*(yk - delta2)))*( (x*(b1 - a1) + (a1 + b1))/2 - alpha1)*((x*(b1 - a1) + (a1 + b1))/2 - beta1)*((x*(b1 - a1) + (a1 + b1))/2 - alpha2)*((x*(b1 - a1) + (a1 + b1))/2 - beta2)*(b2-a2)/2;
end