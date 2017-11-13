function v2 = acal2(y,xk,yk,alpha2,beta2,gamma1,gamma2,delta1,delta2,a2,b2)
v2 = (1/((xk - alpha2)*(xk - beta2)*(yk - gamma2)*(yk - delta2)))*((y*(b2-a2) + (a2 + b2))/2 - gamma1)*((y*(b2-a2) + (a2 + b2))/2 - delta1)*((y*(b2-a2) + (a2 + b2))/2 - gamma2)*((y*(b2-a2) + (a2 + b2))/2 - delta2)*(b2 - a2)/2;
end