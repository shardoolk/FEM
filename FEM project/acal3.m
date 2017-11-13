function v3 = acal3(y,xj,yj,alpha1,beta1,gamma1,gamma2,delta1,delta2,a1,a2,b1,b2)
v3 = (1/((xj - alpha1)*(xj - beta1)*(yj - gamma1)*(yj - delta1)))*(y*(b2 - a2) + (a2 + b2) - (gamma1 + delta1))*(y*(b2-a2) + (a2 + b2) - (gamma2 + delta2))*(b1-a1)/2;
end