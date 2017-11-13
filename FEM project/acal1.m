function v1 = acal1(x,xj,yj,alpha1,beta1,alpha2,beta2,gamma1,delta1,a1,b1)
v1 = (1/((xj - alpha1)*(xj - beta1)*(yj - gamma1)*(yj - delta1)))*(x*(b1-a1) + (a1+b1) - (alpha1 + beta1))*(x*(b1-a1) + (a1 + b1) -(alpha2 + beta2))*(b1-a1)/2;
end
