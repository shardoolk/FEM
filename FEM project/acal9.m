function v10= acal9(x,a1,b1,alpha,beta,gamma,delta,xj,yj,f2)

v10=f2*1/((xj-alpha)*(xj-beta)*(yj-gamma)*(yj-delta))...
    *((x*(b1-a1)+(a1+b1))/2-alpha)*((x*(b1-a1)+(a1+b1))/2-beta) ...
    *(b1-a1)/2;
end