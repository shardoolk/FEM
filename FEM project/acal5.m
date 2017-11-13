function v6=acal5(x,xj,yj,alpha,beta,gamma,delta,a1,b1)

v6=(1/((xj-alpha)*(xj-beta)*(yj-gamma)*(yj-delta)))*(x*(b1-a1)+a1+b1-(alpha+beta))*(b1-a1)/2;
end
