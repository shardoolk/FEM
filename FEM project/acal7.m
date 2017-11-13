function v8=acal7(x,b1,a1,xj,yj,alpha,beta,gamma,delta)
v8=1/((xj-alpha)*(xj-beta)*(yj-gamma)*(yj-delta))*((x*(b1-a1)+(a1+b1))/2-alpha)*((x*(b1-a1)+a1+b1)/2-beta)*(b1-a1)/2;
end
