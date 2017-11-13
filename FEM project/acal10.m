function v11 =acal10(y,b2,a2,gamma,delta)
v11=((y*(b2-a2)+(a2+b2))/2-gamma)*((y*(b2-a2)+(a2+b2))/2-delta)...
    *(b2-a2)/2;
end