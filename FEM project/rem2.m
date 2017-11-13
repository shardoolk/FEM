%node to coordinate 
function v = rem2(x,nx)
v = rem(x,2*nx+1);
if v == 0;
    v = 2*nx + 1;
else
end

