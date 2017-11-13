function b = SourceAssembler1D(x,f,kappa,g)
% n = length(x) - 1;
b = LoadAssembler1D(x,f);
% for i = 1:n
%     h = x(i+1) - x(i);
%     b(i) = b(i) + f(x(i))*h/2;
%     b(i+1) = b(i+1) + f(x(i+1))*h/2;
% end
b(1) = b(1) + kappa(1)*g(1);
b(end) = b(end) + kappa(2)*g(2);
