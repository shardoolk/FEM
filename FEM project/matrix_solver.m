function theta = matrix_solver(A,F)

%eps = %Give the value of h squared

Afinal = A;

%Initial Solution guess
init_guess = ones(length(Afinal(:,1)),1);

%Under relaxation factor
omega = 1.6;

new_Guess = zeros(length(Afinal(:,1)),1);



    new_Guess(1,1) = (1/Afinal(1,1))*(F(1,1) - Afinal(1,2:length(Afinal(:,1)))*init_guess(2:length(init_guess),1));
for i=2:length(Afinal(:,1))
    new_Guess(i,1) = (1/Afinal(i,i))*(F(i,1)  - Afinal(i,1:i-1)*new_Guess(1:i-1) - Afinal(i,i+1:length(Afinal(:,1)))*init_guess(i+1:length(init_guess),1));
    
end

while abs(sum(new_Guess - init_guess))>1e-6
    init_guess = new_Guess;
    
    new_Guess(1,1) = (1-omega)*init_guess(1,1) + omega*(1/Afinal(1,1))*(F(1,1) - Afinal(1,2:length(Afinal(:,1)))*init_guess(2:length(init_guess),1));
    for i=2:length(Afinal(:,1))
    
           new_Guess(i,1) = (1-omega)*init_guess(i,1) + omega*(1/Afinal(i,i))*(F(i,1) - Afinal(i,i+1:length(Afinal(:,1)))*init_guess(i+1:length(init_guess),1) - ...
               Afinal(i,1:i-1)*new_Guess(1:i-1));
    end
end

theta=new_Guess;



end
