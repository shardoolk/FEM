
%Dirichelt Nodes
global a b c d nx ny N hx hy nnodes A1 B1 F1 Nc C Afinal vel;
vel= -0.0001;

DiriN = zeros(4*(nx+ny),1);
for i = 1:(2*nx+1)
    DiriN(i,1) = i;
end
for j = 1:(2*ny)
    DiriN(2*nx+1+j,1) = DiriN((2*nx + j),1) + (2*nx+1);
end
for k = 1:2*nx
    DiriN(2*(nx +ny) + 1 +k,1 ) = DiriN(2*(nx + ny) + k,1) - 1;
end
for m = 1:(2*ny-1)
    DiriN(4*nx + 2*ny + 1 + m,1) = DiriN(4*nx + 2*ny + m,1) - (2*nx+1);
end

Afinal=zeros(2*nnodes,2*nnodes);
p1=1;
k1=1;
for i=1:length(A1)
    for j=1:length(A1)
        Afinal(p1,k1)=A1(i,j);
        Afinal(p1+1,k1+1)=A1(i,j);
        k1=k1+2;
    end
    k1=1;
    p1=p1+2;
end

 Sorted = sort(DiriN,'descend');
% 
 F2=zeros(2*nnodes,2*length(Sorted));
 for i = 1:length(Sorted)
     F2(:,2*i) = Afinal(:,2*Sorted(i));
     F2(:,(2*i-1)) = Afinal(:,2*Sorted(i)-1);
     Afinal(:,2*Sorted(i)) = [];
     Afinal(:,2*Sorted(i)-1)=[];
%      Afinal(2*Sorted(i),:)=[];
%      Afinal(2*Sorted(i)-1,:)=[];
 end
 
 for i = 1:length(Sorted)
     Afinal(2*Sorted(i),:)=[];
     Afinal(2*Sorted(i)-1,:)=[];
     F2(2*Sorted(i),:) = [];
     F2(2*Sorted(i)-1,:) = [];
 end

  
% 
 D1 = zeros(N,2*length(DiriN));
% B1 = xlsread('B1matrix.xlsx');
% 
 for i = 1:length(Sorted)
     D1(:,2*i) = B1(:,2*Sorted(i));
     D1(:,2*i-1) = B1(:,2*Sorted(i)-1);
     B1(:,2*Sorted(i)) = [];
     B1(:,2*Sorted(i)-1)=[];
 end
 %xlswrite('BTmatrix.xlsx',B1);
 %xlsread('F1matrix.xlsx');
   for i = 1:length(Sorted)
   F1(2*Sorted(i)) = [];
   F1(2*Sorted(i)-1)=[];
   end
 
 
zz = 2*nx + 1;
  for i = 1:zz
      for j = 1:length(F1)
      F1(j)=F1(j)-vel*F2(j,2*i-1);
      end
  end
  
  %xlswrite('Ffmatrix.xlsx',F1);
  D=zeros(N,1);
   for i = 1:zz
      for j = 1:length(D)
      D(j)=D(j)-vel*D1(j,2*i-1);
      end
   end
