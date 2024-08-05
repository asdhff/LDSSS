function [Z,C,C1,S] = LDSSS_Hessian(X, Par, mu)
%Input is p*n
TanParam.NCoordDim=2;
TanParam.DimGiven=1;

nn = 5;
v = size(X,2);
B = cell(1,v);
for i = 1:v
    B{i} = ConstructHessian(X{i}',GetKNN(X{i}',nn),TanParam);
end

n = size(X{1},2);

muMax=1e8;
iterMax=Par.maxIter;
toler  = 1e-3;
iter = 2;

Z = cell(1, v);
C = cell(1, v);
C1 = cell(1, v);
Y = cell(1, v);
S = zeros(n, n);
M = cell(1, v);

for i = 1:v
    M{i} = X{i}'*X{i};
    Z{i} = zeros(n, n);
    C{i} = zeros(n, n);
    C1{i} = zeros(n, n);
    Y{i} = zeros(n, n);
end

con=zeros(1,iterMax);

C_store = cell(1, iterMax+1);
S_store = cell(1, iterMax+1);
C_store{1}=C{1};
S_store{1}=S;

while  ( iter<(iterMax+1) )
     Z = Z_update(M, C, B, Par);
     C = C_update(Z,C1,Y,Par,mu);
     C1 = C1_update(C,S,Y,Par,mu);
     Y = Y_update(C,C1,Y,mu);
     S = S_update(C1,Par);
     mu = min(mu*Par.rho,muMax);
     C_store{iter}=C{1};
     S_store{iter}=S;
     C_diff = norm(C{1} - C_store{(iter-1)},inf);
     C1_diff = norm(C1{1} - C{1},inf);
     S_diff = norm(S - S_store{(iter-1)},inf);
     con(1, iter-1)= max([C_diff,C1_diff,S_diff]);
     if con(1, iter-1) <= toler
        break
     end
     iter = iter + 1;
end

end
