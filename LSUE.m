function [value] = LSUE( tr_cell,tr_TF, B, Lamda, te_cell, te_TF )



[~,~,C_ID] = unique(tr_cell);
clear cell
[~,~,T_ID] = unique(tr_TF);
clear TF

c_num = max(C_ID);

A = zeros(length(cell),max(C_ID)+max(T_ID));

for i = 1:size(A,1)
    A(i, C_ID(i)) = 1;
    A(i, c_num+T_ID(i)) = 1;
end

l = size(B,2);

Num = 4 ;

Interval = ceil(l/Num);

Z_Index = cell(1,1);
Z = cell(1,1);

[U1,D1,V1] = svd(A,'econ');
D1 = diag(D1);

U1 = U1(:,D1>max(D1)*1e-6);
V1 = V1(:,D1>max(D1)*1e-6);
D1 = D1(D1>max(D1)*1e-6);


B = U1'*B;

results = [];
T = B*B';
[P,D] = eig(T);

D = diag(D);

index = D>mean(D)*1e-6;

P = P(:,index);

B = P'*B;

W_L = A1*V1*diag(1./D1)*P;



for i = 1:3
    Z_Index{i} = [(i-1)*Interval+1:min((i+1)*Interval,l)];
    Z{i} = zeros(size(B,1),length(Z_Index{(i+1)/2}));
end

W_A = zeros(1,l);
for i = 1:length(Z_Index)
    W_A(Z_Index{i}) = W_A(Z_Index{i}) + 1;
end


Y = Z;

mu = 2e-2;
rol = 1.25;

iter = 30;

for i = 1:iter
    tic
    T = repmat(1./(1 + mu*W_A/2),size(B,1),1);
    X = 2*B;
    for j = 1:length(Z_Index)
        X(:,Z_Index{j}) = X(:,Z_Index{j}) - Y{j} + mu*Z{j}; 
    end
    X = T.*X/2;
    
    rank = 0;
    
    for j = 1:length(Z_Index)
        T = X(:,Z_Index{j}) + Y{j}/mu;
%         size(T)
        T1 = T*T';
        [U,D] = eig((T1+T1')/2);
        D = diag(D);
        D(D<0) = 0;
        D = D.^0.5;
        D = D - Lamda/mu;
        D(D<0) = 0;
        U = U(:,D>0);
        D = D(D>0);
        rank = rank + length(D);
        V = T'*U;
        V = normc(V);
        Z{j} = U*diag(D)*V';
    end
    
    clear T V U
    
    for j = 1:length(Z_Index)
        Y{j} = Y{j} + mu*(X(:,Z_Index{j})-Z{j});
    end
        
    mu = mu*rol;
    
    value = W_L*X; 
end

