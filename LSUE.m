function [value] = LSUE( tr_Cell,tr_TF, O, Group_Index, Lamda, te_Cell, te_TF )

% Input format: O is n_o*n_p signal matrix with each of its rows corresponding to a chip-seq
% dataset. tr_Cell and tr_TF are two n_p*1 cell arrays of strings, which records the name of TF and Cell for the
% corresponding row of O. Group_Index is a n_g*1 cell array, with its
% i-th element recording the index of genomic pisitions arranged into the
% i-th group. te_Cell and te_TF are two strings which records the name of
% TF and Cell for which we need to impute the corresponding profile.
%Output format: value is a 1*n_p matrix, which returns the predicted
%profile of te_TF in te_Cell

[Cell,~,C_ID] = unique(tr_Cell);
clear cell
[TF,~,T_ID] = unique(tr_TF);
clear TF

ng = length(Group_Index);
np = size(O,2);

nc = max(C_ID);
nm = max(T_ID);

A = zeros(length(cell),max(C_ID)+max(T_ID));

for i = 1:size(A,1)
    A(i, C_ID(i)) = 1;
    A(i, nc+T_ID(i)) = 1;
end

tmp = strcmp(te_Cell,Cell);
index_c = tmp>0;

tmp = strcmp(te_TF,TF);
index_t = find(tmp>0);

te_A = zeros(1,size(A,2));
te_A(index_c) = 1;
te_A(index_t+nc) = 1;

l = size(O,2);


[U1,D1,V1] = svd(A,'econ');
D1 = diag(D1);

U1 = U1(:,D1>max(D1)*1e-6);
V1 = V1(:,D1>max(D1)*1e-6);
D1 = D1(D1>max(D1)*1e-6);


O = U1'*O;

T = O*O';
[P,D] = eig(T);

D = diag(D);

index = D>mean(D)*1e-6;
P = P(:,index);
O = P'*O;
W_L = te_A*V1*diag(1./D1)*P;

Z = cell(1,ng);

for i = 1:ng
    Z{i} = zeros(size(O,1),length(Group_Index{i}));
end

W_A = zeros(1,l);
for i = 1:length(Group_Index)
    W_A(Group_Index{i}) = W_A(Group_Index{i}) + 1;
end


Y = Z;

mu = 1e-2;
rol = 1.1;

iter = 30;
% Compute X using ADMM
for i = 1:iter
    T = repmat(1./(1 + mu*W_A/2),size(O,1),1);
    X = 2*O;
    for j = 1:length(Group_Index)
        X(:,Group_Index{j}) = X(:,Group_Index{j}) - Y{j} + mu*Z{j}; 
    end
    X = T.*X/2;
    
    for j = 1:length(Group_Index)
        T = X(:,Group_Index{j}) + Y{j}/mu;
        T1 = T*T';
        [U,D] = eig((T1+T1')/2);
        D = diag(D);
        D(D<0) = 0;
        D = D.^0.5;
        D = D - Lamda/mu;
        D(D<0) = 0;
        U = U(:,D>0);
        D = D(D>0);
        V = T'*U;
        V = normc(V);
        Z{j} = U*diag(D)*V';
    end
    
    clear T V U
    
    for j = 1:length(Group_Index)
        Y{j} = Y{j} + mu*(X(:,Group_Index{j})-Z{j});
    end        
    mu = mu*rol;
end

% Compute B using ADMM
mu = 1e-2;
rol = 1.1;
iter = 30;

Oe = mean(O,2);

B = zeros(nm,nc);
B1 = zeros(nm,nc);
Z = zeros(nm,nc);

OB = B;

Omega = zeros(nm,nc);

for i = 1:size(A,1)
    Omega(T_ID(i),C_ID(i)) = 1;
    OB(T_ID(i),C_ID(i)) = Oe(i);
end

for i = 1:iter
    B(Omega==0) = B1(Omega==0) - Z(Omega==0)/mu;
    B(Omega==1) = (B(Omega==1)*mu/2 + OB(Omega==1) - Z(Omega==1)/2)/(1+mu/2);
    
    B1 = B + Z/mu;
    
    [U,D,V] = svd(B1,'econ');
    D = diag(D);
    D = D - Lamda/mu;
    D(D<0) = 0;
    B1 = U*diag(D)*V';
    
    Z = Z + mu*(B-B1);
    
    mu = mu*rol;
end

value = W_L*X + B(index_t,index_c); 


