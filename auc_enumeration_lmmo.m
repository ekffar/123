function [PWM, PFM, sSn, nCC] = auc_enumeration_range_tr_artificial_reco(  P_Train_Overall_Data_Numerical_Sub, P_Train_Overall_Data_Numerical_Sub_Index, N_Train_Overall_Data_Numerical_Sub, N_Train_Overall_Data_Numerical_Sub_Index,...
    P_Tr_Seq_length, N_Tr_Seq_length, w0, ratio, true_tfbs_starting_position)

if ~exist('w0', 'var') || isempty(w0)
    w0 = randn(size(P_Train_Overall_Data_Numerical_Sub,2),1);
end

P_Train_Overall_Data_Numerical_Sub = double(P_Train_Overall_Data_Numerical_Sub);
N_Train_Overall_Data_Numerical_Sub = double(N_Train_Overall_Data_Numerical_Sub);

Logical_P_Train_Overall_Data_Numerical_Sub = logical(P_Train_Overall_Data_Numerical_Sub);
Logical_N_Train_Overall_Data_Numerical_Sub = logical(N_Train_Overall_Data_Numerical_Sub);

w0 = double(w0);

no_dims = length(w0);

np = max(P_Train_Overall_Data_Numerical_Sub_Index);
nn = max(N_Train_Overall_Data_Numerical_Sub_Index);

PWM = cell(1,20);
PFM = cell(1,20);

sSn = [];
nCC = [];

P_Score = zeros(np,1);
N_Score = zeros(nn,1);

[np nn]

P_Indicator = zeros(size(P_Train_Overall_Data_Numerical_Sub,1),1);
N_Indicator = zeros(size(N_Train_Overall_Data_Numerical_Sub,1),1);

for j = 1:np
    Index = [P_Tr_Seq_length(j)+1:P_Tr_Seq_length(j+1)];
    P_Indicator(Index) = j;
end
for j = 1:nn
    Index = [N_Tr_Seq_length(j)+1:N_Tr_Seq_length(j+1)];
    N_Indicator(Index) = j;
end

P_Indicator_1 = cell(1,length(w0));
P_Indicator_2 = cell(1,length(w0));
N_Indicator_1 = cell(1,length(w0));
N_Indicator_2 = cell(1,length(w0));

 for i = 1:length(w0)
     P_Indicator_1{i} = find(P_Train_Overall_Data_Numerical_Sub(:,i)==0)-1;
     P_Indicator_2{i} = find(P_Train_Overall_Data_Numerical_Sub(:,i)==1)-1;
     
     N_Indicator_1{i} = find(N_Train_Overall_Data_Numerical_Sub(:,i)==0)-1;
     N_Indicator_2{i} = find(N_Train_Overall_Data_Numerical_Sub(:,i)==1)-1;
 end



for iter = 1:1
    %     iter
    %     for i = 1:length(w0)/4
    %         w0(1+(i-1)*4:i*4) = w0(1+(i-1)*4:i*4) - mean(w0(1+(i-1)*4:i*4));
    %     end
    iter
    position_auc = zeros(size(w0));
    position_w = position_auc;
    for i = 1:1
        
        
        w = w0;
        
        if i == 1
            P_Value = (P_Train_Overall_Data_Numerical_Sub)*w;
        end
       
        P_Dual_Score= ComputeDual_Score_1(P_Value,P_Indicator, P_Indicator_1{i}, P_Indicator_2{i}, np, length(P_Indicator_1{i}), length(P_Indicator_2{i}));
        
        P_Dual_Score = P_Dual_Score';
        
        P_Coord = P_Dual_Score;
        
        P_Coord(:,1) = P_Dual_Score(:,2) - P_Dual_Score(:,1);
        
        if i == 1
            N_Value = N_Train_Overall_Data_Numerical_Sub*w;
        end
        
        N_Dual_Score= ComputeDual_Score_1(N_Value,N_Indicator,  N_Indicator_1{i}, N_Indicator_2{i}, nn, length(N_Indicator_1{i}), length(N_Indicator_2{i}));
        N_Dual_Score = N_Dual_Score';
        
        N_Coord = N_Dual_Score;
             
        N_Coord(:,1) = N_Dual_Score(:,2) - N_Dual_Score(:,1);
        
        X = [P_Coord(:,1);N_Coord(:,1)];
        [~,min_index] = sort(abs(X));
        
        tl = min(X(min_index(1:ceil(length(min_index)*ratio))));
        tr = max(X(min_index(1:ceil(length(min_index)*ratio))));
        
        P_Coord(P_Coord(:,1)<tl,2) = P_Coord(P_Coord(:,1)<tl,2) + tl - P_Coord(P_Coord(:,1)<tl,1);
        P_Coord(P_Coord(:,1)<tl,1) = tl;
        
        N_Coord(N_Coord(:,1)<tl,2) = N_Coord(N_Coord(:,1)<tl,2) + tl - N_Coord(N_Coord(:,1)<tl,1);
        N_Coord(N_Coord(:,1)<tl,1) = tl;
        
        a = tl - 0*rand(length(find(P_Coord(:,1)<tl)),1)*10;       
        
        P_Coord(P_Coord(:,1)<tl,2) = P_Coord(P_Coord(:,1)<tl,2) + a - P_Coord(P_Coord(:,1)<tl,1);
        P_Coord(P_Coord(:,1)<tl,1) = a;
        
        a = tl - 0*rand(length(find(N_Coord(:,1)<tl)),1)*10;
        
        N_Coord(N_Coord(:,1)<tl,2) = N_Coord(N_Coord(:,1)<tl,2) + a - N_Coord(N_Coord(:,1)<tl,1);
        N_Coord(N_Coord(:,1)<tl,1) = a;
        
        P_Coord(P_Coord(:,1)>tr,1) = tr+0*rand(size(P_Coord(P_Coord(:,1)>tr,1),1),1)*10;
        N_Coord(N_Coord(:,1)>tr,1) = tr+0*rand(size(N_Coord(N_Coord(:,1)>tr,1),1),1)*10;
        
        Coord = [N_Coord ones(size(N_Coord,1),1);P_Coord zeros(size(P_Coord,1),1)];
        [Output] = RangeTree(Coord);
        x0 = Output(:,1);
        y = -Output(:,2);
        Score = [0;y];
        x = Score;
        
        a = rand(length(x0)-1,1);
        
        x(2:end-1) = x0(1:end-1).*a + x0(2:end).*(1-a);
        
        x(1) = 2*x(2) - x(3);
        
        x(end) = 2*x(end-1) - x(end-2);
        x1 = x;
        
        [~, min_index] = min(abs(x1));        
        
        Score = Score - Score(min_index);
        
        ll = ceil(length(Score)/2)';
        
        order = randperm(length(Score));
        order = order(1:ll);
        
        order = sort(order);        
        
        Score = Score(order);
        x = x(order);
                
        [position_auc(i), k] = max(Score + 1e-20*rand(size(Score)));
        position_auc(i) = position_auc(i)/np/nn;
        position_w(i) = x(k);
        
    end
%     position_auc
    [~, index] = max(position_auc);
    
    if iter >=2
        w0(index) = w0(index) + position_w(index);
    end
    
    P_Value(Logical_P_Train_Overall_Data_Numerical_Sub(:,index)) = P_Value(Logical_P_Train_Overall_Data_Numerical_Sub(:,index)) + position_w(index);    
    N_Value(Logical_N_Train_Overall_Data_Numerical_Sub(:,index)) = N_Value(Logical_N_Train_Overall_Data_Numerical_Sub(:,index)) + position_w(index);
    
    a = reshape(w0,4,[]);
    for j = 1:size(a,2)
        a(:,j) = a(:,j) - max(a(:,j));
        a(:,j) = exp(a(:,j));
        a(:,j) = a(:,j)/sum(a(:,j));
    end
        
%     a([1 4 2 3],:) = a;
    PWM{iter} = a;
    
    a
    
    P_Instance = zeros(np,length(w0));
    Value = P_Train_Overall_Data_Numerical_Sub*w0;
    for j = 1:np
        index = find(P_Train_Overall_Data_Numerical_Sub_Index == j);
        [P_Score(j), index1] = max(Value(index));
        P_Instance(j,:) = P_Train_Overall_Data_Numerical_Sub(index(index1),:);
    end
    P_Instance = mean(P_Instance);
    %     w0 = log(P_Instance'+eps);
    a = reshape(P_Instance,4,[]);
    a([1 4 2 3],:) = a;
    PFM{iter} = a;
     
    P_Value = P_Train_Overall_Data_Numerical_Sub*w0;
    
    sTP = 0;
    sFP = 0;
    sFN = 0;
    
    nTP = 0;
    nFP = 0;
    nTN = 0;
    nFN = 0;
    
    P_Value = P_Train_Overall_Data_Numerical_Sub*w0;
    
    for j = 1:np
        Value = P_Value(P_Tr_Seq_length(j)+1:P_Tr_Seq_length(j+1));
        P_Score(j) = max(Value);
    end
    
    N_Value = N_Train_Overall_Data_Numerical_Sub*w0;
    
    for j = 1:nn
        Value = N_Value(N_Tr_Seq_length(j)+1:N_Tr_Seq_length(j+1));
        N_Score(j) = max(Value);
    end
    
%     [f1,xi1] = ksdensity(P_Score);
%     [f2,xi2] = ksdensity(N_Score);
% 
%     subplot(2,1,1)
%     plot(xi1,f1,xi2,f2);
    
    Diff = repmat(P_Score,1,length(N_Score)) - repmat(N_Score',length(P_Score),1);
    
    Diff(Diff>0) = 1;
    Diff(Diff==0) = 0.5;
    Diff(Diff<0) = 0;
    
    auc0 = mean(Diff(:));
    
    auc0
    
    Train_AUC(iter) = auc;
    
    if iter >= 5
        if (Train_AUC(iter) - Train_AUC(iter-1)) < 1e-4
            break;
        end
    end
    
    auc
    
    for j = 1:np
        if true_tfbs_starting_position(j,2) >=0
            l = P_Tr_Seq_length(j+1) - P_Tr_Seq_length(j);
            %         l
            predicted_tfbs = zeros(l/2+length(w0)/4-1,1);
            true_tfbs = zeros(l/2+length(w0)/4-1,1);
            true_tfbs(true_tfbs_starting_position(j,1):true_tfbs_starting_position(j,1)+length(w0)/4-1) = 1;
            Value = P_Value(P_Tr_Seq_length(j)+1:P_Tr_Seq_length(j+1));
            [P_Score(j), index] = max(Value+randn(size(Value))*1e-30);
            %         index
            %         l
            if index <= l/2
                predicted_tfbs(index:index+length(w0)/4-1) = 1;
                
            else
                index = index - l/2;
                %             index:index+length(w0)-1
                predicted_tfbs(index:index+length(w0)/4-1) = 1;
                predicted_tfbs = predicted_tfbs(end:-1:1);
            end
            %         [length(predicted_tfbs) length(true_tfbs)]
            nTP = nTP + sum(predicted_tfbs(true_tfbs>0));
            nTN = nTN + sum((1-predicted_tfbs).*(1-true_tfbs));
            nFP = nFP + sum(predicted_tfbs.*(1-true_tfbs));
            nFN = nFN + sum((1-predicted_tfbs).*true_tfbs);
            
            if sum(predicted_tfbs.*true_tfbs) >= 0.25
                sTP = sTP + 1;
            else
                sFP = sFP + 1;
                sFN = sFN + 1;
            end
        end
    end
    nCC_1 = (nTP*nTN - nFP*nFN)/((nTP+nFN)*(nTN+nFP)*(nTP+nFP)*(nTN+nFN))^0.5;
    sSn_1 = sTP/(sTP + sFN);
    
    sSn(iter) = sSn_1;
    nCC(iter) = nCC_1;
    
    
    [sSn_1 nCC_1]
end
% close(h)
