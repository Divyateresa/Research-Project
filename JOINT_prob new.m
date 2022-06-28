clc;
clear all;
close all;
count = 0;

% for function F = A*B 
B = [0,1,1;1,0,1]
A = ones(2,3)
f = @(x) x
X = 1:10
for i = 1:length(B)-1
    for j = 1:length(B)
        b1(i,j) = B(i,j)*A(i,j)
if b1(i,j) == 0
    count = count+1;
end
    end 
end
%%% total entities transfer from A
total = numel(A)-count 
for x = 1:length(X)
%%% Probability distribution table 1 
pa_b1 = [f(x) 1/2-f(x); 1/2-f(x) f(x)]
A1_0 = sum(pa_b1(1,:))
A1_1 = sum(pa_b1(2,:))
P1_A = A1_0+A1_1; %% total probability of A

B1_0 = sum(pa_b1(:,1))
B1_1 = sum(pa_b1(:,2))
P1_B = B1_0+B1_1 %% total probability of B

%%% Entropy of H(X)
H1_A(x,:) = A1_0*log2(A1_0)+1-(A1_1)*log2(1-(A1_1))
H1_B(x,:) = B1_0*log2(B1_0)+1-(B1_1)*log2(1-(B1_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%

%%% Probability distribution table 2
pa_b2 = [f(x) 1/2; 1/2-2*f(x) f(x)]
A2_0 = sum(pa_b2(1,:))
A2_1 = sum(pa_b2(2,:))
P2_A = A2_0+A2_1; %% total probability of X

B2_0 = sum(pa_b2(:,1))
B2_1 = sum(pa_b2(:,2))
P2_B = B2_0+B2_1 %% total probability of Y

%%% Entropy of H(X)
H2_A(x,:) = (A2_0*log2(A2_0)+1-(A2_1)*log2(1-(A2_1)))
H2_B(x,:) = (B2_0*log2(B2_0)+1-(B2_1)*log2(1-(B2_1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%

%%% probability distribution table 3
pa_b3 = 0.5*pa_b1+((0.5)*pa_b2)
A3_0 = sum(pa_b3(1,:))
A3_1 = sum(pa_b3(2,:))
P3_A = A3_0+A3_1; %% total probability of X

B3_0 = sum(pa_b3(:,1))
B3_1 = sum(pa_b3(:,2))
P3_B = B3_0+B3_1 %% total probability of

%%% Entropy of H(X)
H3_A(x,:) = (A3_0*log2(A3_0)+1-(A3_1)*log2(1-(A3_1)))
H3_B(x,:) = abs(-B3_0*log(B3_0)-(B3_1)*log(B3_1))

%%% Entropy of H(A/B =1 )
H322(x,:)  = pa_b3(2,2)
H_AB(x,:)  = (H322(x,:))./(B3_1)
% H_AB3 = (H_AB).*log(H_AB)



%H_AB2(x,:) = (P_a_1_b_1)*log(P_a_1_b_1)
%%%%%%%%%%%%%%%5%5%5%5%%% % 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(X,total*H1_A,'-o','MarkerIndices',1:20:length(H1_A),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
hold on
plot(X,H2_A,'-o','MarkerFaceColor',[1,1,1],'MarkerIndices',1:20:length(H2_A),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
plot(X,H_AB,'-o','MarkerFaceColor',[1,1,1],'MarkerIndices',1:20:length(H_AB),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
hold on
xlabel('number of iteration')
ylabel('entropy')
title('communication cost of A for function A*B')
legend('H(Pab(1))','H(Pab(2))','H(Pab(3))')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


B1 = [0,1,1;1,0,1]
A1 = ones(2,3)
f = @(x) x
count1 = 0
% for function F = A+B 
for i = 1:length(B1)-1
    for j = 1:length(B1)
        B1(i,j) = A1(i,j)
if B1(i,j) == A1(i,j)
    count1 = count1+1;
end
    end 
end

%%it 
%%% total entities transfer from A
total1 = numel(A)-count1 
for x = 1:length(X)
%%% Probability distribution table 1 
pa_b1 = [f(x) 1/2-f(x); 1/2-f(x) f(x)]
A1_0 = sum(pa_b1(1,:))
A1_1 = sum(pa_b1(2,:))
P1_A = A1_0+A1_1; %% total probability of A

B1_0 = sum(pa_b1(:,1))
B1_1 = sum(pa_b1(:,2))
P1_B = B1_0+B1_1 %% total probability of B

%%% Entropy of H(X)
H1_A(x,:) = A1_0*log(A1_0)+1-(A1_1)*log(1-(A1_1))
H1_B(x,:) = B1_0*log(B1_0)+1-(B1_1)*log(1-(B1_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%

%%% Probability distribution table 2
pa_b2 = [f(x) 1/2; 1/2-2*f(x) f(x)]
A2_0 = sum(pa_b2(1,:))
A2_1 = sum(pa_b2(2,:))
P2_A = A2_0+A2_1; %% total probability of X

B2_0 = sum(pa_b2(:,1))
B2_1 = sum(pa_b2(:,2))
P2_B = B2_0+B2_1 %% total probability of Y


%%% Entropy of H(X)
H2_A(x,:) = (A2_0*log(A2_0)+1-(A2_1)*log(1-(A2_1)))
H2_B(x,:) = (B2_0*log(B2_0)+1-(B2_1)*log(1-(B2_1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%

%%% probability distribution table 3
pa_b3 = 0.5*pa_b1+((0.5)*pa_b2)
A3_0 = sum(pa_b3(1,:))
A3_1 = sum(pa_b3(2,:))
P3_A = A3_0+A3_1; %% total probability of X

B3_0 = sum(pa_b3(:,1))
B3_1 = sum(pa_b3(:,2))
P3_B = B3_0+B3_1 %% total probability of Y

%%% Entropy of H(X)
H3_A(x,:) = (A3_0*log(A3_0)+1-(A3_1)*log(1-(A3_1)))
H3_B(x,:) = (B3_0*log(B3_0)+1-(B3_1)*log(1-(B3_1)))

%%% Entropy of H(A/B =1 )
H3B1(x,:)  = pa_b3(2,2)
H_AB(x,:)  =(H3B1(x,:))./(B3_1)
%%%%%%%%%%%%%%%5%5%5%5%%% % 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(X,H1_A,'-o','MarkerIndices',1:20:length(H1_A),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
hold on
plot(X,H2_A,'-o','MarkerFaceColor',[1,1,1],'MarkerIndices',1:20:length(H2_A),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
plot(X,H_AB,'-o','MarkerFaceColor',[1,1,1],'MarkerIndices',1:20:length(H_AB),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
    
xlabel('number of iteration')
ylabel('entropy')
title('communication cost of A for function A=B')
legend('H(Pab(1))','H(Pab(2))','H(Pab(3))')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B1 = [0,1,1;1,0,1]
A1 = ones(2,3)
f = @(x) x
count2 = 0
% for function F = A+B 
for i = 1:length(B1)-1
    for j = 1:length(B1)
        b2(i,j) = B1(i,j)+A1(i,j)
if b2(i,j) == 0
    count2 = count2+1;
end
    end 
end

%%it 
%%% total entities transfer from A
total2 = numel(A)-count2
for x = 1:length(X)
%%% Probability distribution table 1 
pa_b1 = [f(x) 1/2-f(x); 1/2-f(x) f(x)]
A1_0 = sum(pa_b1(1,:))
A1_1 = sum(pa_b1(2,:))
P1_A = A1_0+A1_1; %% total probability of A

B1_0 = sum(pa_b1(:,1))
B1_1 = sum(pa_b1(:,2))
P1_B = B1_0+B1_1 %% total probability of B

%%% Entropy of H(X)
H1_A(x,:) = A1_0*log(A1_0)+1-(A1_1)*log(1-(A1_1))
H1_B(x,:) = B1_0*log(B1_0)+1-(B1_1)*log(1-(B1_1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%

%%% Probability distribution table 2
pa_b2 = [f(x) 1/2; 1/2-2*f(x) f(x)]
A2_0 = sum(pa_b2(1,:))
A2_1 = sum(pa_b2(2,:))
P2_A = A2_0+A2_1; %% total probability of X

B2_0 = sum(pa_b2(:,1))
B2_1 = sum(pa_b2(:,2))
P2_B = B2_0+B2_1 %% total probability of Y

%%% Entropy of H(X)
H2_A(x,:) = (A2_0*log(A2_0)+1-(A2_1)*log(1-(A2_1)))
H2_B(x,:) = (B2_0*log(B2_0)+1-(B2_1)*log(1-(B2_1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%

%%% probability distribution table 3
pa_b3 = 0.5*pa_b1+((0.5)*pa_b2)
A3_0 = sum(pa_b3(1,:))
A3_1 = sum(pa_b3(2,:))
P3_A = A3_0+A3_1; %% total probability of X

B3_0 = sum(pa_b3(:,1))
B3_1 = sum(pa_b3(:,2))
P3_B = B3_0+B3_1 %% total probability of Y

%%% Entropy of H(X)
H3_A(x,:) = (A3_0*log(A3_0)+1-(A3_1)*log(1-(A3_1)))
H3_B(x,:) = (B3_0*log(B3_0)+1-(B3_1)*log(1-(B3_1)))

%%% Entropy of H(A/B =1 )
H3B1(x,:)  = pa_b3(1,2)
H_AB(x,:)  = (H322(x,:))./(B3_1)
%%%%%%%%%%%%%%%5%5%5%5%%% % 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(X,H1_A,'-o','MarkerIndices',1:20:length(H1_A),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
hold on
plot(X, H2_A,'-o','MarkerFaceColor',[1,1,1],'MarkerIndices',1:20:length(H2_A),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')
plot(X,total2*H_AB,'-o','MarkerFaceColor',[1,1,1],'MarkerIndices',1:20:length(H_AB),'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','b')

xlabel('number of iteration')
ylabel('entropy')
title('communication cost of A for function A+B')
legend('H(Pab(1))','H(Pab(2))','H(Pab(3))')
