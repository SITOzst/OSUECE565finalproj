%ECE 565 final project%
%Medical Imaging%
clear;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Baisc setting 
%%% parameter:
% % m: block number
% % n: observation number
% % A: probability matrix that denote how many particle may across the area
% % theta: denote the partical across the block after take log(other domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take an easy example - 9 block area
n = 16;  
m = 9;  
%normA = [1, 2, 5, 10,12, 15,20]; % initialize the             
%normA = [50,100,200,500,1000,5000,10000]; % initialize the     
normA = [0.01,0.05,0.1,0.5,1,5,10];%
mses = ones(1,7); % initialize set of MSEs in each iteration (to normA)--container %
mses_mle = ones(1,7); % --container
vars = ones(1,7);% --container                                                   %
crlbs = ones(1,7);% initialize set of CRLBs in each iteration (to normA)  --container %  
ii = 1; % index of the position of value in each iteration                %
A = rand(n,m); % random form matrix A    
%theta = rand(m,1); % random form theata
theta = [1,0,1,0,1,0,1,0,1]'%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% loop start 
% % using 200 Monte-Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_norm = normA
% % %     write the general form of theta-j
    Ait = A * i_norm;   
    y = poissrnd(Ait*theta) % random generate observation
    theta0 = ones(m,1); %init
    errs = zeros(1,200); %init
    errsmle = zeros(1,200); %init
    fims = zeros(1,200); %init
    for kp = 1:200 
        for j = 1:m
            sum1 = 0;
            for i = 1:n 
                sum2 = 0;
                for k = 1:m % sum of all blocks
                    tem2 = Ait(i,k) * theta0(k);
                    sum2 = sum2 + tem2;
                end       
                tem = (y(i)*Ait(i,j))/sum2;
                sum1 = sum1 + tem;
            end
            theta0(j) = (theta0(j)/sum(Ait(:,j)))*sum1;
        end
        errs(kp) = immse(theta, theta0);

        fim_tem = 0;
        %{
        for icr = 1:n
            sum3 = 0;
            for jcr = 1:j
                tem3 = Ait(icr,jcr) * theta0(jcr);
                sum3 = sum3 + tem3;
            end
            Atem = Ait';
            tem4 = (Atem(icr,:) * transpose(Atem(icr,:)))/sum3;
            fim_tem = fim_tem + tem4;
        end
       %}
        for icr = 1:m
            Atem = Ait';
            tem5 = A*theta;
            fim_tem = Atem(icr,:)*Atem(icr,:)'/tem5(icr);
        end
        
       %{
        fim_tem = Ait'*diag(1/(Ait*theta0))*Ait;
        %}
        fims(kp) = fim_tem;
        
        
        cvx_begin quiet
        variable thetamle(m)
        minimize( sum(Ait*thetamle)-y'*log(Ait*thetamle))
        cvx_end;
        errsmle(kp) = immse(theta, thetamle);
    end
theta0
%figure()
plot(errs)        
mse_emperical = mean(errs);
mse_mle = mean(errsmle);
fim = mean(fims);
vars(ii) = var(theta);
mses(ii) = mse_emperical;
mses_mle(ii) = mse_mle;
crlbs(ii) = inv(fim);
ii = ii + 1;
end
mses
figure;
plot(normA,crlbs,'-s',normA,mses,'-d',normA,vars,'-p',normA,mses_mle,'-r');
xlabel('coefficient of A ');ylabel('CRLB and MSE-EM and VAR and MSE-ML');
leg2=legend('CRLB','MSE-EM','Varience','MSE-ML','Location','SouthEast');
grid on
figure;
loglog(normA,crlbs,'-s',normA,mses,'-d',normA,vars,'-p',normA,mses_mle,'-r');
xlabel('coefficient of A');ylabel('CRLB and MSE and VAR and MSE-ML');
leg2=legend('CRLB','MSE','Varience','MSE-ML','Location','SouthEast');
grid on


% figure;
% loglog(normA,crlbs,'-d');
% xlabel('coefficient of A');ylabel('CRLB and MSE and VAR ');
% grid on

%{
f=imread('images.jpg');
I = rgb2gray(f);
sizeI = size(I);
PPPP = I(1,:);

for iim = 2:sizeI(1)
    PPPP = [PPPP,I(iim,:)];
end
n = 10;
m = 16384;
A = rand(n,m); % random form matrix A
theta = double(PPPP);
theta0 = ones(m,1)*100;
errs = zeros(1,5);
for kp = 1:5
    for j = 1:m
        sum1 = 0;
        for i = 1:n
            sum2 = 0;
            for k = 1:m
                tem2 = A(i,k) * theta0(k);
                sum2 = sum2 + tem2;
            end       
            tem = (y(i)*A(i,j))/sum2;
            sum1 = sum1 + tem;
        end
        theta0(j) = (theta0(j)/sum(A(:,j)))*sum1;
    end
    errs(kp) = immse(theta', theta0);
end
%}