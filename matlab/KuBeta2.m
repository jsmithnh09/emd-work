clear; close all;clc

B=.1:0.0001:5; % Theoretical Beta range with small step size of 0.0001
Ku=((gamma(5./B).*gamma(1./B))./(gamma(3./B).^2))-3; % GGD theoretical Excess Kurtosis


N = 2e3; % I used here 2000 samples
mu = 0;
alpha = 1;
sigma = 1;

beta = 1:0.1:4; % shape parameter range to be tested

for i = 1:length(beta)
    i
    X = rand_ggd(1,N,mu,alpha,beta(i)); % I use this GGD generator
    Y= ggdrnd(mu, sigma, beta(i), [N,1]); % GGD generator from your GitHub

    TT=length(X); % X and Y have the same length

    KK4X=TT*(sum((X-mean(X)).^4)/sum((X-mean(X)).^2)^2)-3;  % Kurtosis of X
    KK4Y=TT*(sum((Y-mean(Y)).^4)/sum((Y-mean(Y)).^2)^2)-3;  % Kurtosis of Y

    shape_NewtonX(i) = gcm_search(X,5,100);
    shape_NewtonY(i) = gcm_search(Y,5,100);


    [~,Kindex]=min(abs(Ku-KK4X));
    shape_KurtosisX(i) = B(Kindex);

    [~,Kindex2]=min(abs(Ku-KK4Y));
    shape_KurtosisY(i) = B(Kindex2);
end

figure,plot(beta,shape_NewtonX),hold, plot(beta,shape_KurtosisX,'r'),title('Signal X')
legend('Newton shape','Kurtosis shape')
figure,plot(beta,shape_NewtonY),hold, plot(beta,shape_KurtosisY,'r'),title('Signal Y')
legend('Newton shape','Kurtosis shape')