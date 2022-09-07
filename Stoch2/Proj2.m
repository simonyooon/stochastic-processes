%Michael Bentivegna, Joya Debi, Simon Yoon
%ECE310 Stochastic Processes Project 2: MMSE Estimation

clc
clear
close all

% This project discusses how to find and estimate the minimum mean square error
% of random variables. Part 1 analyzes how to find both the linear MMSE and
% Bayesian MMSE estimators when analyzing one observation from a summation of two
% uniformly distributed random variables. Part 2 finds the linear estimator
% for a similar situation except with an arbitrary number of observations.
% Clearly, as the number of noisy observations for a given Y value increase, so
% does the accuracy of the guess. This subsequently causes a decrease in
% MSE value of the MMSE estimator.  The variance also plays a role in the MSE
% as a low variance provides smaller discrepancies in estimates and thus
% a lower MSE. This information is corroborated in the generated graph in Part 2.

%% Part 1 - Bayes MMSE and linear MMSE 8.5 and 8.6

%Declarations
N = 10^6;
Y = -1 + (2)*rand(1,N);
W = -2 + (4)*rand(1,N); 
X = Y + W;
yHat = X; %Temporary 

%Estimates for yHat based on X (Equation 8.30)
yHat(X<-1) = 0.5+X(X<-1)*.5;
yHat(X>1) = -0.5+X(X>1)*.5;
yHat(X>-1 & X<1) = 0;

%Bayes MMSE estimate
bayes = mean((Y-yHat).^2);
bayesMSE = mean(bayes);

%LMMSE estimate (Equation 8.51)
linear = mean(Y) + (var(Y)/var(X)).*(X - mean(X));
linearMSE = mean((Y-linear).^2);

%Given theoretical values
theobayesMSE = 1/4;
theolinearMSE = 4/15;

%Show on table
table = table([theobayesMSE;theolinearMSE],[bayesMSE;linearMSE],'VariableNames', ...
    {'Theoretical','Experimental'},'RowNames',{'Bayes MSE','Linear MSE'});
disp(table)

%% Part 2 - Estimation from an arbitrary number of noisy measurements 8.8

%Declarations
observations = 20;
sigmaR = [.5, 1, sqrt(2)];
sigmaY = [.5, 1, sqrt(2)];

%Each row has a unique variance
%Each column has a unique # of observations
theoretical = zeros(3, observations);
experimental = zeros(3, observations);

for j = 1:3
    for i = 1:observations  
        [theoretical(j,i), experimental(j,i)] =  twoNoisyEstimates(i, sigmaR(j), sigmaY(j));
    end
end

%Plotting
figure(1);
hold on;
for j = 1:3
    plot(1:observations, theoretical(j, :))
    scatter(1:observations, experimental(j, :), 'X')
end
title("Minimum Mean Squared Error Plot")
xlabel("# of Observations for a given Y");
ylabel("MMSE");
legend("THEORETICAL: Var = .25", "EXPERIMENTAL: Var = .25", "THEORETICAL: Var = 1", "EXPERIMENTAL: Var = 1", "THEORETICAL: Var = 2", "EXPERIMENTAL: Var = 2")

%Function for finding both theoretical and experimental values of the model
function [theoretical, experimental] = twoNoisyEstimates(observations, sigmaR, sigmaY)
   
    %Number of trials for a specific number of observations
    N = 10^4;
    
    %Theoretical MMSE (Equation 8.80)
    theoretical = (sigmaY^2 * sigmaR^2) / (observations * sigmaY^2 + sigmaR^2);

    %Experimental Values
    Y = normrnd(1, sigmaY, [N 1]);
    R = normrnd(0, sigmaR, [N observations]);
    X = zeros(N, observations);
    for i = 1:observations
       X(:,i) = R(:,i) + Y; 
    end

    %Find MMSE using experimental data (Equation 8.79)
    varR = mean(var(R'));
    yH = (1 / (observations * var(Y) + varR)) * (varR * mean(Y) + var(Y) * sum(X, 2));
    experimental = mean((Y-yH).^2);
end

