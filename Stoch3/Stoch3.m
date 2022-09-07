% Michael Bentivegna, Simon Yoon, Joya Debi
% ECE310 Stochastic Processes Project 3: Maximum Likelihood Estimation

% This project utilizes maximum likelihood (ML) estimators to predict the
% value of the alpha and lambda parameters in the rayleigh and exponential
% distributions, respectively.  In part 1, data is generated from randomly
% drawing samples from ground truth distributions. The formula to find the
% ML estimator for each distribution was derived by hand before applying it
% to the data.  Once the estimator was calculated, it's corresponding mean
% squared error (MSE), bias, and variance was plotted for differing numbers of
% observations.  Note that for each number of observations 100000 trials
% was completed. In part 2, data was provided from either a rayleigh or
% exponential distribution.  Using the calculations in part 1, the
% parameter was estimated for each potential distribution. The log
% likelihood of these corresponding distributions was then found.  Since
% the rayleigh distribution had a greater log likelihood, the data
% most likely came from sampling a rayleigh pdf.

% DERIVED FORMULAS IN ATTACHED PDF

clear;
clc;
close all;

%% Part 1: Exponential and Rayleigh Samples
alphaOrLambda = [1, 2, 5]; % Ground Truth

% Arrays for each category (20 observations for each)
rayleighMSE = zeros([3, 20]);
rayleighBIAS = zeros([3, 20]);
rayleighVAR = zeros([3, 20]);
exponentialMSE = zeros([3, 20]);
exponentialBIAS = zeros([3, 20]);
exponentialVAR = zeros([3, 20]);

% Fill each matrix element with function value corresponding to number of
% observations and the value of alpha or lambda
for j = 1:3
    for i = 1:20
        [rayleighMSE(j, i), rayleighBIAS(j, i), rayleighVAR(j, i)] = rayleighDist(i, alphaOrLambda(j));
        [exponentialMSE(j, i), exponentialBIAS(j, i), exponentialVAR(j, i)] = exponentialDist(i, alphaOrLambda(j)); 
    end
end

% Plot MSE, Bias, and Variance
graphing(rayleighMSE, exponentialMSE, "Mean Squared Error")
graphing(rayleighBIAS, exponentialBIAS, "Bias")
graphing(rayleighMSE, exponentialMSE, "Variance")

%% Part 2: Find Best Graph for given data
givenData = load('data.mat').data;

% Rayleigh alpha estimator
alphaGuess = sqrt(.5.*(mean(givenData.^2)));

% Exponential lambda estimator
lambdaGuess =  1000 / (sum(givenData));

% Check each one's log likelihood value
rayLogLikelihood = sum(log(raylpdf(givenData, alphaGuess)));
expLogLikelihood = sum(log(exppdf(givenData, 1/lambdaGuess)));

% Since the Rayleigh Distribution's likelihood is greater that is most
% likely the distribution
fprintf('The Rayleigh Distribution Log Likelihood Value was: %f \n', rayLogLikelihood);
fprintf('The Exponential Distribution Log Likelihood Value was: %f \n', expLogLikelihood);
fprintf('Thus, these samples likely came from a Rayleigh Distribution! \n');

%% Functions
function graphing(data1, data2, Title)
    x = 4:23; % # of observations
    figure;
    subplot(1, 2, 1)
    plot(x, data1(1, :), x, data1(2, :), x, data1(3, :))
    title("Rayleigh Distribution " + Title);
    ylabel(Title);
    xlabel("# of Observations");
    xlim([0, 25])
    legend("alpha = 1", "alpha = 2", "alpha = 5")
    
    subplot(1, 2, 2)
    plot(x, data2(1, :), x, data2(2, :), x, data2(3, :))
    title("Exponential Distribution " + Title);
    ylabel(Title);
    xlabel("# of Observations");
    xlim([0, 25])
    legend("lambda = 1", "lambda = 2", "lambda = 5")
end

function [MSE, bias, variance] = rayleighDist(obs, alpha)
    N = 100000; % Number of trials
    ray = raylrnd(alpha, N, obs + 3); % Get randomly generated rayleigh samples
    avgObs = mean(ray.^2, 2); 
    alphaGuesses = sqrt(.5.*avgObs); % Get the alpha guesses for each trial
    
    % Analysis
    MSE = mean((alpha-alphaGuesses).^2);
    bias = mean(alphaGuesses) - alpha;
    variance = var(alphaGuesses);
end

function [MSE, bias, variance] = exponentialDist(obs, lambda)
    N = 100000; % Number of trials
    exp = exprnd(1/lambda, N, obs + 3); % Get randomly generated rayleigh samples
    sumObs = sum(exp, 2);
    lambdaGuesses = obs ./ (sumObs); % Get the lambda guesses for each trial
    
    % Analysis
    MSE = mean((lambdaGuesses-lambda).^2);
    bias = mean(lambdaGuesses) - lambda;
    variance = var(lambdaGuesses);
end
