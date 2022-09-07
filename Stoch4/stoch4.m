% Michael Bentivegna, Simon Yoon, Joya Debi
% ECE302 Stochastic Processes Project 4: Detection

clear;
clc;
close all;

% Summary of Project

% Question 1 looks at how to model a radar detector using MAP decision
% rules.  Part A finds the optimal decision boundary for detection and
% finds the theoretical probability of error given the decision rule.
% These values were then compared to the theoretical values that were found
% mathematically.  Part B then plots the receiver operating curve for
% differing signal to noise ratios.  As expected, the higher the signal to
% noise ratio, the better the detection algorithm performed.  Part C
% specifies missing a target is 10 times worse than a false positive
% judgement.  This updated our part A answers and the differing decision
% boundaries were displayed on a new ROC graph.  The new decision boundary
% had a higher tolerance for false positive in comparison to Part A, which
% logically agrees with the updated cost. Part E causes the target to no longer
% have a unique mean, but rather a unique variance.  This caused there to
% be two decision boundaries, with the target classified in the middle region.
% A new decision boundary, error, and ROC function were created to settle
% this case.  The experimental values again agreed with the theoretical
% results.  Question 2 introduces a new dataset that utilized
% likelihood and prior information for classification.  The data
% was first split up into a testing and training dataset. The training
% dataset used it's known classes to make a 4D pdf for each case.  The
% testing data was then input into each pdf and multiplied by the corresponding
% bayesian prior to get it's probability of being a member of that class.
% The class with the highest probability was selected and the corresponding
% confusion matrix was plotted.

%% Question 1: Detection Algorithm

% ***********Part A*************

% Theoretical Values
fprintf('Theoretical Decision Boundary for Part A = 1.693. For any samples greater than this value the target will be predicted to be present.\n')
fprintf('Theoretical Probability of Error for Part A = .112 \n\n')

% The theoretical decision boundary was found by plotting the two gaussian functions
% with their prior likelihood and seeing where they intersected.

% The theoretical probability of error was found by finding the area under each curve
% where either a false positive or false negative could occur. This was
% done through definite integral calculation.

% Known Values
A = 2;
Z = 0;
var = 1;
probWithout = .8;

% Call Decision Region Function
decisionRegion = findDecisionRegion(Z, var, A, var, probWithout, 1);
fprintf('Experimental Decision Region for Part A = %f. For any samples greater than this value the target will be predicted to be present.\n', decisionRegion)

% Call Experimental Probability of Error Function
error = experimentalProbError(Z, var, A, var, probWithout, decisionRegion);
fprintf('Experimental Probability of Error for Part A = %f \n', error)

% ***********Part B*************

%Call function for different ROC values (note the zero in the last parameter this will be used in part c)
[falsePos1, truePos1] = ROCcurve(Z, .5*var, A, .5*var, probWithout, 0);
[falsePos2, truePos2] = ROCcurve(Z, var, A, var, probWithout, 0);
[falsePos3, truePos3] = ROCcurve(Z, 2*var, A, 2*var, probWithout, 0);

% Plotting ROC for different SNR
figure(1)
hold on;
plot(falsePos1, truePos1)
plot(falsePos2, truePos2)
plot(falsePos3, truePos3)
legend('SNR = 4', 'SNR = 2', 'SNR = 1')
xlabel("False Positive Rate")
ylabel("True Positive Rate")
title("ROC Curve for Different SNRs")

% ***********Part C*************

% Now input a cost other than 1 in the decision region function
decisionRegionC = findDecisionRegion(Z, var, A, var, probWithout, 10);
[falsePosC, truePosC, specificSampleCoord] = ROCcurve(Z, var, A, var, probWithout, [decisionRegionC, decisionRegion]);


% Plotting ROC w/ optimal Coordinate Labelled
figure(2)
hold on;
plot(falsePosC, truePosC)
plot(specificSampleCoord(1,1), specificSampleCoord(1,2), 'o')
plot(specificSampleCoord(2,1), specificSampleCoord(2,2), 'o')
xlabel("False Positive Rate")
ylabel("True Positive Rate")
title("ROC Curve")
legend("ROC Curve w/ SNR = 2", "Optimal Marker w/ C_{01} = 10*C_{10}", "Optimal Marker w/ C_{01} = C_{10}")

% ***********Part D*************

% Omitted

% ***********Part E*************

%Theoretically Calculated Values (found the same was as in part A)
fprintf('\n***************************\n\n')
fprintf('Theoretical Decision Region for Part E = -1.187 to 1.187.  For any samples in this region the target will be predicted to be present. \n')
fprintf('Theoretical Probability of Error for Part E = .1414 \n\n')
% Known Values
Ae = 0;
varWithout = 8;
varTarget = 1;
probWithoute = .8;

% Call decision region function
decisionRegionE = findDecisionRegionE(Ae, varWithout, varTarget, probWithoute);
fprintf('Experimental Decision Region for Part E = %f to %f. For any samples in this region the target will be predicted to be present.\n', decisionRegionE(1), decisionRegionE(2))

% Call error function
errorE = experimentalProbErrorE(Ae, varWithout, varTarget, probWithoute, decisionRegionE);
fprintf('Experimental Probability of Error for Part E = %f \n', errorE)

%Call ROC curve functions
[falsePosE1, truePosE1] = ROCcurveE(Ae, 16, 1, probWithoute);
[falsePosE2, truePosE2] = ROCcurveE(Ae, 8, 1, probWithoute);
[falsePosE3, truePosE3] = ROCcurveE(Ae, 4, 1, probWithoute);

%Plotting
figure(3)
hold on;
plot(falsePosE1, truePosE1)
plot(falsePosE2, truePosE2)
plot(falsePosE3, truePosE3)
legend('Var Ratio = 16', 'Var Ratio = 8', 'Var Ratio = 4')
xlabel("False Positive Rate")
ylabel("True Positive Rate")
title("ROC Curve for Different Variances")
hold off;


%% Question 2:

% Load in Iris data (comes with features and labels variables)
load('Iris.mat');

% 150 4D samples
dimOfData = size(features);
samples = dimOfData(1);

% 3 Unique labels (1, 2, 3)
uniqueLabels = length(unique(labels));

% Split data into training and testing
halfSamples = samples/2;
shuffling = transpose(randperm(samples));

training = shuffling(1:halfSamples);
testing = shuffling(halfSamples+1:samples);

%Get training and testing data
trainFeatures = features(training, :);
trainLabels = labels(training, :);

testFeatures = features(testing, :);
testLabels = labels(testing, :);

% Classify each 4D sample with its correct label in the training data and get its mean and covariance
featuresForOne = trainFeatures((trainLabels == 1),:);
featuresForTwo = trainFeatures((trainLabels == 2),:);
featuresForThree = trainFeatures((trainLabels == 3),:);

meanOne = mean(featuresForOne);
meanTwo = mean(featuresForTwo);
meanThree = mean(featuresForThree);

covOne = cov(featuresForOne);
covTwo = cov(featuresForTwo);
covThree = cov(featuresForThree);

% Get prior values for each label (since the labels are 1,2,3 they are also in the same spot as their index
bayesianPriors = histcounts(testLabels) / length(testLabels);

% Probability of test features to come from each of the three labels using training data's mean and covariance
probOfEach = zeros([halfSamples, uniqueLabels]);
probOfEach(:,1) = mvnpdf(testFeatures, meanOne, covOne) * bayesianPriors(1);
probOfEach(:,2) = mvnpdf(testFeatures, meanTwo, covTwo) * bayesianPriors(2);
probOfEach(:,3) = mvnpdf(testFeatures, meanThree, covThree) * bayesianPriors(3);

% Choose the best label for each samples and see the error rate
[maximum, columnIndex] = max(probOfEach, [], 2);
errorProb = 1 - mean(columnIndex == testLabels);
fprintf("\nThe Probability of Classification Error: %f", errorProb);

%Plot Confusion Matrix
figure(4);
confusionchart(confusionmat(columnIndex, testLabels));
title('Confusion Matrix for Iris Data');

%% Functions for Question 1

% Find border of decision regions for Part 1A and 1C
function decisionRegion = findDecisionRegion(mu0, var0, mu1, var1, probWithout, costMiss)
    % Setup
    spacing = 10000;
    x = linspace(-10, 10, spacing);
    probWith = 1 - probWithout;
    
    % Get pdf when target is present and when it is not present
    yWithoutDetection = normpdf(x, mu0, var0)*probWithout;
    yWithDetection = normpdf(x, mu1, var1)*probWith*costMiss;
    
    % Find where the threshold value is by seeing where the functions intersect
    for i = 1:spacing
        if (yWithoutDetection(i) <= yWithDetection(i))
            decisionRegion = x(i);
            break;
        end
    end
end

%Find experimental probability of error for the decision region in Part 1A
function error = experimentalProbError(mu0, var0, mu1, var1, probWithout, decisionRegion)
    % Setup
    totalSamples = 100000;
    
    % Generate proper amount of samples from each distribution
    nonTargetSamples = normrnd(mu0, var0, [1, round(probWithout*totalSamples)]);
    targetSamples = normrnd(mu1, var1, [1, round((1-probWithout)*totalSamples)]);
    

    % Check to see how often samples are correctly categorized
    error = 1 - ((length(targetSamples(targetSamples > decisionRegion)) + length(nonTargetSamples(nonTargetSamples < decisionRegion))) / totalSamples);
end

% Find ROC arrays for Part 1A conditions
function [falsePositiveRate, truePositiveRate, specificSampleCoord] = ROCcurve(mu0, var0, mu1, var1, probWithout, decisionRegion)
    % Setup
    totalSamples = 10000;
    spacing = 1000;
    x = linspace(-10, 10, spacing);
    falsePositiveRate = zeros([1, 100]);
    truePositiveRate = zeros([1, 100]);
    specificSampleCoord = zeros([length(decisionRegion), 2]);
    
    % Generate proper amount of samples from each distribution
    nonTargetSamples = normrnd(mu0, var0, [1, round(probWithout*totalSamples)]);
    targetSamples = normrnd(mu1, var1, [1, round((1-probWithout)*totalSamples)]);
    
    %Run through all possible decision region choices to get ROC curve
    for i = 1:spacing
            falsePositiveRate(i) = length(nonTargetSamples(nonTargetSamples > x(i)))/ (round(probWithout*totalSamples));
            truePositiveRate(i) = length(targetSamples(targetSamples > x(i))) / (round((1-probWithout)*totalSamples));
    end
    
    %Use specific coordinate for the decision regions (mostly for part c)
    for i = 1:length(decisionRegion)
        specificSampleCoord(i, 1) = length(nonTargetSamples(nonTargetSamples > decisionRegion(i)))/ (round(probWithout*totalSamples));
        specificSampleCoord(i, 2) = length(targetSamples(targetSamples > decisionRegion(i))) / (round((1-probWithout)*totalSamples));
    end
end

% Find border of decision regions for Part 1E
function decisionRegion = findDecisionRegionE(mu0, var0, var1, probWithout)
    % Setup
    spacing = 10000;
    x = linspace(-10, 10, spacing);
    probWith = 1 - probWithout;
    
    % Get pdf when target is present and when it is not present
    yWithoutDetection = normpdf(x, mu0, var0)*probWithout;
    yWithDetection = normpdf(x, mu0, var1)*probWith;
    
    % Find where the threshold value is by seeing where the functions intersect
    decisionRegion = zeros([1,2]);
    for i = 1:spacing
        if (yWithoutDetection(i) <= yWithDetection(i))
            decisionRegion(2) = x(i);
        end
        
        if (yWithoutDetection(spacing + 1 - i) <= yWithDetection(spacing + 1 - i))
            decisionRegion(1) = x(spacing + 1 - i);
        end
    end
end

%Find experimental probability of error for the decision region in Part 1E
function error = experimentalProbErrorE(mu0, var0, var1, probWithout, decisionRegion)
    % Setup
    totalSamples = 100000;
    
    % Generate proper amount of samples from each distribution
    nonTargetSamples = normrnd(mu0, var0, [1, round(probWithout*totalSamples)]);
    targetSamples = normrnd(mu0, var1, [1, round((1-probWithout)*totalSamples)]);
    

    % Check to see how often samples are correctly categorized
    error = (length(nonTargetSamples((nonTargetSamples < decisionRegion(2)) & (nonTargetSamples > decisionRegion(1)))) + length((targetSamples(targetSamples > decisionRegion(2)))) + length((targetSamples(targetSamples < decisionRegion(1))))) / totalSamples;

end

% Find ROC arrays for Part 1E conditions
function [falsePositiveRate, truePositiveRate] = ROCcurveE(mu0, var0, var1, probWithout)
    % Setup
    totalSamples = 10000;
    spacing = 1000;
    x = linspace(-10, 10, spacing);
    falsePositiveRate = zeros([1, 100]);
    truePositiveRate = zeros([1, 100]);
    
    % Generate proper amount of samples from each distribution
    nonTargetSamples = normrnd(mu0, var0, [1, round(probWithout*totalSamples)]);
    targetSamples = normrnd(mu0, var1, [1, round((1-probWithout)*totalSamples)]);
    
    %Run through all possible decision region choices to get ROC curve
    for i = 1:spacing
            falsePositiveRate(i) = length(nonTargetSamples(abs(nonTargetSamples) < abs(x(i))))/ (round(probWithout*totalSamples));
            truePositiveRate(i) = length(targetSamples(abs(targetSamples) < abs(x(i)))) / (round((1-probWithout)*totalSamples));
    end
    
end
