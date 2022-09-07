Question 1 looks at how to model a radar detector using MAP decision
rules.  Part A finds the optimal decision boundary for detection and
finds the theoretical probability of error given the decision rule.
These values were then compared to the theoretical values that were found
mathematically.  Part B then plots the receiver operating curve for
differing signal to noise ratios.  As expected, the higher the signal to
noise ratio, the better the detection algorithm performed.  Part C
specifies missing a target is 10 times worse than a false positive
judgement.  This updated our part A answers and the differing decision
boundaries were displayed on a new ROC graph.  The new decision boundary
had a higher tolerance for false positive in comparison to Part A, which
logically agrees with the updated cost. Part E causes the target to no longer
have a unique mean, but rather a unique variance.  This caused there to
be two decision boundaries, with the target classified in the middle region.
A new decision boundary, error, and ROC function were created to settle
this case.  The experimental values again agreed with the theoretical
results.  Question 2 introduces a new dataset that utilized
likelihood and prior information for classification.  The data
was first split up into a testing and training dataset. The training
dataset used it's known classes to make a 4D pdf for each case.  The
testing data was then input into each pdf and multiplied by the corresponding
bayesian prior to get it's probability of being a member of that class.
The class with the highest probability was selected and the corresponding
confusion matrix was plotted.
