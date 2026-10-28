The pairplot and boxplot show that there are some differences. The distribution of flavanoids is useful for separating groups, although the classes can still overlap. The colours and clusters in the pairplot suggest that combinations of features are more informative than any one measurement.

Some of the strongest positive correlations between features:
- Total phenols and flavanoids: r = 0.865.
- Flavanoids and OD280/OD315 of diluted wines: r = 0.787.
- Total phenols and OD280/OD315: r = 0.700.
- Alcohol and proline: r = 0.644.

negative correlations include:
- Malic acid and hue: r = -0.561.
- Flavanoids and nonflavanoid phenols: r = -0.538.
- Color intensity and hue: r = -0.522.

Principal component analysis was performed after standardising all 13 features. The first principal component explains 36.20% of the variance and the second explains 19.21%.


All three models performed better using all 13 features than using only two principal components. The random forest achieved the highest test accuracy on this particular split. Cross-validation would help to assess stability.

For the decision tree, the most important variables were flavanoids (0.408), color intensity (0.400), and proline (0.111). For the random forest, the leading variables were color intensity (0.188), flavanoids (0.160), proline (0.147), alcohol (0.118), and hue (0.102). The exact importance scores depend on the training data and model, and do not by themselves establish causal influence.


Research questions and hypotheses:
Question 1: Can phenolic compounds predict wine class?

Hypothesis: Flavanoids, total phenols and OD280/OD315 together should provide useful class information. They have strong correlations with one another, and flavanoids are important for both tree models. A future analysis could compare accuracy using only these features with accuracy using the full dataset.

