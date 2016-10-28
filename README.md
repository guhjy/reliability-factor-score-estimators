# Reliability Estimates for Three Factor Score Estimators
The scripts in this repository accompany our paper (Beauducel, Harms & Hilger, 2016) and reflect the appendices A and B of the publication.

## Example for R script
Calculates the reliability estimates for a given factor analysis solution:
```R
Loadings <- matrix(c(
	 0.50,-0.10, 0.10,
 	 0.50, 0.10, 0.10,
 	 0.50, 0.10,-0.10,
 	-0.10, 0.50, 0.15,
 	 0.15, 0.50, 0.10,
 	-0.15, 0.50, 0.10,
 	 0.10, 0.10, 0.60,
 	 0.10,-0.10, 0.60,
 	 0.10, 0.10, 0.60
  	),
 	nrow=9, ncol=3,
 	byrow=TRUE)
 InterCorr <- matrix(c(
 	1.00, 0.30, 0.20,
  	0.30, 1.00, 0.10,
  	0.20, 0.10, 1.00
 	),
 	nrow=3, ncol=3,
 	byrow=TRUE)

reliabilities <- factor.score.reliability(Lambda = Loadings, Phi = InterCorr, Estimators = c("Regression", "Bartlett", "McDonald"))
lapply(reliabilities, round, 3)
```

**Remark:** If a factor solution was obtained through the `fa` package, the factor loadings from the solution can be passed through `Lambda`.

## Reference
Beauducel, A., Harms, C., & Hilger, N. (2016). Reliability estimates for three factor score predictors. International Journal of Statistics and Probability, 5(6), 94–107. http://doi.org/10.5539/ijsp.v5n6p94
