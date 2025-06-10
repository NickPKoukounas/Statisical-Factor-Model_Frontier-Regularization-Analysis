# Statistical Factor Models and Portfolio Regularization in Mathematica

This repository contains a full implementation of a statistical factor modeling and portfolio optimization framework using Wolfram Mathematica. It includes ETF return data processing, missing data imputation, factor model fitting, mean and covariance shrinkage via cross-validation, efficient frontier construction, and in-sample vs. out-of-sample evaluation.

---

## Investment Universe

The investment universe is a diversified set of ETFs

| Ticker | Description                      |
| ------ | -------------------------------- |
| EEM    | Emerging Markets Equity          |
| EFA    | Developed Markets ex-US & Canada |
| EWJ    | Japan Equity                     |
| IEF    | Intermediate US Treasuries       |
| IEV    | Europe Equity                    |
| IVV    | S\&P 500 Index Fund              |
| RWR    | Real Estate (REITs)              |
| SHY    | Short-term US Treasuries         |
| TLT    | Long-term US Treasuries          |
| VTI    | Total US Market Equity           |

Data is collected using `FinancialData` between April 2000 and March 2025.

---

##  Data Processing

1. **Daily Close Prices** are retrieved for all ETFs.
2. **Monthly Returns** are computed by selecting end-of-month prices and calculating simple returns:

   * Returns: $r_t = \frac{P_t}{P_{t-1}} - 1$
   * Aligned to a shared monthly calendar using `xExpandCalendar`.
3. **Missing Data Handling**:

   * Missing returns are estimated using `xMeanCovMissingMLE`, which outputs:

     * Completed return matrix
     * Estimated mean vector and covariance matrix
     * Log-likelihood progression and run ID
     * Must import this file to use its function
---

## Factor Model Estimation

* Factor model is fit using `xFactorFitMLE`, where the number of factors is selected via BIC minimization.
* Must import this file to use its function
* Initialization via `xInitializeFactorModel`.
* BIC values are computed and plotted.
* The best model is chosen

---

##  Training/Testing Split

* Random 70% of observations selected as training set.
* Remaining 30% used as testing set.
* 5-fold cross-validation on the training set for hyperparameter tuning.

---

## Regularization using Shrinkage

### Shrinkage of Mean Vector

* Shrinkage target: Grand mean of sample means.
* Shrinkage form: $\mu_{\text{shrink}} = (1 - \gamma) \mu + \gamma \bar{\mu}$
* 5-fold CV across $\gamma \in [0, 1]$ by 0.05 steps.
* Optimal $\gamma$ selected via minimal validation MSE.

### Shrinkage of Covariance Matrix

* Shrinkage target: Scaled identity $\lambda \cdot \bar{\sigma}^2 I$
* Shrinkage form: $\Sigma_{\text{shrink}} = (1 - \lambda) \Sigma + \lambda \bar{\sigma}^2 I$
* 5-fold CV over $\lambda \in [0, 1]$.
* Optimal $\lambda = 0.05$ selected.

---

## Efficient Frontier Construction

### Quadratic Optimization

* Solves: $\min_x \frac{1}{2} x^T \Sigma x$, s.t. $x \geq 0, \sum x = 1, x^T \mu = \tau$
* Frontier built by sweeping $\tau$ from min to max achievable return.

### Portfolios Computed:

* Minimum variance portfolio
* Maximum return portfolio
* Interpolated efficient portfolios
* Repeated for:

  * Regular sample estimates
  * Shrinkage estimates

---

##  In-Sample vs Out-of-Sample Comparison

* Mean and covariance estimated from testing set.
* Each frontier portfolio evaluated out-of-sample:

  * Risk (std dev): $\sqrt{x^T \Sigma x}$
  * Return: $x^T \mu$

Plots compare:

* In-sample frontiers (regular vs. shrinkage)
* Out-of-sample performance of these portfolios

---

##  Visualizations

* ETF time series plots
* Heatmap of missing data
* Log-likelihood and BIC progression
* Cross-validation MSE vs. $\gamma$ and $\lambda$
* Efficient frontiers and out-of-sample frontier comparison

---

##  Files

* `xMeanCovMissingMLE`: MLE for missing data
* `xFactorFitMLE`: MLE for factor model fitting
