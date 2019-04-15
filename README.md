# modeFinder

## Mode estimation via Bernstein polynomials

The R function is based on the draft of:

Bowen Liu and Sujit K. Ghosh, *”On empirical estimation of mode based on weakly dependent samples”*, in preparation.


## Simple empirical mode estimation

It is natural to estimate the population distribution function by the empirical distribution function, $F_n$ given by
$$
F_n = \frac{1}{n} \sum^n_{i = 1} I \{X_i \leq x\}
$$
and for any $X^*$ such that $X^{(i)} < X^* < X^{(i + 1)}$, $i = 1, 2, \dots, n$, the slope of the empirical distribution function can be approximated by 
$$
\frac{ F_n (X^{(i + 1)}) - F_n (X^{(i)})}{X^{(i + 1)} - X^{(i)}}.
$$
Consequently, we can give a simple mode estimator, $u^{(\hat{k})}$ based on the empirical distribution function by
$$
\hat{k} = \underset{0 \leq i \leq m^* - 1}{argmax} \frac{ F_n (u^{(i + 1)}) - F_n (u^{(i)})}{u^{(i + 1)} - u^{(i)}}.
$$
where $\{u^{(1)}, u^{(2)}, \dots, u^{(m^*)}\} = \{ \frac{1}{m^*}, \frac{2}{m^*}, \dots, \frac{m^*}{m^*}\}$. However, the simple empirical mode estimation does not collaborate with the possible smoothness characteristics of the underlying distribution function. 






### Note: 

macOS requires Xcode developer tools otherwise *xcrun: error* flags. Run the following: *xcode-select --install*.

