### Fused LASSO using alternating direction method of multipliers (ADMM)

Few R packages exist with flexible implementations of the fused LASSO. Here, I showcase my implementation of the fused LASSO for discovering differential partial correlations. For example, consider the gene expression of treatment $(1)$ versus control $(0)$ groups. Here, we have a response vector
$Y=[Y_1,Y_0]$. The design matrix is $X =  \text{BlockDiag}(X_1, X_0) $. In this setting, we essentially run two separate regressions, as the beta coefficients are paired yet unconstrained. The fused LASSO imposes both a penalty on the beta coefficients, and the differences between beta coefficients of corresponding variables. 

# Objective function
In many biological applications, we are interested in the association of a response variable $Y \in \mathbb{R}^{N \times 1}$, with respect to some predictor variables $X \in \mathbb{R}^{N \times p}$, across varying conditions. In the simple 2-condition scenario, we have the following model: 
$$Y = X \beta + \varepsilon $$
$$ \varepsilon \sim \text{N}(0, \sigma^2) $$
$$ Y = [Y_1, Y_0]^T $$
$$ X = \text{BlockDiag}(X_1, X_0) $$
\\
To choose relevant predictor variables, we can impose an L1 penalty on the beta coefficient vector, $\beta$. Additionally, when we expect the relationship between the predictors and response to be mostly the same across groups, it is useful to impose a penalty on the beta coefficients of paired coefficients, $\beta_j$ and $\beta_{j+p}$, between treatment and control groups, respectively. We then minimize: 
$$ \frac{1}{2N} || Y - X \beta ||_2^2 + \lambda_1 || \beta ||_1 + \lambda_2 || D \beta ||_1 $$
However, because there are multiple nondifferentiable L1 penalty factors, the objective function cannot be optimized using standard gradient based methods. To efficiently solve the optimization, we introduce two auxiliary variables, $z$ and $d$, in place of $\beta$ and $D \beta$:
$$z = \beta $$
$$d = D \beta$$ 
Our new objective function becomes: 
$$ \frac{1}{2N} ||Y - X \beta ||_2^2 + \lambda_1 ||z||_1 + \lambda_2 ||d||_1 ; $$
$$ \text{subject to } \beta - z = 0, \ \  D \beta -d = 0$$
As this is a constrained optimization problem, we use the augmented Lagrangian as our objective function (See appendix for details on the derivation):

$$\begin{aligned}
    \mathcal{L} (\beta, z, d, u, w) =& \frac{1}{2N} ||Y - X \beta ||_2^2 + \lambda_1 ||z||_1 + \lambda_2 ||d||_1 \\
    & \quad  + \frac{\rho_1}{2} || \beta - z + u ||_2^2 + \frac{\rho_2}{2} || D \beta - d + w ||_2^2
\end{aligned}$$

\section*{Optimization}

We can now differentiate with respect to $\beta$ and then update $z$ and $d$ at each step. This procedure resembles a block coordinate descent algorithm, where we alternate between optimizing over $\beta$, $(z, d)$, and updating the dual variables $(u, w)$. The key is that we have removed the non-differentiable terms with respect to $\beta$ and exchanged them for auxiliary variables so that the function is differentiable over $\beta$. \\\\
\textit{Updating $\beta$} --- We differentiate the augmented Lagrangian with respect to $\beta$ and set it equal to zero. This results in the following linear system: 
$$ 
\begin{aligned}
    \Big( \frac{1}{N} X^T X + \rho_1 I + \rho_2 D^T D \Big) =& \frac{1}{N}X^T Y + \rho_1 (z^{(k)}-u^{(k)} ) \\
    & + \rho_2 (z^{(k)} - u^{(k)} ) + \rho_2 D^T (d^{(k)} - w^{(k)} )
\end{aligned}$$
Which we solve by Cholesky factorization of the LHS using the R function \texttt{chol}. \\\\
\textit{Update $z$ and $d$} --- At the $k$-th step, $z^{(k+1)}$ and $d^{(k+1)}$ are updated as follows: 
$$ z^{(k+1)} = \text{S} \Big( \beta^{(k+1)} + u^{(k)} , \frac{\lambda_1}{\rho_1}\Big) $$
$$ d^{(k+1)} = \text{S} \Big( D \beta^{(k+1)} + w^{(k)}, \frac{\lambda_2}{\rho_2}\Big) $$
Where $\text{S}(x, t)$ is the soft threshold operator: 
$$ \text{S} (x, t) = \text{sign} (x) \cdot \max ( |x| - t, 0 ) $$
The threshold, $t_j = \lambda_j / \rho_j $ is based on the Karush-Kuhn-Tucker conditions of the augmented Lagrangian (see appendix for details). 

\textit{Update $u$ and $w$} --- $u$ and $w$ are updated as follows: 
$$ u^{(k+1)} = u^{(k)} + \beta^{(k+1)} - z^{(k+1)} $$
$$ w^{(k+1)} = w^{(k)} + D \beta^{(k+1)} - d^{(k+1)} $$
We iterate through these steps until the objective function decreases by less than $10^{-5}$ for ten consecutive iterations. 
\end{document}
