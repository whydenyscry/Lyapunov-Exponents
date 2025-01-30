# Lyapunov-Exponents
MATLAB script for solving augmented IVP with variational equation, calculating Lyapunov exponents and Kaplan—Yorke dimension.

## Table of Contents
- [Notes on Nonautonomous Systems](#notes-on-nonautonomous-systems)

## Notes on Nonautonomous Systems
Consider a nonautonomous IVP 

$$ \dot{\mathbf{z}}\left(t\right) = \mathbf{f}\left(t, \mathbf{z}\right), \quad \left.\mathbf{z}\left(t\right)\right|_{t=t_0}=\mathbf{z}_0, $$

where $\mathbf{z} \in \mathbb{R}^m, ~\mathbf{f}: \mathbb{R}\times\mathbb{R}^m \to \mathbb{R}^m$.

In order to compute the Lyapunov exponents, we first need to convert the nonautonomous system to an autonomous system as follows:

$$ z_{m+1} = t,\\ \dot{\mathbf{y}}\left(t\right) = \mathbf{g}\left(\mathbf{y}\right), \quad \left.\mathbf{y}\left(t\right)\right|_{t=t_0}=\mathbf{y}_0, $$