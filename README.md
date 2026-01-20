# Lyapunov Exponents
MATLAB script for solving augmented IVP with variational equation, calculating Lyapunov exponents and Kaplan—Yorke dimension.

## Table of Contents
- [Notes on Nonautonomous Systems](#notes-on-nonautonomous-systems)
- [The Variational Equation](#the-variational-equation)
- [The Lyapunov Exponents](#the-lyapunov-exponents)
- [The Kaplan—Yorke Dimension](#the-kaplanyorke-dimension)
- [Example](#example)

## Notes on Nonautonomous Systems
Consider a nonautonomous IVP 
```math
\dot{\mathbf{z}}\left(t\right) = \mathbf{f}\left(t, \mathbf{z}\right), \quad \left.\mathbf{z}\left(t\right)\right|_{t=t_0}=\mathbf{z}_0, 
```
where $`\mathbf{z}\left(t\right):\mathbb{R} \mapsto \mathbb{R}^m, ~\mathbf{f}\left(t, \mathbf{z}\right): \mathbb{R}\times\mathbb{R}^m \to \mathbb{R}^m`$.

In order to compute the Lyapunov exponents, we first need to convert the nonautonomous system to an autonomous system as follows:
```math
z_{m+1} = t, \\
\dot{\mathbf{y}}\left(t\right) = \mathbf{g}\left(\mathbf{y}\right), \quad \left.\mathbf{y}\left(t\right)\right|_{t=t_0}=\mathbf{y}_0
```
where $\mathbf{y}:\mathbb{R}\mapsto \mathbb{R}^{m+1}, ~\mathbf{g}\left(\mathbf{y}\right): \mathbb{R}^{m+1} \to \mathbb{R}^{m+1}$:
```math
\mathbf{y} = \begin{bmatrix}
			\mathbf{z}\\
			z_{m+1}	
		\end{bmatrix}, \quad  \mathbf{g}\left(\mathbf{y}\right) = \begin{bmatrix}
		\mathbf{f}\left(z_{m+1}	, \mathbf{z}\right)\\
		1
		\end{bmatrix}, \quad \mathbf{y}_0 = \begin{bmatrix}
		\mathbf{z}_0\\
		t_0
		\end{bmatrix}.
```
## The Variational Equation

Consider an autonomous IVP 
```math
\dot{\mathbf{z}}\left(t\right) = \mathbf{f}\left(\mathbf{z}\right), \quad \left.\mathbf{z}\left(t\right)\right|_{t=t_0}=\mathbf{z}_0, 
```
where $`\mathbf{z}\left(t\right):\mathbb{R} \mapsto \mathbb{R}^m, ~\mathbf{f}\left(\mathbf{z}\right): \mathbb{R}^m \to \mathbb{R}^m`$.

For this system the variational equation has the following form:
```math
\dot{\bm{\Phi}}\left(t\right) = \mathbf{J}_\mathbf{f}\left(\mathbf{z}\right) \bm{\Phi}\left(t\right), \quad \left.\bm{\Phi}\left(t\right)\right|_{t=t_0} = \mathbf{I}
```
	
where $`\mathbf{J}_\mathbf{f}\left(\mathbf{z}\right): mathbb{R} \mapsto \mathbb{R}^{m\times m}`$ is Jacobian of the function $`\mathbf{f}\left(\mathbf{z}\right)`$ 
```math
\mathbf{J}_\mathbf{f}\left(\mathbf{z}\right) = \frac{\text{d}\mathbf{f}}{\text{d}\mathbf{z}} = \begin{bmatrix}
			\boldsymbol{\nabla}^\top f_1 
			\\
			\vdots
			\\
			\boldsymbol{\nabla}^\top f_m
		\end{bmatrix}=
		\begin{bmatrix}
			\dfrac{\partial f_1}{\partial z_1} &\cdots& \dfrac{\partial f_1}{\partial z_m}\\
			\vdots & \ddots & \vdots \\
			\dfrac{\partial f_m}{\partial z_1} &\cdots& \dfrac{\partial f_m}{\partial z_m} 
		\end{bmatrix},
```
and $`\bm{\Phi}: \mathbb{R} \mapsto \mathbb{R}^{m\times m}`$ is variational matrix. To find out what happens to the variations, you need to solve the variational equation and the system equation simultaneously. To do this, you work with a new augmented state vector of length $m + m^2$:
```math
\mathbf{z}_\ast = \begin{bmatrix}
			\mathbf{z}\\
			\text{vec}\left(\bm{\Phi}\right)
		\end{bmatrix}
```
		
The augmented IVP is solved using [General Algorithm for the Explicit Runge—Kutta Method](https://github.com/whydenyscry/General-algorithm-of-the-explicit-Runge-Kutta-method). 

To test the algorithm, an example from [here](https://home.cs.colorado.edu/~lizb/chaos/variational-notes.pdf) was used. The results can be seen in the [ExampleOfUse.mlx](ExampleOfUse/ExampleOfUse.pdf) file.

## The Lyapunov Exponents

The calculation of the Lyapunov exponent was based on the QR decomposition method, the application of which can be viewed via the script [odeExplicitGeneralLE.m](Scripts/odeExplicitGeneralLE.m).

### Syntax
`[t, zsol, lyap_exp] = odeExplicitSolversLyapunovExponents(odefun, tspan, tau, incond)`
`[t, zsol, lyap_exp, dzdt_eval] = odeExplicitSolversLyapunovExponents(..., "Method", method_name, "SafeRegime", flag)`

### Input Arguments
- `odefun`: function handle that defines the **augmented** system of ODEs to be integrated. This function must compute the derivatives for both the original system state and the flattened variational matrix;
- `tspan`: interval of integration $[t_{start}, t_{end}]$, specified as a two-element vector;
- `tau`: fixed time discretization step $\tau$;
- `incond`: vector of initial conditions $\mathbf{z}_0$. This vector must consist of the initial state vector concatenated with the elements of the initial orthogonal matrix (usually the flattened identity matrix);
- `Method` (Name-Value Pair): string specifying the explicit Runge-Kutta method to be used. Available options:
  - `"RK3"`, `"RK4"` (default), `"RKB5"`, `"RKN5"`, `"RKB6"`, `"RKB7"`, `"RKCV8"`;
- `SafeRegime` (Name-Value Pair): logical flag (`true` or `false`). If set to `true`, the solver performs a check for `NaN` or `Inf` values at every integration step. Default is `false`.

### Output Arguments
- `t`: column vector of evaluation points used to perform the integration;
- `zsol`: matrix in which each row corresponds to the full solution vector (state + variations) at the value returned in the corresponding row of `t`;
- `lyap_exp`: matrix of Lyapunov exponents in which each row contains the spectrum $`[\lambda_1, \dots, \lambda_m]`$ calculated at the value returned in the corresponding row of `t`;
- `dzdt_eval`: (optional) matrix of the evaluated derivative vectors at each time step.

## The Kaplan—Yorke Dimension

Let the Lyapunov exponents be sorted in descending order $\lambda _{1}\geq \lambda _{2}\geq \dots \geq \lambda _{m}$, then

$$
D_\text{KY}=k+{\frac {\sum _{i=1}^{k}\lambda_{i}}{|\lambda_{k+1}|}},
$$

where for $k$
$$
\sum _{i=1}^{k}\lambda _{i}\geq 0, \quad \sum _{i=1}^{k + 1}\lambda _{i}<0.
$$

## Example

The Rössler Attractor in the [ExampleOfUse.mlx](ExampleOfUse/ExampleOfUse.pdf) was chosen as an example:
 
$$ 
\begin{cases}
			\frac{\mathrm{d}x}{\mathrm{d}t} =-y-z,\\
			\frac{\mathrm{d}y}{\mathrm{d}t} = x+\alpha y, \\
			\frac{\mathrm{d}z}{\mathrm{d}t} = \beta+z\left(x-\varsigma\right),
		\end{cases}
$$

$$ 
\begin{bmatrix}
			\alpha\\
			\beta\\
			\varsigma
		\end{bmatrix}=\begin{bmatrix}
		0.2\\
		0.2\\
		5.7
		\end{bmatrix}.
$$


<p align="center">
  <img src="ExampleOfUse/images/The_Rossler_Attractor.png"/>
</p>

<p align="center">
  <img src="ExampleOfUse/images/time_series.png"/>
</p>

<p align="center">
  <img src="ExampleOfUse/images/LEs_plot.png"/>
</p>
