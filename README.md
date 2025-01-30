# Lyapunov-Exponents
MATLAB script for solving augmented IVP with variational equation, calculating Lyapunov exponents and Kaplan—Yorke dimension.

## Table of Contents
- [Notes on Nonautonomous Systems](#notes-on-nonautonomous-systems)
- [The Variational Equation](#the-variational-equation)

## Notes on Nonautonomous Systems
Consider a nonautonomous IVP 

$$ \dot{\mathbf{z}}\left(t\right) = \mathbf{f}\left(t, \mathbf{z}\right), \quad \left.\mathbf{z}\left(t\right)\right|_{t=t_0}=\mathbf{z}_0, $$

where $\mathbf{z} \in \mathbb{R}^m, ~\mathbf{f}: \mathbb{R}\times\mathbb{R}^m \to \mathbb{R}^m$.

In order to compute the Lyapunov exponents, we first need to convert the nonautonomous system to an autonomous system as follows:

$$ z_{m+1} = t,$$

$$\dot{\mathbf{y}}\left(t\right) = \mathbf{g}\left(\mathbf{y}\right), \quad \left.\mathbf{y}\left(t\right)\right|_{t=t_0}=\mathbf{y}_0, $$

where $\mathbf{y} \in \mathbb{R}^{m+1}, ~\mathbf{g}: \mathbb{R}^{m+1} \to \mathbb{R}^{m+1}$:

$$\mathbf{y} = \begin{bmatrix}
			\mathbf{z}\\
			z_{m+1}	
		\end{bmatrix}, \quad  \mathbf{g}\left(\mathbf{y}\right) = \begin{bmatrix}
		\mathbf{f}\left(z_{m+1}	, \mathbf{z}\right)\\
		1
		\end{bmatrix}, \quad \mathbf{y}_0 = \begin{bmatrix}
		\mathbf{z}_0\\
		t_0
		\end{bmatrix}.$$
		
## The Variational Equation

Consider an autonomous IVP 

$$\dot{\mathbf{z}}\left(t\right) = \mathbf{f}\left(\mathbf{z}\right), \quad \left.\mathbf{z}\left(t\right)\right|\_{t=t\_0}=\mathbf{z}\_0, $$

where $\mathbf{z} \in \mathbb{R}^m, ~\mathbf{f}:\mathbb{R}^m \to \mathbb{R}^m$.

For this system the variational equation has the following form:

$$\dot{\boldsymbol{\updelta}}\left(t\right) = \mathbf{J}\_\mathbf{f}\left(\mathbf{z}\right) \boldsymbol{\updelta}\left(t\right), \quad \left.\boldsymbol{\updelta}\left(t\right)\right|\_{t=t_0} = \mathbf{I},$$
	
where $\mathbf{J}_\mathbf{f} \in \mathbb{R}^{m\times m}$ is Jacobian of the function $\mathbf{f}\left(\mathbf{z}\right)$ and $\boldsymbol{\updelta} \in \mathbb{R}^{m\times m}$ is variational matrix:
	
$$\mathbf{J}_\mathbf{f}\left(\mathbf{z}\right) = \frac{\text{d}\mathbf{f}}{\text{d}\mathbf{z}} = \begin{bmatrix}
			\boldsymbol{\nabla}^\text{T} f_1 
			\\
			\vdots
			\\
			\boldsymbol{\nabla}^\text{T} f_m
		\end{bmatrix}=
		\begin{bmatrix}
			\dfrac{\partial f_1}{\partial z_1} &\cdots& \dfrac{\partial f_1}{\partial z_m}\\
			\vdots & \ddots & \vdots \\
			\dfrac{\partial f_m}{\partial z_1} &\cdots& \dfrac{\partial f_m}{\partial z_m} 
		\end{bmatrix},$$	

$$\boldsymbol{\updelta}\left(t\right) = 
		\begin{bmatrix}
			\delta_{z_1 z_1}\left(t\right)&\cdots& \delta_{z_m z_1}\left(t\right)\\
			\vdots & \ddots & \vdots \\
			\delta_{z_1 z_m}\left(t\right) &\cdots& \delta_{z_m z_m}\left(t\right)
		\end{bmatrix}.$$
		
To find out what happens to the variations, you need to solve the variational equation and the system equation simultaneously. To do this, you work with a new augmented state vector of length $m + m^2$:

$$\mathbf{z}\_\ast = \begin{bmatrix}
			\mathbf{z}\\
			\boldsymbol{\updelta}\_{z\_1}\\
			\vdots\\
			\boldsymbol{\updelta}\_{z\_m}
		\end{bmatrix},$$
		
where $\forall m$

$$\boldsymbol{\updelta}\_{z_m} = \begin{bmatrix}
			\delta_{z_m z_1}\\
			\vdots\\
			\delta_{z_m z_m}
		\end{bmatrix}.$$
		
The augmented IVP is solved using [General Algorithm for the Explicit Runge—Kutta Method] (https://github.com/whydenyscry/General-algorithm-of-the-explicit-Runge-Kutta-method).