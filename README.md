# TwoStageOptimalControl.jl
[![Build Status](https://github.com/michaelfreiberger/TwoStageOptimalControl.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/michaelfreiberger/TwoStageOptimalControl.jl/actions/workflows/ci.yml?query=branch%3Amaster)


This is a Julia-package which allows for the solution of two-stage optimal control problems and vintage-structured optimal control problems.

# Installation


The package is currently in submission to the Julia repository. Until acceptance, you are able to directly install the package from my Github-Repo by using the following commands in Julia:

> `using Pkg`
>
> `Pkg.add(url="https://github.com/michaelfreiberger/TwoStageOptimalControl.jl")`

# Documentation

For information on how to use the toolbox and some examples, please take a look into the documentation (ADD LINK). 

# Problem class

This package is designed to solve problems of the form

$$ 
\max_{C_1: [0,T]\to\mathbb{R}^{n_1},\quad C_2:[0,T]^2\to\mathbb{R}^{n_2}}\int_{0}^{T} \exp(-\rho t)\Big[F_1(t,C_1(t),S_1(t)) + Q(t) \Big]dt +
$$

$$
\hspace{4cm} + D_1(S(T)) + \int_0^T D_2(S_2(T,s))ds
$$

$$
s.t. \qquad\dot{S_1}(t) = g_1(t,C_1(t),S_1(t)) ,\hspace{2cm} S_1(0) = S_1^0
$$

$$
\hspace{3cm} \frac{d}{dt} S_2(t,s) = g_2(t,s,C_2(t,s),S_2(t,s)) ,\qquad S_2(s,s) = h(s,C_1(s),S_1(s))
$$

$$
\qquad \qquad Q(t) = \int_0^t F_2(t,s,C_2(t,s),S_2(t,s)) ds
$$

$$
\underline{C_1}(i) \leq C_{1,i}(t) \leq \overline{C_1}(i) \qquad\forall i = 1,\ldots,n_1 \quad,\quad\forall t\in [0,T]
$$

$$
\underline{C_2}(i) \leq C_{2,i}(t,s) \leq \overline{C_2}(i) \qquad\forall i = 1,\ldots,n_2 \quad,\quad\forall s\in [0,T]\quad,\quad\forall t\in [s,T]
$$


## **Variable and function explanations**

The model consist of three different types of variables:

* concentated state and control variables (depending on time $t$)
* distributed state and control variables (depending on time $t$ and vintage $s$)
* aggregated state variables (depending on time $t$)

### Concentated variables

The concentated variables consist of the $n_1$-dimensional control variable $C_1$ and the $m_1$-dimensional state variable $S_1$. The initial value of $S_1$ is given by $S_1^0$ and its dynamics are affected by $C_1$. Furthermore $C_1$ also affects the initial values of the distributed variables (see below).

### Distributed variables

Analogous to the concentated variables, the distributed variables consist of the $n_2$-dimensional control variable $C_2$ and the $m_2$-dimensional state variable $S_2$. The strictly speaking partial differential equation can be solved along the characteristic lines with fixed vintage $s=\overline{s}$. While the dynamics of $S_2$ are again affected by the distributed controls $C_2$, the initial values of $S_2$ at the beginning of each vintage $s$ are determined by the value of the concentated state $S_1$ at time $s$ and the concentated control $C_1$.

### Aggregated variables

Finally we also have an aggregated variable $Q$ which aggregates at each point in time $t$ a functional form $F_2$ over all vintages which startet before $t$. 

## **Algorithm principles**

The algorithm uses a gradient based approach following a sequence of steps:

1. Start with a given guess for the optimal solution of the control variables.
2. Calculate the corresponding profiles of the state variables and co-state variables (based on the maximum principle)
3. Calculate the gradient based on the Hamiltonian.
4. Adjust the guess for the control in the direction of the gradient.
5. Find the optimal adjustment step in the direction of the gradient.
6. Define new currently best solution.
7. Iterate from step 2 until no improvement can be found anymore.

## **Two-stage optimal control problems with random switching time**

This package is based on a project about the solution of two-stage optimal control problems with random switching time. This problem class takes the following form:

$$ 
\max_{C_1: [0,T]\to\mathbb{R}^{n_1}}\mathbb{E}_{\tau}\int_{0}^{\tau} \exp(-\rho t) F_1(t,C_1(t),S_1(t)) dt + \exp(-\rho\tau)V^*(\tau,C_1(\tau),S_1(\tau))
$$

$$
s.t. \qquad\dot{S_1}(t) = g_1(t,C_1(t),S_1(t)) ,\hspace{2cm} S_1(0) = S_1^0
$$

$$
\underline{C_1}(i) \leq C_{1,i}(t) \leq \overline{C_1}(i) \qquad\forall i = 1,\ldots,n_1
$$

with $V^*(\tau,S_1(\tau))$ being the value-function of the second stage problem

$$
V^*(\tau,S_1(\tau)) = \max_{C_2:[\tau,T]\to\mathbb{R}^{n_2}}\int_{\tau}^T \exp(-\rho(t-\tau))F_2(t,\tau,C_2(t),S_2(t)) dt
$$

$$
\hspace{3cm} \dot{S_2}(t) = g_2(t,\tau,C_2(t),S_2(t)) ,\qquad S_2(\tau) = h(\tau,C_1(\tau),S_1(\tau))
$$

$$
\underline{C_2}(i) \leq C_{2,i}(t) \leq \overline{C_2}(i) \qquad\forall i = 1,\ldots,n_2 \quad,\quad\forall t\in [\tau,T]
$$

and $\tau$ being the random switching time between the two stages with the stochastic properties being best described by the $Z(t):= \mathbb{P}\left[\tau > t\right]$, the probability of surviving in the first stage. The hazard rate $\eta$ of $Z$ thereby can depend on time $t$ the state variables $S_1$ and the control variablee $C_1$.

$$
\dot{Z}(t) = -\eta(t,C_1,S_1)\cdot Z(t)
$$

Overall this problem type can be transformed into a vintage-structed optimal control problem presented above. For details see the work of [Wrzaczek, Kuhn and Frankovic (2020)](https://link.springer.com/article/10.1007/s10957-019-01598-5)

# Planned extensions of the toolbox

### 1. Generalise to the age-structured optimal control problems with initial distribution of state variables, intial and boundary controls, and more generalised version of aggregated variables.

### 2. Allow for more general/complicated feasible regions for the control variables.

### 3. Include the option to define equilibrium conditions, which the solution has to fulfill.