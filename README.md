# Schrodinger-Equation
Using verlet method to solve the Schrodinger for an arbitrary potential and bisection method to solve for eigen-energies

## Background

### Verlet Method

The Verlet method is a method for numerically solving a second order differential equation (typically of a specific form, which is what is used here).

It achieves a greater degree of numerical accuracy than Euler's method typically does, and is relatively simple in comparison to Runge-Kutta.

Let us say that we have a differential equation of the form:

```math
\ddot y = f(y(t))
```
Or, in the case here:

```math
\frac{d^2}{dx^2} y = f(y(x))
```

Then, for example, with initial conditions $y(0) = c_1$ and $y'(0) = c_2$, Verlet can be intialized with:

```math
y(\Delta x) = y(0) + y'(0) \Delta x + \dfrac 1 2 f(y(0)) \Delta x^2
```

And then, after initialization, verlet can run with:
```math
y(x + \Delta x) = 2y(x) - y(x - \Delta x) + f(y(x)) \Delta x^2
```
