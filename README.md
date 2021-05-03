# Fast Poisson Solver and Integral Equations for 2D Elliptical PDEs

## Project Aims

We aim to look at a faster approach for solving the typical finite difference system of elliptical two dimensional PDEs. The method uses the Fast Fourier Transform (FFT) algorithm to quickly solve the finite difference system. We will also look into if this method can be applied to non-uniform and/or non-square meshes. Following looking into this approach we will see how the results of the former method compare to a approach using integral equations. The two PDEs we will apply the two aforementioned methods to are the Laplace and Poisson equations.

**Laplace's Equation:**

![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Baligned%7D%20-%5Cnabla%5E2%20%5Cpsi%20%26%3D%200%20%5Cquad%20%5Ctext%7Bor%7D%20%5Cquad%20%5Ctext%7Binside%7D%20%5Cquad%20%5COmega%20%5Cin%20%5Cmathbb%7BR%7D%5E2%20%5C%5C%20%5Cpsi%28x%2C%20y%29%20%26%3D%20g%28x%2C%20y%29%20%5Cquad%20%5Ctext%7Bon%7D%20%5Cquad%20%5Cpartial%20%5COmega%20%5Cend%7Baligned%7D)

**Poisson's Equation:**

![equation](https://latex.codecogs.com/gif.latex?%5Cbegin%7Baligned%7D%20%26-%5Cnabla%5E2%20%5Cvarphi%20%3D%20f%28x%2C%20y%29%20%5Cquad%20%5Ctext%7Binside%7D%20%5Cquad%20%5COmega%20%5Cin%20%5Cmathbb%7BR%7D%5E2%20%5C%5C%20%26%5Cvarphi%28x%2C%20y%29%20%3D%20g%28x%2C%20y%29%20%5Cquad%20%5Ctext%7Bon%7D%20%5Cquad%20%5Cpartial%20%5COmega%20%5Cend%7Baligned%7D)

We will look at a multitude of examples and compare the accuracy between the two methods.

## About this rep

This repo contains the report and presentation for our APPM 4660: Numerical Anaylsis II final project. The code for the project can be found under the matlab directory.
