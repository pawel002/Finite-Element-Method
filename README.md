## Solving electromagnetic potential DE using Finite Elements Method
Equation on domain (0, 3):
<p align="center">
<img src="https://github.com/pawel002/fem-diff-equation/blob/main/images/equation.png"
     alt="Equation."
     style="display: inline-block; margin: 0 auto; max-width: 300px">
</p>

Numerical solution for n=1000:
<p align="center">
<img src="https://github.com/pawel002/fem-diff-equation/blob/main/images/graph.png"
      style="display: inline-block; margin: 0 auto; max-width: 300px">
</p>

## Description
This repository contains a project made for the Differential Equations course at AGH University in Kraków. The derivation of the weak formulation for the equation can be seen [here](https://github.com/pawel002/fem-diff-equation/blob/main/images/FEM.pdf) in Polish. To compile the code, run:
```
g++ -I .\Eigen .\main.cpp -o main
```
Integration method: [Gauss-Legendre Quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature). <br />
Sparse matrix operations: [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page). <br />
Plotting using short matplotlib and python programm. <br />
In linux version, mathgl is used for plotting. ( -lmgl needed when compiling!)  <br />
Eigen is included in the files.
