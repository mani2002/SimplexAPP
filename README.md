# Simplex Method Calculator

## Overview
This Python code is an implementation of the simplex method to solve linear programming problems. The calculator is built using the Python programming language and utilizes the numpy and tkinter libraries for numerical computations and graphical user interface (GUI) design, respectively.

## Features
The calculator allows the user to input a linear programming problem in standard form and it will then solve the problem using the simplex method. Specifically, the calculator will perform the following steps:
- Convert the problem into tableau form.
- Identify the entering and leaving variables using the minimum ratio test.
- Update the tableau by pivoting on the entering and leaving variables.
- Repeat steps 2 and 3 until an optimal solution is found.

The calculator will display the optimal solution, Objective value and nature of the solution(feasible, infeasible and Unbounded).

## Dependencies
This calculator requires the following libraries to be installed:
- `numpy`: for numerical computations.
- `tkinter`: for GUI design.

## How to use
1. Clone this repository to your local machine.
2. Open the command prompt and navigate to the cloned repository directory.
3. Run the following command to install the required libraries: 
   - `pip install numpy tkinter`
4. Run the following command to start the calculator:
   - `py SimplexAPP.py`
5. Enter the linear programming problem.
6. Click on the `Solve` button to solve the problem.
7. The optimal solution and objective function value and solution nature will be displayed.

## Example 
 **Maximize Z = 6x1 + 4x2  
subject to the constraints  
(i) 2x1 + 3x2 ≤ 30,  
(ii) 3x1 + 2x2 ≤ 24,  
(iii) x1 + x2 ≥ 3  
and x1, x2 ≥ 0**  

We can enter this problem into the calculator. After clicking on the Solve button, the calculator will display the following output:  

![example](https://github.com/mani2002/SimplexAPP/blob/main/Example.png)
