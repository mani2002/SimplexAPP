import tkinter as tk
import numpy as np
'''
The input parameters of the simplex function are:

c: a numpy array representing the coefficients of the objective function to be maximized or minimized
A: a numpy array representing the coefficients of the constraints
b: a numpy array representing the right-hand side of the constraints
maximize: a boolean parameter indicating whether to maximize or minimize the objective function. The default is True (maximization).
sign: a list of strings representing the sign of the constraints. If sign is not provided, the function adds slack variables to convert the inequalities to equalities.
is_unrestricted: a list of boolean values indicating whether the variables are unrestricted. If is_unrestricted is not provided, the variables are assumed to be restricted.

The function returns a tuple of two numpy arrays:

x: the optimal solution
obj: the value of the objective function at the optimal solution
s: solution status indicating feasibility, unboundedness, or infeasibility
'''
def simplex(c, A, b, maximize=True,sign=None,is_unrestricted = None):
    m, n = A.shape
    assert b.shape[0] == m
    assert c.shape[0] == n
    art = np.zeros((m), dtype=int)
    l = 0
    h = 0
    if sign is None:
        # Adding slack variables to convert inequalities(<=) to equalities
        A = np.hstack([A, np.eye(m)])
        c = np.hstack([c, np.zeros(m)])
    else:
        # Adding surplus and artifitial variable to convert inequalities(>= or =) to equalities
        c = np.hstack([c, np.zeros(m)])
        I = np.eye(m)
        if maximize:
            M = -1e100
        else:
            M = 1e10
        for i,s in enumerate(sign):
            if s == ">=":
                art[i] = len(sign) - i + l
                c = np.hstack([c, M])
                I[i][i] = -1
                N =np.zeros(m)
                N[i] = 1
                N = N.reshape(m,1)
                I = np.hstack([I,N])
                l += 1
            elif s == "=":
                c[n+i] = M
        A = np.hstack([A, I])
 
    # Initializing the tableau
    tableau = np.vstack([np.hstack([A, b.reshape(-1, 1)]), np.hstack([c, 0])])
    
    # Adding non-negative variables for unrestricted in sign variable
    if is_unrestricted is not None:
        for i in range(len(is_unrestricted)):
            if is_unrestricted[i]:
                d = -np.copy(tableau[:,i])
                d.shape = (m+1,1)
                tableau = np.hstack((tableau[:,:i+1+h],d,tableau[:,i+1+h:]))
                d = -np.copy(A[:,i])
                d.shape = (m,1)
                A = np.hstack((A[:,:i+1+h],d,A[:,i+1+h:]))
                h += 1
    
    # Indices of basic and nonbasic variables
    basic_vars = np.arange(n, n+m)
    nonbasic_vars = np.arange(n)
    if sign is not None:
        basic_vars += art
        for i in range(n,n+m):
            if i not in basic_vars:
                nonbasic_vars = np.hstack([nonbasic_vars,[i]])
    if is_unrestricted is not None:
        basic_vars += np.ones(m,dtype=int)*h
        for i in range(n,n+m):
            if i not in basic_vars:
                nonbasic_vars = np.hstack([nonbasic_vars,[i]])
   

    # Executing the simplex algorithm
    while True:
        # Selecting the entering variable
        zj_cj = np.dot(tableau[-1, basic_vars],tableau[:-1,:-1]) - tableau[-1, np.arange(n+m+l+h)] 
        if maximize:
            entering_var = np.argmin(zj_cj)
        else:
            entering_var = np.argmax(zj_cj)
        # Verifying that the solution is unbounded
        if False not in (tableau[:-1,entering_var] <= 0):
            s = "Solution to the given problem is unbounded"
            return None,None,s
        
        # Determining whether an optimal solution has been achieved
        if zj_cj[entering_var] >= 0 and maximize:
             break
        elif zj_cj[entering_var] <= 0 and not maximize:
            break
        
        # Computing the ratios
        ratios = tableau[:-1, -1] / tableau[:-1, entering_var]
        ratios[tableau[:-1, entering_var] <= 0] = np.inf
        ratios[np.logical_not(np.isfinite(ratios))] = 1e16
        
        # Selecting the leaving variable
        leaving_var = np.argmin(ratios)
        
        # Updating the tableau
        pivot = tableau[leaving_var, entering_var]
        tableau[leaving_var, :] /= pivot
        for i in range(m+1):
            if i != leaving_var:
                tableau[i, :] -= tableau[i, entering_var] * tableau[leaving_var, :]  
        evi = np.where(nonbasic_vars == entering_var)[0][0]
        lv = basic_vars[leaving_var]
        basic_vars[leaving_var] = entering_var
        nonbasic_vars[evi] = lv
        
    # Extracting the Basic feasible solution and Optimal solution from the tableau
    x = np.zeros(n+h)
    for i in range(m):
        if basic_vars[i] < n+h:
            x[basic_vars[i]] = tableau[i, -1]
    obj = tableau[-1, -1]
    obj = -obj
    s = "Solution to the given problem is feasible" 
    
    # Verifying that the solution is infeasible
    for i in range(m):
        if sign[i] == "<=":
            if np.dot(A[i][0:len(x)], x) > b[i]:
                s = "Solution to the given problem is infeasible"
        elif sign[i] == ">=":
            if np.dot(A[i][0:len(x)], x) < b[i]:
                s = "Solution to the given problem is infeasible" 
        elif sign[i] == "=":
            if np.dot(A[i][0:len(x)], x) != b[i]:
                s = "Solution to the given problem is infeasible" 
                
    return x, obj,s

'''
The input fields include:

"Objective function coefficients": a single line entry field where the user enters the coefficients of the objective function.
"Constraint coefficients": a multi-line entry field where the user enters the coefficients of the constraint matrix, with one row per constraint.
"Constraint values": a single line entry field where the user enters the right-hand side values of the constraints.
"Enter Constraint types": a single line entry field where the user enters the type of each constraint, with "<=", ">=", or "=" symbols.
"unrestricted Variable or Not": a single line entry field where the user enters whether each variable is unrestricted or not, with "T" for unrestricted and "F" for restricted.
"Maximize": a checkbox that allows the user to select whether to maximize or minimize the objective function.

When the user clicks the "Solve" button, the program reads the values from the input fields, converts them to numpy arrays, and calls the simplex function to solve the linear programming problem. If a solution is found, the program displays it in a label on the GUI. If there is an error, the program displays an error message in the label.

'''

class SimplexGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Simplex Method Calculator")
        self.master.geometry("600x400")

        self.c_label = tk.Label(self.master, text="Objective function coefficients:")
        self.c_label.grid(row=0, column=0, padx=10, pady=10)
        self.c_entry = tk.Entry(self.master, width=30)
        self.c_entry.grid(row=0, column=1, padx=10, pady=10)

        self.A_label = tk.Label(self.master, text="Constraint coefficients (one row per constraint):")
        self.A_label.grid(row=1, column=0, padx=10, pady=10)
        self.A_entry = tk.Text(self.master, width=30, height=3)
        self.A_entry.grid(row=1, column=1, padx=10, pady=10)

        self.b_label = tk.Label(self.master, text="Constraint values:")
        self.b_label.grid(row=2, column=0, padx=10, pady=10)
        self.b_entry = tk.Entry(self.master, width=30)
        self.b_entry.grid(row=2, column=1, padx=10, pady=10)

        self.s_label = tk.Label(self.master, text="Enter Constraint types:")
        self.s_label.grid(row=3, column=0, padx=10, pady=10)
        self.s_entry = tk.Entry(self.master, width=30)
        self.s_entry.grid(row=3, column=1, padx=10, pady=10)
        self.s_entry.insert(0,"<=,>=,=")

        self.r_label = tk.Label(self.master, text="unrestricted Variable or Not(if yes:T, if not : F):")
        self.r_label.grid(row=4, column=0, padx=10, pady=10)
        self.r_entry = tk.Entry(self.master, width=30)
        self.r_entry.grid(row=4, column=1, padx=10, pady=10)
        self.r_entry.insert(0,"F,F,F")

        self.maximize_var = tk.BooleanVar()
        self.maximize_var.set(True)
        self.maximize_check = tk.Checkbutton(self.master, text="Maximize", variable=self.maximize_var)
        self.maximize_check.grid(row=5, column=0, padx=10, pady=10)

        self.solve_button = tk.Button(self.master, text="Solve", command=self.solve)
        self.solve_button.grid(row=5, column=1, padx=10, pady=10)

        self.result_label = tk.Label(self.master, text="")
        self.result_label.grid(row=6, column=0, columnspan=2, padx=10, pady=10)

    def solve(self):
        c = np.fromstring(self.c_entry.get(), sep=',')
        A = np.fromstring(self.A_entry.get("1.0", tk.END), sep=',').reshape(-1, len(c))
        b = np.fromstring(self.b_entry.get(), sep=',')
        signs = np.array((self.s_entry.get()).split(','))
        iu = np.array((self.r_entry.get()).split(','))
        is_unrestricteds = (iu == 'T')
        maximize = self.maximize_var.get()

        try:
            x, obj,s = simplex(c, A, b, maximize,sign=signs,is_unrestricted=is_unrestricteds)
            print(f"Solution: {x}\nObjective value: {obj}\n{s}")
            if x is not None:
                self.result_label.config(text=f"Solution: {x}\nObjective value: {obj}\n{s}")
            else:
                self.result_label.config(text=s)
        except Exception as e:
            self.result_label.config(text=str(e))

root = tk.Tk()
gui = SimplexGUI(root)
root.mainloop()
