import numpy as np
import sympy as sp
from IPython.display import display, Math
import sys

def write_equations_to_latex_file(filename, equations,bmat):
    # Open the file in write mode
    with open(filename, 'w') as file:
        # Write the LaTeX header
        file.write("\\documentclass{article}\n")
        file.write("\\usepackage{amsmath}\n")
        file.write("\\begin{document}\n")
        file.write("For the following b and corresponding Q:")
        file.write(sympy_matrix_to_latex(bmat))
        # Write the equations in aligned LaTeX form
        file.write("\\begin{align*}\n")
        for equation in equations:
            latex_equation = sp.latex(equation).replace("=", "&=")
            file.write(f"\t{latex_equation} \\\\\n")
        file.write("\\end{align*}\n")

        # Write the LaTeX footer
        file.write("\\end{document}\n")

#assume we're working in k/4
print("how many columns of b do you want to specify?")
columns = int(input())
print("What is the m we're working in for osp(k|2m)?")
m = int(input())
b = np.empty((2*m,columns),dtype = np.complex_)
for i in range(2*m):
    print("please write out the values of row " + str(i+1) + " with a space in between each value.")
    text = input()
    templist = text.split()
    #for i in range(columns):
    #    templist[i] = complex(templist[i])
    b[i] = templist
print(b)

Im = np.empty((m,m),dtype = np.complex_)
for i in range(m):
    Im[i,i] = 1
    
Im = sp.Matrix(Im)
#print(Im)
matJ = sp.Matrix(sp.BlockMatrix([[sp.Matrix.zeros(m),Im],[-Im,sp.Matrix.zeros(m)]]))
#print(matJ)
matJ = np.array(matJ.tolist())
print(matJ)

jb = np.matmul(matJ,b)
#print(jb)
btjb = np.matmul(np.transpose(-b),jb)
#
if np.any(btjb):
    print("-bTJb not a zero matrix, thus given b and corresponding Q is not nilpotent")
    print(btjb)
    sys.exit("Program stopped")
#check if Q is correct
jb = np.matmul(-b,np.transpose(b))
btjb = np.matmul(jb,matJ)
if np.any(btjb):
    print("-bbTJ not a zero matrix, thus given b and corresponding Q is not nilpotent")
    print(btjb)
    sys.exit("Program stopped")

print("The given b and corresponding Q is indeed in the nilpotent")
#check if Q is correct

#Actually, it seems that using the number of rows we specify as k is an easily extrapolated difference.
#now lets create alpha_1, a columns x columns matrix
alpha1_symbols = sp.symbols('m1:%d' % ((columns*columns-columns)/2+1))
alpha_1 = sp.Matrix.zeros(columns,columns)
#print(alpha1_symbols[1])
index = 0
for i in range(columns):
    for j in range(i+1,columns):
        alpha_1[i,j] = alpha1_symbols[index]
        alpha_1[j,i] = -alpha1_symbols[index]
        index = index + 1
#print(alpha_1)

matJ_var = sp.Matrix(matJ)

b_var = sp.Matrix(b)

#create the general sp4 matrix

a_symbols = sp.symbols('a1:%d' % 5)
b_symbols = sp.symbols('b1:%d' % 4)
c_symbols = sp.symbols('c1:%d' % 4)

a1, a2, a3, a4 = a_symbols
b1, b2, b3 = b_symbols
c1, c2, c3 = c_symbols

a_matrix = sp.Matrix([[a1,a2],[a3,a4]])
a2_matrix = sp.Matrix([[-a1,-a3],[-a2,-a4]])
b_matrix = sp.Matrix([[b1,b2],[b2,b3]])
c_matrix = sp.Matrix([[c1,c2],[c2,c3]])

sp_matrix = sp.Matrix(sp.BlockMatrix([[a_matrix,b_matrix],[c_matrix,a2_matrix]]))
 
print(sp_matrix)

#lets form the second equality first: αδTJ = δTJγ 
print("We will first look at the equality so5 * bT * J = bT * J * sp4, the left side yields:")
bTJ = b_var.T * matJ_var
so5bTJ = alpha_1 * bTJ
print("so5bTJ:")
print(so5bTJ) 

bTJsp = bTJ * sp_matrix
print("bTJsp:")
print(bTJsp)

def sympy_matrix_to_latex(matrix):
    # Convert the SymPy matrix to its LaTeX representation
    latex_str = sp.latex(matrix)

    # Enclose the LaTeX string with the necessary LaTeX matrix environment
    latex_matrix = "\\begin{bmatrix}\n" + latex_str + "\n\\end{bmatrix}"

    return latex_matrix

#look at the latex print of the matrix
print(sympy_matrix_to_latex(bTJsp))

#form the system of equations

equations = []
for i in range(columns):
    for j in range(4):
        # Pair corresponding elements from both matrices
        element1 = so5bTJ[i, j]
        element2 = bTJsp[i, j]

        # Form the equation by setting both elements equal to each other
        equation = sp.Eq(element1, element2)

        # Append the equation to the list
        equations.append(equation)
simpleeq1 = [sp.simplify(eq) for eq in equations]
simpleeq1 = [sp.simplify(eq) for eq in simpleeq1]

#equations
print("\\begin{align*}")
for equation in simpleeq1:
    latex_equation = sp.latex(equation).replace("=", "&=")
    print(f"\t{latex_equation} \\\\")
print("\\end{align*}")

#second equation
spb = sp_matrix * b_var
balpha = b_var * alpha_1
 
for i in range(4):
    for j in range(columns):
        # Pair corresponding elements from both matrices
        element1 = spb[i, j]
        element2 = balpha[i, j]

        # Form the equation by setting both elements equal to each other
        equation = sp.Eq(element1, element2)

        # Append the equation to the list
        equations.append(equation)
simpleeq2 = [sp.simplify(eq) for eq in equations]
simpleeq2 = [sp.simplify(eq) for eq in simpleeq2]
eq_set = set(simpleeq2)
simpleeq2 = list(eq_set)
#equations 2
print("for the following b:")
print(sympy_matrix_to_latex(b_var))

solutions = sp.solve(simpleeq2, dict=True)
'''
# Isolate each individual 'm' prefix variable in the equations
equations_with_m_prefix = []
for equation in simpleeq2:
    m_variables = [var for var in equation.free_symbols if str(var).startswith('m')]
    isolated_eq = equation
    for var in m_variables:
        isolated_eq = sp.solve(isolated_eq, var)[0]
    equations_with_m_prefix.append(isolated_eq)
'''
print("\\begin{align*}")
for equation in simpleeq2:
    latex_equation = sp.latex(equation).replace("=", "&=")
    print(f"\t{latex_equation} \\\\")
print("\\end{align*}")
print(sympy_matrix_to_latex(alpha_1))
print(sympy_matrix_to_latex(sp_matrix))
#matrix rep test
equations_moved_left = [eq.lhs - eq.rhs for eq in simpleeq2]
# Simplify the equations to collect like terms and express in standard form
simplified_equations = [eq.simplify() for eq in equations_moved_left]
# Get all the unique variables from the simplified equations
variables = set().union(*[eq.free_symbols for eq in simplified_equations])

# Create the ordered list of variables with the desired order
variable_info = [(str(var)[0], int(str(var)[1:])) for var in variables]

# Sort variables based on prefix and then suffix numbers
ordered_variables = []
seen_variables = set()

for prefix in ['m', 'a', 'b', 'c']:
    prefix_variables = [var for var, (p, num) in zip(variables, variable_info) if p == prefix]
    sorted_prefix_variables = sorted(prefix_variables, key=lambda var: next(num for p, num in variable_info if p == prefix and var == sp.symbols(p + str(num))))
    
    for var in sorted_prefix_variables:
        if var not in seen_variables:
            ordered_variables.append(var)
            seen_variables.add(var)
#
num_equations = len(simplified_equations)
num_variables = len(ordered_variables)
coeff_matrix = sp.zeros(num_equations, num_variables)

for i, equation in enumerate(simplified_equations):
    for j, variable in enumerate(ordered_variables):
        coeff_matrix[i, j] = equation.coeff(variable)

print(sympy_matrix_to_latex(coeff_matrix))

# Create a vertical matrix of ordered variables
ordered_variables_matrix = sp.Matrix(ordered_variables)
print("\nOrdered Variables Matrix:")
print(sympy_matrix_to_latex(ordered_variables_matrix))

# Perform row reduction (Gaussian elimination) on the coefficient matrix
rref_matrix, pivot_columns = coeff_matrix.rref()


# Output the row-echelon form of the coefficient matrix
print("Row-Echelon Form:")
print(sympy_matrix_to_latex(rref_matrix))
