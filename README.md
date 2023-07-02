# flexible_lp
A pure python3 library to solve linear programs using, for example, fractions.

The main interface is the `linprog` function in simplex.py.
Here is its signature and docstring:
```
def linprog(c, d = 0, A_g = None, b_g = None, A_e = None, b_e = None, A_l = None, b_l = None, maximize = False, verbose = False, value_map = lambda x: x):
	'''
	Parameters:
		The first parameters define a linear program:
			objective function: c^T x + d
			constraints: A_g x >= b_g, A_e x = b_e, A_l x <= b_l, x >= 0.
		the A's (if not None) should all be lists of lists of length len(c)
		the b's (if not None) should all be lists with the length of their corresponding A

		maximize: will minimize by default, but will maximize if maximize set to True.
		verbose: set True to print more detailed information during the solve.
		value_map: will apply a map to all numbers involved before solving.
					For example, you might want to solve over rationals but enter the program in integers,
					in which case you could import fractions and set "value_map = lambda x: Fraction(x)",
					or just "value_map = Fraction".

	Return:
		usual solution :: (optimal_value, x_vector_obtaining_this_value)
		infeasible problem :: (None, "infeasible program, bad row encountered in phase 1.")
		unbounded problem :: (float("inf"} or float("-inf"), x_vector_obtaining_this_value)

	Note:
		All the numbers involved, from c, d, etc, can be any class that supports the following:
			-  *,  /, +, unary -, str [for Tableau pivotting and printing]
			- /, > (with own type and 0), < (with own type and 0),  [for the simplex algorithm]
	
	'''
```

# Examples
Suppose we had the following lp
$$3x + 5y;\ x + y \leq 10,\ x,y\geq0.$$

Here is a complete python3 file to solve using fractions:
```
from simplex import linprog
from fractions import Fraction

opt_val, opt_vec = linprog(c = [3,5], 
```
