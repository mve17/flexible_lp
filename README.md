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
		usual problem :: (optimal_value, x_vector_obtaining_this_value)
		infeasible problem :: (None, "infeasible program, bad row encountered in phase 1.")
		unbounded problem :: (float("inf"} or float("-inf"), x_vector_obtaining_this_value)

	Note:
		All the numbers involved, from c, d, etc, can be any class that supports the following:
			-  *,  /, +, unary -, str [for Tableau pivotting and printing]
			- /, > (with own type and 0), < (with own type and 0),  [for the simplex algorithm]
	
	'''
```
# Example 0: Solving an LP
Suppose we had the following lp
$$\text{maximize}\ 3x + 5y;\ x + y \leq 10,\ x,y\geq0.$$
Its optimal value is 50, with $x = 0, y = 10$.
Here is a python3 file to solve it:
```
from simplex import linprog
from fractions import Fraction

opt_val, opt_vec = linprog(c = [3,5], A_l = [[1, 1]], b_l = [10], maximize = True)
print(opt_val, opt_vec)
```
This outputs the following:
```
>>> 50.0 [0, 10.0]
```
Note it will solve using floating point numbers by default. The next example shows how to solve exactly with fractions.

# Example 1: Solving an LP with Fractions
Suppose we have the same lp
$$\text{maximize}\ 3x + 5y;\ x + y \leq 10,\ x,y\geq0.$$
Here is a python3 file to solve using fractions. 
The `value_map` parameter converts the supplied integers to the Fraction type before solving:
```
from simplex import linprog
from fractions import Fraction

opt_val, opt_vec = linprog(c = [3,5], A_l = [[1, 1]], b_l = [10], maximize = True, value_map = Fraction)
print(opt_val, opt_vec)
```
This outputs the following:
```
>>> 50 [0, Fraction(10, 1)]
```
# Example 2: Solving an LP with Unbounded Value (and with Fractions)
$$\text{minimize}\ -2x + 7y + z$$

$$x + y + z \geq 10,$$

$$5x+3y+1\geq 9,$$

$$x,y,z\geq0.$$
```
from simplex import linprog
from fractions import Fraction
print(linprog(c = [-2,7,1], A_g = [[5,1,1], [1,3,1]], b_g = [5,9], value_map = Fraction))
```

```
>>> (-inf, [inf, 0, 0])
```

