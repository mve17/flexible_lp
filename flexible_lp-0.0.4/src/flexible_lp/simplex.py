'''
Matthew Ellison (m.gr at dartmouth.edu)
Summer 2023

Solve linear programs exactly over the rationals! (or if you're creative, some other number class you define)

Implements the two-phase simplex algorithm using Tableau objects. My reference was Vaserstein's Introduction to Linear Programming.
I have not worked to prevent 'cycling', but my understanding is that this doesn't really arise in practice.

You are free, with attribution, to use this code in any way.
'''
from .tableau import Tableau

def phase_1(tableau, verbose = False):
	'''
	This is the first phase of the two-phase simplex method. 
	It attempts to perform pivots so that the basic solution (all top variables set to 0) is feasible.

	returns:
		("success", []) when the tableau is successfully modified to to have feasible basic solution.
		("bad row", [bad_row_idx]) when the lp is infeasible (i.e. the constraints cannot be satisfied).
	'''
	if verbose:
		print("----\nStarting Phase 1 of Simplex. The goal of this phase is to find a vertex in the feasible region.\n----")
	while True:
		#check if already feasible (i.e. b>=0, where the tableau is [[A,b],[c,d]]])	
		feasible = True
		negative_b_idx = None
		for i in range(tableau.m - 1):
			if tableau.array[i][-1]< 0:
				feasible = False
				negative_b_idx = i 
				break
		if feasible:
			if verbose:
				print("The tableau if now feasible. End of phase 1.")
			return "success", []
		#check for bad row (i.e. b_i < 0 and all a_ij <= 0)
		for i in range(tableau.m - 1):
			if tableau.array[i][-1] < 0:
				bad_row = True
				for j in range(tableau.n -1):
					if tableau.array[i][j] > 0:
						bad_row = False
						break
				if bad_row:
					if verbose:
						print("Encountered a bad row, so phase 1 cannot be completed.")
					return "bad row", [i]
		#pivot
		pivot_column = None
		for i in range(tableau.m - 1):
			if tableau.array[negative_b_idx][i] > 0:
				pivot_column = i
				break
		pivot_row = None
		best_ratio = float("-inf")
		for i in range(negative_b_idx + 1):
			if tableau.array[i][pivot_column] < 0 or i == negative_b_idx:
				ratio = tableau.array[i][-1] / tableau.array[i][pivot_column]
				if ratio > best_ratio:
					best_ratio = ratio
					pivot_row = i	
		if verbose:
			print(f"Pivoting to activate row var {tableau.right_labels[pivot_row]} and deactivate column var {tableau.top_labels[pivot_column]}.")

		tableau.pivot(pivot_row, pivot_column)	

		if verbose:
			print(tableau)

def phase_2(tableau, verbose = False):
	'''
	starting from a tableau where the basic solution is feasible, attempts to perform pivots so that the basic solution is optimal.

	returns:
		("success", []) when the tableau is pivoted so that the basic solution is optimal.
		("bad row", [bad_row_idx]) when the initial tableau does not have feasible basic solution.
		("bad column", [bad_col_idx]) when the problem is unbounded. The variable labelling column bad_col_idx in the modified tableau may be taken to infinity.
	'''
	if verbose:
		print("----\nStarting Phase 2 of Simplex. The goal of this phase is move to an optimal vertex of the feasible region.\n----")
		print("the initial tableau is")
		print(f"{tableau}")
	#verify feasible (i.e. b>=0, where the tableau is [[A,b],[c,d]]])
	for i in range(tableau.m - 1):
		last_entry = tableau.array[i][-1]
		if last_entry < 0:
			if verbose:
				print(f"Tableau is infeasible for phase 2, the last entry of row {i}, {last_entry}, is negative.")
			return "bad row", [i]
	#phase_2 loop
	while True:
		#check if optimal
		optimal = True
		negative_c_idx = None
		for i in range(tableau.n - 1):		
			if tableau.array[-1][i] < 0:
				optimal = False
				negative_c_idx = i
				break
		if optimal:
			if verbose:
				print("The tableau is optimal. Phase 2 complete.")
			return "success", []
		#check for bad column
		for i in range(tableau.n - 1):
			if tableau.array[-1][i]< 0:
				bad_column = True
				for j in range(tableau.m - 1):
					if tableau.array[j][i] < 0:
						bad_column = False
						break
				if bad_column:
					if verbose:
						print(f"Bad column  in phase 2: {i}, program unbounded.")
					return "bad column", [i]
		#do the pivot step
		max_ratio = float("-inf")
		best_row_idx = None
		for i in range(tableau.m - 1):
			if tableau.array[i][negative_c_idx] < 0:
				ratio = tableau.array[i][-1] / tableau.array[i][negative_c_idx]
				if ratio > max_ratio:
					max_ratio = ratio
					best_row_idx = i
		if verbose:
			print(f"Pivoting to activate row var {tableau.right_labels[best_row_idx]} and deactive col var {tableau.top_labels[negative_c_idx]}.")

		tableau.pivot(best_row_idx, negative_c_idx)		

		if verbose:
			print(tableau)

def simplex(tableau, verbose = False):
	result_1 = phase_1(tableau, verbose = verbose)
	if result_1[0] != "success":
		return result_1
	result_2 = phase_2(tableau, verbose = verbose)
	if result_2[0] != "success":
		if result_2[0] == "bad row":
			raise Exception("An infeasible tableau was passed to phase 2. This should never happen, please file a bug report.")
		return result_2
	return "success", []

def get_basic_solution(tableau):
	optimal_value = tableau.array[-1][-1]				 
	primal_sol_dict = dict()
	for item in tableau.top_labels[:-1]:
		primal_sol_dict[item] = 0
	for row_idx in range(tableau.m - 1):
		primal_sol_dict[tableau.right_labels[row_idx]] = tableau.array[row_idx][-1]	
	dual_sol_dict = dict()
	for item in tableau.right_labels[:-1]:
		dual_sol_dict[item + '_dual'] = 0
	for col_idx in range(tableau.n - 1):
		dual_sol_dict[tableau.top_labels[col_idx] + '_dual'] = tableau.array[-1][col_idx]	
	return optimal_value, primal_sol_dict, dual_sol_dict

def get_unbounded_solution(tableau, bad_col_idx):
	optimal_value = float("-inf")
	primal_sol_dict = dict()
	for i, item in enumerate(tableau.top_labels[:-1]):
		if i == bad_col_idx:
			primal_sol_dict[item] = float("inf")
		else:
			primal_sol_dict[item] = 0
	for row_idx in range(tableau.m - 1):
		if tableau.array[row_idx][bad_col_idx] == 0:
			primal_sol_dict[tableau.right_labels[row_idx]] = tableau.array[row_idx][-1]
		elif tableau.array[row_idx][bad_col_idx] > 0:
			primal_sol_dict[tableau.right_labels[row_idx]] = float("inf")
		else:
			raise Exception("bad column is not actually bad! This should never happen, please file a bug report.")
	#i haven't though through what happens with the dual variables here -- and I don't think it matters.
	return optimal_value, primal_sol_dict

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
	#check all the arguments...
	assert isinstance(c, list), "parameter c must be a list"
	assert len(c) > 0, "parameter c must have length greater than 0"
	def assert_valid_A_b(A, b, A_name, b_name):
		if (A is None) and (b is None):
			return
		assert isinstance(A, list) and isinstance(b, list) and len(A) == len(b), f"{A_name} and {b_name} must be lists, and of the same length."
		for i, row in enumerate(A):
			assert isinstance(row, list), f"each element of {A_name} must be a list. the element at idx {i} is not."
			assert len(row) == len(c), f"each row of {A_name} must have length equal to the number of variables (i.e. {len(c)}). The row at idx {i} has length {len(row)}."
	assert_valid_A_b(A_g, b_g, "A_g", "b_g")
	assert_valid_A_b(A_e, b_e, "A_e", "b_e")
	assert_valid_A_b(A_l, b_l, "A_l", "b_l")
	assert isinstance(maximize, bool)
	assert isinstance(verbose, bool)

	#build tableau
	tableau_array = []
	num_vars = len(c)
	if A_g is not None:
		for i, row in enumerate(A_g):
			new_constraint = []
			for item in row:
				new_constraint.append(item)
			new_constraint.append(-b_g[i])
			tableau_array.append(new_constraint)
	if A_e is not None:
		for i, row in enumerate(A_e):
			new_constraint_1 = []
			new_constraint_2 = []
			for item in row:
				new_constraint_1.append(item)
				new_constraint_2.append(-item)
			new_constraint_1.append(-b_e[i])
			new_constraint_2.append(b_e[i])
			tableau_array.append(new_constraint_1)	
			tableau_array.append(new_constraint_2)
	if A_l is not None:
		for i, row in enumerate(A_l):
			new_constraint = []
			for item in row:
				new_constraint.append(-item)
			new_constraint.append(b_l[i])
			tableau_array.append(new_constraint)
	obj_row = []
	max_min_mult = -1 if maximize else 1
	for item in c:
		obj_row.append(max_min_mult * item)
	obj_row.append(max_min_mult * d)
	tableau_array.append(obj_row)
	
	tableau_top_labels = [f"x_{i}" for i in range(num_vars)] + ['1']
	tableau_right_labels = [f"u_{i}" for i in range(len(tableau_array) - 1)] + [f"V -->{'max' if maximize else 'min'}"] 

	tableau = Tableau(tableau_array, list(tableau_top_labels), tableau_right_labels)
	tableau.map_over_array(value_map)
	#solve with simplex
	out = simplex(tableau, verbose = verbose)
	if out[0] == "success":
		opt_val, primal_sol_dict, dual_sol_dict = get_basic_solution(tableau)
		return max_min_mult * opt_val, [primal_sol_dict[v_name] for v_name in tableau_top_labels[:-1]]
	elif out[0] == "bad row":
		return None, "infeasible program, bad row encountered in phase 1."
	elif out[0] == "bad column":
		#the program is unbounded
		optimal_value, primal_sol_dict = get_unbounded_solution(tableau, out[1][0])
		return max_min_mult * optimal_value, [primal_sol_dict[v_name] for v_name in tableau_top_labels[:-1]]
	else:
		raise Exception(f"Unknown output message from simplex method: '{out[0]}'. Please file a bug report.")
	
	
if __name__ == '__main__':
	from fractions import Fraction
	test = []
	if -1 in test: #bad arguments
		print(linprog(c = 5))
	if -2 in test:
		print(linprog(c = [7], A_l = [[2,1], [3,4,5]], b_l = [3,4]))	
	if 0 in test: #zero linear program, no constraints
		print(linprog(c = [0,0,0]))
		print("should be (0, [0, 0, 0])")
	if 1 in test: #a regular program with A_l, maximization
		print(linprog(c = [4,11], A_l = [[1,1], [2,1]], b_l = [3,4], maximize = True))
		print("should have val 33, sol [0,3]")
	if 2 in test: #a regular program with A_g
		print(linprog(c = [2,7], A_g = [[5,1], [1,3]], b_g = [5,9], verbose = True))
		print("should be 18, [9,0]")
	if 3 in test: #unbounded program
		print(linprog(c = [5,7], maximize = True))
		print("should have inf value")
		print(linprog(c = [-2,7], A_g = [[5,1], [1,3]], b_g = [5,9]))
		print(linprog(c = [-2,7], A_g = [[5,1], [1,3]], b_g = [5,9]))
		print("should have -inf value.")
	if 4 in test: #infeasible
		print(linprog(c = [1], A_g = [[1]], b_g = [3], A_l = [[1]], b_l = [0], verbose = True))
	if 5 in test: #equality constraints and d
		pass
	if 5 in test: #testing different sizes by maximizing a random linear function in a d-dimensional cube
		import random
		import time
		ds = []
		times = []
		d = 1
		e = 1
		while d < 5000:
			c = [random.random()  for i in range(d)] 
			start = time.time()
			opt_val, vars = linprog(c, A_l = [[1 if i == j else 0 for j in range(d)] for i in range(d)], b_l = d * [1], maximize = True, value_map = Fraction)	
			print(opt_val)
			end = time.time()
			ds.append(d)
			times.append(round(end - start,2))
			print(f"time for d = {d}: {round(end - start,2)} seconds")
			#print(f"{c=}")
			#print(opt_val, vars)
			d *= 2
		
		import matplotlib.pyplot as plt
		plt.plot(ds, times)
		plt.show()