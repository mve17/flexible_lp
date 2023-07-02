'''
Matthew Ellison (m.gr at dartmouth.edu)
Summer 2023

Tableaus for linear programming (following Vaserstein's Introduction to Linear Programming)

You are free, with attribution, to use this code in any way.
'''
def is_zero(x, ftol):
	if isinstance(x, float):
		return -ftol < x < ftol
	else:
		return x == 0

class Tableau():
	'''
	A Tableau is an array of 'numbers', together with labels on the top and right margins. 
	the numbers must support the following operations
		- 
	'''
	def __init__(self, array, top_labels, right_labels, ftol = 1e-6):
		self.array = array 
		self.top_labels = top_labels
		self.right_labels = right_labels
		assert isinstance(self.right_labels, list)
		assert isinstance(self.top_labels, list)
		self.m = len(self.right_labels)
		self.n = len(self.top_labels)
		#check dimensions
		assert isinstance(self.array, list), "array parameter must be a list."
		assert len(self.array) == self.m, "array and right_labels must have the same length."
		for i, row in enumerate(self.array):
			assert len(row) == self.n, f"array row with idx {i} has length {len(row)}. This is not equal to the number of top labels, {self.n}."
		self.ftol = ftol
		assert isinstance(self.ftol, float) and self.ftol > 0, f"self.ftol: {self.ftol}"	
	def map_over_array(self, f):
		self.array = [[f(self.array[i][j]) for j in range(self.n)] for i in range(self.m)]
	def pivot(self, i, j, max_i = None, max_j = None):
		if max_i is not None:
			assert i <= max_i
		if max_j is not None:
			assert j <= max_j
		assert i < self.m and j <= self.n
		pivot_entry = self.array[i][j]
		assert not is_zero(pivot_entry, self.ftol), "cannot pivot on zero entry!"
		#swap margin labels
		self.top_labels[j], self.right_labels[i] = self.right_labels[i], self.top_labels[j]
		#replace pivot entry by its multiplicative inverse
		self.array[i][j] = 1 / pivot_entry
		#replace every entry 'b' in the pivot row (but not in the pivot entry) by '-b/pivot_entry'
		for k in range(self.n):
			if k != j:
				self.array[i][k] = -self.array[i][k] / pivot_entry		
		#replace every other entry 'd' by 'd - b * g / pivot_entry', where 'b' is in the pivot row and same column,
		# and 'g' in the pivot column and same row.	
		for k in range(self.m):
			if k != i:
				for l in range(self.n):
					if l != j:
						self.array[k][l] += self.array[i][l] * self.array[k][j] 
		#replace every entry 'b' in the pivot column (but not the pivot entry) by 'b/pivot_entry'
		for k in range(self.m):
			if k != i:
				self.array[k][j] = self.array[k][j] / pivot_entry
	def __str__(self):
		#decide on field width
		fw = 0
		for item in self.top_labels:
			fw = max(fw, len(item))
		for i in range(self.m):
			for j in range(self.n):
				fw = max(fw, len(str(self.array[i][j])))
		fw += 1
		def pad(string, fw = fw, c = ' ', j = 'r'):
			if j == 'r':
				return (fw - len(string)) * c + string
			elif j == 'l':
				return string + (fw - len(string)) * c
			else:
				raise Exception("Oh no!")

		out = ''
		#header row and line below
		header_row = ''
		for item in self.top_labels:
			header_row += pad(item) 
		header_row += '\n'
		out += header_row
		out += len(header_row) * '-' + '\n'		
		#data rows
		for i, row in enumerate(self.array):
			new_row = ''
			for item in row:
				new_row += pad(str(item))
			new_row += '| '
			new_row += pad(self.right_labels[i], j = 'l')
			new_row += '\n'
			out += new_row
		return out[:-1]	

if __name__ == '__main__':
	from fractions import Fraction
	array = [[1,-1,0,1,0],[-1,1,-2,-1,-2],[0,1,2,-2,-2]]
	array = [list(map(Fraction, row)) for row in array]
	top_labels = ['x_1','x_2','x_3','x_4','1']
	right_labels = ['e','f','->min']							
	tableau = Tableau(array, top_labels, right_labels)
	print(tableau)
	print('\n')
	tableau.pivot(1,2)
	print(tableau)