author = 'Surabhi S Nath'
import numpy as np
import pytest

"""
This module implements the local spatial complexity measure based on Javaheri Javid, M. A. (2019). Aesthetic Automata: Synthesis and Simulation of Aesthetic Behaviour in Cellular Automata (Doctoral dissertation, Goldsmiths, University of London).
Refer to Nath et al., 2023, Relating Objective Complexity, Subjective Complexity and Beauty in Binary Pixel Patterns, PACA, Section 2.2 (3) and AII for details.
"""

def get_i_j_(dirn, i, j):
	if dirn == 1:
		return i, j + 1
	if dirn == 2:
		return i, j - 1
	if dirn == 3:
		return i + 1, j
	if dirn == 4:
		return i - 1, j
	if dirn == 5:
		return i + 1, j + 1
	if dirn == 6:
		return i + 1, j - 1
	if dirn == 7:
		return i - 1, j + 1
	if dirn == 8:
		return i - 1, j - 1
	return None, None

def get_tuples(col1, col2, dirn, grid, grid_size):
	"""
	This function calculates the numerator of the conditional and joint distributions and the denominator of the conditional distribution

	Args:
		col1 (integer): state of the first pixel
		col2 (integer): state of the second pixel
		dirn (integer): direction of the second pixel with respect to the first  
		grid (array): grid for local spatial complexity computation
		grid_size (integer): grid size of the grid

	Returns:
		tuple: numerator of the conditional and joint distributions, and the denominator of the conditional distribution
	"""
	count_tuples = 0
	count_neigh = 0
	for i in range(grid_size):
		for j in range(grid_size):
			i_, j_ = get_i_j_(dirn, i, j)
			if i_ >= 0 and i_ < grid_size and j_ >= 0 and j_ < grid_size:
				if grid[i_, j_] == col2:
					count_neigh += 1
				if grid[i, j] == col1 and grid[i_, j_] == col2:
					count_tuples += 1
	return count_neigh, count_tuples

def get_joint_conditional(col1, col2, dirn, grid, grid_size):
	"""
	This function evaluates the condinational and joint probabilities given the 2 colours as stated in eq. 4

	Returns:
		tuple: conditional and joint probabilities
	"""
	den1 = grid_size * (grid_size - 1)
	den2, num = get_tuples(col1, col2, dirn, grid, grid_size)
	if den2 == 0:
		return 0, 0
	conditional = num/den2
	joint = num/den1
	return conditional, joint

def calculate_local_spatial_complexity(grid, grid_size):
	"""
	This function provides an implementation for the calculation of Spatial Complexity and Local Asymmetry as described in Javid, 2019.

	Args:
		grid (array): pattern whose calculate local spatial complexity is being calculated
		grid_size (integer): grid size of the pattern

	Returns:
		tuple: local spatial complexity and local asymmetry
	"""
	grid = grid.reshape(grid_size, grid_size)
	states = np.unique(grid)
	dirns = [1, 2, 3, 4, 5, 6, 7, 8]	# set of 8 directions
	
	# Calculate spatial complexity
	Gs = []
	muG = 0
	for dirn in dirns:
		G = 0
		for state1 in states:
			for state2 in states:
				conditional, joint = get_joint_conditional(state1, state2, dirn, grid, grid_size)
				if conditional != 0:
					G += -1 * joint * np.log2(conditional)
		muG += G
		Gs.append(G)

	local_spatial_complexity = muG/len(dirns)

	# Calculate local asymmetry
	horizontal_asymm = max(Gs[0], Gs[1]) - min(Gs[0], Gs[1])
	vertical_asymm = max(Gs[2], Gs[3]) - min(Gs[2], Gs[3])
	secdiag_asymm = max(Gs[5], Gs[6]) - min(Gs[5], Gs[6])
	primdiag_asymm = max(Gs[4], Gs[7]) - min(Gs[4], Gs[7])

	local_asymm = (horizontal_asymm + vertical_asymm + secdiag_asymm + primdiag_asymm)/4		# local asymmetry is considered to be the mean of asymmetries along the 4 directions	

	return local_spatial_complexity, local_asymm

if __name__ == "__main__":

	grid_size = 27

	# calculate density for some example patterns:
	#    1) Full black pattern
	#    2) White grid with one central black cell
	#    3) checkerboard
	#    4) Random pattern

	all_black = np.ones(grid_size * grid_size, dtype=int)
	assert calculate_local_spatial_complexity(all_black, grid_size) == (0.0, 0.0)
	print(calculate_local_spatial_complexity(all_black, grid_size))

	grid_with_one_centre = np.zeros(grid_size * grid_size, dtype=int)
	grid_with_one_centre[(grid_size * grid_size) //2] = 1
	assert calculate_local_spatial_complexity(grid_with_one_centre, grid_size) == (0.015480986966185604, 0.0)
	print(calculate_local_spatial_complexity(grid_with_one_centre, grid_size))

	checkerboard = np.zeros(grid_size * grid_size, dtype=int)
	checkerboard[1::2] = 1
	assert calculate_local_spatial_complexity(checkerboard, grid_size) == (0.0, 0.0)
	print(calculate_local_spatial_complexity(checkerboard, grid_size))

	random = np.random.choice([0, 1], size=(grid_size, grid_size))
	print(calculate_local_spatial_complexity(random, grid_size))