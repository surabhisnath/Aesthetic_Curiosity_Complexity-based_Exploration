author = 'Surabhi S Nath'
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
import imageio
import os

def matrix_to_pattern(p):
    """
    This function converts the pattern array to an image sized array
    
    Returns:
        array: scaled up image-sized 2D array
    """
    ratio = 40      # scales up the matrix `ratio` times to save as an image
    new_data = np.zeros(np.array(p.shape) * ratio)
    for j in range(p.shape[0]):
        for k in range(p.shape[1]):
            if p[j, k] == 0:    # white
                new_data[j * ratio: (j+1) * ratio, k * ratio: (k+1) * ratio] = 1    # 1 for white colour
            else:               # black
                new_data[j * ratio: (j+1) * ratio, k * ratio: (k+1) * ratio] = 0    # 0 for black colour
    new_data = np.uint8(np.round(new_data*255))
    return new_data

def neighbourhood_sum(i, j, grid, totalistic, outer_totalistic, nine_neighbourhood, five_neighbourhood):
    """
    This function evaluates the neighbourhood sum depending on the neighbourhood defined using tot/Otot and 9-neigh/5-neigh
   
    Returns:
        integer: neighbourhood sum of cell (i, j)
    """
    sum = 0
    if nine_neighbourhood:      # consider 8 neighbours
        if totalistic:      # for totalistic rules, the cell itself and neighbouring cells are considered together - i.e., the cell itself is part of the neighbourhood sum
            tuples = [(i, j), (i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1), (i + 1, j + 1), (i + 1, j - 1), (i - 1, j + 1), (i - 1, j - 1)]
        elif outer_totalistic:      # for outer-totalistic rules, the cell itself and neighbouring cells are considered separately - i.e., the cell itself is not part of the neighbourhood sum
            tuples = [(i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1), (i + 1, j + 1), (i + 1, j - 1), (i - 1, j + 1), (i - 1, j - 1)]
    elif five_neighbourhood:        # consider 4 neighbours
        if totalistic:
            tuples = [(i, j), (i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)]
        elif outer_totalistic:
            tuples = [(i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)]
    
    # calculate the sum of the cells in the neighbourhood
    for t in tuples:
        if t[0] >= 0 and t[0] < grid_size and t[1] >= 0 and t[1] < grid_size:
            sum += grid[t[0], t[1]] > 0
    
    return sum

def get_pixel(i, j, grid, totalistic, outer_totalistic, nine_neighbourhood, five_neighbourhood, powers):
    """
    This function calculates the new state of cell (i, j)

    Returns:
        integer: new state of the cell
    """
    neighsum = neighbourhood_sum(i, j, grid, totalistic, outer_totalistic, nine_neighbourhood, five_neighbourhood)
    if totalistic:
        if neighsum in powers:      # `powers` is derived from the rule number. For totalistic rules, it is an array of neighbourhood sums for which the cell is active in the next iteration.
            return 1
    elif outer_totalistic:
        if 2 * neighsum + grid[i, j] in powers:     # `powers` is derived from the rule number. For outer-totalistic rules, it is an array of (2 * neighbourhood sums + cell state) for which the cell is active in the next iteration.
            return 1
    return 0

def evolve(evolution, code, totalistic, outer_totalistic, nine_neighbourhood, five_neighbourhood, IC, initial, s1, s2, powers):
    """
    This function takes in all generation parameters such as number of evolution steps, rule code, totalistic/outer-totalistic, neighbourhood size, IC, initial configuration and evolves the CA and saves every 8th pattern in stimuli_path and the whole evolution gif in gif_path
    """

    grid = initial
    patterns = []
    stimuli = {}
    stimuli_path = "./stimuli/"
    gif_path = "./gifs/"

    for time in range(evolution):
        previous = grid.copy()

        # save every 5th pattern
        if (time + 1) % int(evolution/8) == 0:
            ret = matrix_to_pattern(grid)
            patterns.append(ret)
            grid_as_string = ','.join(str(item) for innerlist in grid for item in innerlist)
            if grid_as_string not in stimuli:
                stimuli[grid_as_string] = 1
                if not os.path.exists(stimuli_path):
                    os.makedirs(stimuli_path)
                plt.imsave(stimuli_path + s1 + "_" + s2 + "_code" + str(code) + "_IC" + str(IC) + "_iter" +  format(time + 1, '02d') + ".png", ret, cmap='gray', vmin=0, vmax=255)
                
        # update the grid
        gridcopy = np.copy(grid)
        for i in range(grid_size):
            for j in range(grid_size):
                gridcopy[i, j] = get_pixel(i, j, grid, totalistic, outer_totalistic, nine_neighbourhood, five_neighbourhood, powers)
        grid = gridcopy.copy()
        
        # if no change in the pattern, save it and break execution
        if (previous == grid).all():
            ret = matrix_to_pattern(grid)
            patterns.append(ret)
            grid_as_string = ','.join(str(item) for innerlist in grid for item in innerlist)
            if grid_as_string not in stimuli:
                stimuli[grid_as_string] = 1
                if not os.path.exists(stimuli_path):
                    os.makedirs(stimuli_path)
                plt.imsave(stimuli_path + s1 + "_" + s2 + "_code" + str(code) + "_IC" + str(IC) + "_iter" +  format(time + 1, '02d') + ".png", ret, cmap='gray', vmin=0, vmax=255)
            break
    
    # Save the whole evolution as a gif
    if not os.path.exists(gif_path):
        os.makedirs(gif_path)
    imageio.mimsave(gif_path + s1 + "_" + s2 + "_code" + str(code) + "_IC" + str(IC) + ".gif", patterns, duration = 0.2)

def generate(grid_size, evolution, rules):
    """
    This function generates and saves all stimuli and evolution gifs into 2 folders - stimuli and gifs in the generator folder.

    Args:
        grid_size (integer): grid size of the stimuli
        evolution (integer): number of evolution steps
        rules (array): list of rules, where each rule is of the form [rule code, tot/Otot, 9-niegh/5-neigh, IC]
    """
    for p in rules:

        # get all parameters from the rule
        code = p[0]
        totalistic = p[1]
        outer_totalistic = not(totalistic)
        nine_neighbourhood = p[2]
        five_neighbourhood = not(nine_neighbourhood)
        IC = p[3]

        # set `initial` configuration based on IC
        if IC == 1:
            initial = np.zeros((grid_size, grid_size), int)
            initial[int(grid_size/2), int(grid_size/2)] = 1
        elif IC == 2:
            initial = np.zeros((grid_size, grid_size), int)
            initial[int(grid_size/2) - 1, int(grid_size/2)] = 1
            initial[int(grid_size/2) + 1, int(grid_size/2)] = 1
            initial[int(grid_size/2), int(grid_size/2) + 1] = 1
            initial[int(grid_size/2) + 1, int(grid_size/2) + 1] = 1
            initial[int(grid_size/2) + 1, int(grid_size/2) - 1] = 1
        if IC == 3:
            initial = np.random.randint(2, size=(grid_size, grid_size))

        # strings for filename while saving
        if nine_neighbourhood:
            s1 = "9neigh"
        elif five_neighbourhood:
            s1 = "5neigh"
        if totalistic:
            s2 = "tot"
        elif outer_totalistic:
            s2 = "outertot"

        # Derive rule from code
        binary = "{0:b}".format(code)
        l = len(binary)
        powers = [l -  1 - m.start() for m in re.finditer('1', binary)]

        # Evolve CA and save patterns
        evolve(evolution, code, totalistic, outer_totalistic, nine_neighbourhood, five_neighbourhood, IC, initial, s1, s2, powers)

if __name__ == "__main__":
    grid_size = 27
    evolution = 40
    states = [0, 1] # 2 states, white and black

    # Rules of the form [rule code, tot/Otot, 9-niegh/5-neigh, IC]
    rules = [
        # IC_1
        [451, False, False, 1], [452, False, False, 1], [453, False, False, 1],
        [454, False, False, 1], [457, False, False, 1], [459, False, False, 1],
        [461, False, False, 1], [462, False, False, 1], [465, False, False, 1],
        [467, False, False, 1], [468, False, False, 1], [469, False, False, 1],
        [470, False, False, 1], [475, False, False, 1], [478, False, False, 1],
        [481, False, False, 1], [483, False, False, 1], [485, False, False, 1],
        [489, False, False, 1], [491, False, False, 1], [493, False, False, 1],
        [497, False, False, 1], [1022, False, False, 1], [510, False, False, 1],
        [374, False, False, 1], [174, False, False, 1], [494, False, False, 1],
        [750, False, False, 1], [15822, False, True, 1], [699054, False, True, 1],
        [191044, False, True, 1], [93737, False, True, 1],[85507, False, True, 1],

        # IC_2
        [451, False, False, 2], [452, False, False, 2], [453, False, False, 2],
        [454, False, False, 2], [457, False, False, 2], [459, False, False, 2],
        [461, False, False, 2], [462, False, False, 2], [465, False, False, 2],
        [467, False, False, 2], [468, False, False, 2], [469, False, False, 2],
        [470, False, False, 2], [475, False, False, 2], [478, False, False, 2],
        [481, False, False, 2], [483, False, False, 2], [485, False, False, 2],
        [489, False, False, 2], [491, False, False, 2], [493, False, False, 2],

        # uncomment for more variations
        # [497, False, False, 2], [451, False, True, 2], [452, False, True, 2],
        # [453, False, True, 2], [454, False, True, 2], [457, False, True, 2],
        # [459, False, True, 2], [461, False, True, 2], [462, False, True, 2],
        # [465, False, True, 2], [467, False, True, 2], [468, False, True, 2],
        # [469, False, True, 2], [470, False, True, 2], [475, False, True, 2],
        # [478, False, True, 2], [481, False, True, 2], [483, False, True, 2],
        # [485, False, True, 2], [489, False, True, 2], [491, False, True, 2],
        # [493, False, True, 2], [497, False, True, 2], [451, True, False, 2],
        # [452, True, False, 2], [453, True, False, 2], [454, True, False, 2],
        # [457, True, False, 2], [459, True, False, 2], [461, True, False, 2],
        # [462, True, False, 2], [465, True, False, 2], [467, True, False, 2],
        # [468, True, False, 2], [469, True, False, 2], [470, True, False, 2],
        # [475, True, False, 2], [478, True, False, 2], [481, True, False, 2],
        # [483, True, False, 2], [485, True, False, 2], [489, True, False, 2],
        # [491, True, False, 2], [493, True, False, 2], [497, True, False, 2],
        
        [750, False, False, 2], [15822, False, True, 2], [699054, False, True, 2],
        [191044, False, True, 2], [93737, False, True, 2], [85507, False, True, 2],
        [256746, False, True, 2], [736, False, True, 2], [291552, False, True, 2],
        [184, False, True, 2], [46, False, True, 2],

        # IC_3
        [24, True, False, 3], [510, True, False, 3], [52, True, False, 3],
        [736, False, True, 3], [196623, False, True, 3], [152822, False, True, 3],
        [143954, False, True, 3], [3276, False, True, 3], [224, False, True, 3],
        [4, True, False, 3], [6, True, False, 3], [30, True, False, 3],
        [38, True, False, 3], [52, True, False, 3], [54, True, False, 3],
        [56, True, False, 3], [58, True, False, 3]
    ]

    generate(grid_size, evolution, rules)