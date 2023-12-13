import pytest
import sys
sys.path.append("./local spatial complexity/")
sys.path.append("./intricacy/")
from calculate_local_spatial_complexity import *
from calculate_intricacy import *

grid_size = 27

def test_LSC():

    all_black = np.ones(grid_size * grid_size, dtype=int)
    assert calculate_local_spatial_complexity(all_black, grid_size) == (0.0, 0.0)

    grid_with_one_centre = np.zeros(grid_size * grid_size, dtype=int)
    grid_with_one_centre[(grid_size * grid_size) //2] = 1
    assert calculate_local_spatial_complexity(grid_with_one_centre, grid_size) == (0.015480986966185604, 0.0)

    checkerboard = np.zeros(grid_size * grid_size, dtype=int)
    checkerboard[1::2] = 1
    assert calculate_local_spatial_complexity(checkerboard, grid_size) == (0.0, 0.0)


def test_intricacy():
    
    all_black = np.ones(grid_size * grid_size)
    assert calculate_intricacy(all_black, grid_size) == (1, 1)

    grid_with_one_centre = np.zeros(grid_size * grid_size)
    grid_with_one_centre[(grid_size * grid_size) //2] = 1
    assert calculate_intricacy(grid_with_one_centre, grid_size) == (2, 2)

    checkerboard = np.zeros(grid_size * grid_size)
    checkerboard[1::2] = 1
    assert calculate_intricacy(checkerboard, grid_size) == (729, 2)