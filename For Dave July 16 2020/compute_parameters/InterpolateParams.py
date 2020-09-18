import PrecomputeParams
import time
import numpy as np
from scipy.interpolate import interp1d
import copy 
import os
import sys
import shelve

def open_spline_shelf():
    return shelve.open(dir_path + os.path.sep + shelf_name) 
    
def close_spline_shelf(shelf):
    shelf.close()


# get cubic interpolating spline
# calculates cubic interpolating spline for pre-computed parameters as needed
# MUST BE ABLE TO OPEN BOTH param_shelf and spline_shelf
def get_cubic_spline(key):
    # open shelves
    param_shelf = PrecomputeParams.open_shelf()
    spline_shelf = open_spline_shelf()
    
    if key in param_shelf.keys():
        known_values = param_shelf[key]
    else:
        raise ValueError("Key not recognized in parameter shelf; could not load precomputed parameters")
    H_known, results_known = PrecomputeParams.check_validity(known_values)
    
    if len(results_known) < 2:
        raise ValueError("Too few precomputed parameter values to fit a spline")

    prev_dat_len = 0
    prev_spline = None
    # if spline already exists, load it along with the corresponding length of the fitted data
    if key in spline_shelf.keys():
        shelved = spline_shelf[key]
        prev_spline_amp = shelved[0]
        prev_spline_angle = shelved[1]
        prev_dat_len = shelved[2]
    
    # if there was no previous spline OR the previous spline was fit to fewer data points, then fit a new spline and save it in the shelf
    if len(results_known) > prev_dat_len:
        print("Calculating new spline")
        spline_log_amp    = interp1d(np.log(H_known), np.log(np.abs(results_known)), kind='linear', assume_sorted=True, fill_value=(np.log(np.abs(results_known[0])), np.log(np.abs(results_known[-1]))))
        spline_angle      = interp1d(np.log(H_known), np.unwrap(np.angle(results_known)), kind='linear', assume_sorted=True, fill_value=(np.angle(results_known[0]), np.angle(results_known[-1])))
        spline_shelf[key] = (spline_log_amp, spline_angle, len(results_known))
    else:
        print("Using old spline")
        spline_log_amp = prev_spline_amp
        spline_angle = prev_spline_angle
    
    # close shelves
    PrecomputeParams.close_shelf(param_shelf)
    close_spline_shelf(spline_shelf)
    
    return spline_log_amp, spline_angle
    

# THIS SHOULD BE THE ONLY FUNCTION NEEDED OUTSIDE OF InterpolateParams (i.e. this is the only function in the API)
# Get parameter values from interpolating spline for U or O
# - extrapolates by assuming constant beyond pre-calculated region
# - automatically checks to see if the spline needs to be updated (i.e. checks whether there are many more calculated points available)
def get_param_from_spline(wf_string, mode_indices, param_type, H):
    key = PrecomputeParams.make_clean_key(wf_string, mode_indices, param_type)
    spline_log_amp, spline_angle = get_cubic_spline(key)
    # note that interp1d in get_cubic_spline handles the logic for extrapolating
    result = np.exp(spline_log_amp(np.log(H)) + 1j * spline_angle(np.log(H)))
    return result


# Global variables
dir_path = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + "Precomputed" + os.path.sep
shelf_name = "EFTSplines"
warning_flagged = False