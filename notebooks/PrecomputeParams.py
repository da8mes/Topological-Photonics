import Wavefunctions
import shelve
import time
import EFTParams
import numpy as np
import copy 
import os
import sys

def open_shelf():
    return shelve.open(dir_path + os.path.sep + shelf_name) 
    
def close_shelf(shelf):
    shelf.close()

def make_clean_key(wf_string, mode_indices, param_type):
    wf_string, mi, pt = validate_and_clean_param_key(wf_string, mode_indices, param_type)
    key = str((wf_string, mi, pt)) 
    return key

def validate_and_clean_param_key(wf_string, mode_indices, param_type):
    mi = copy.deepcopy(mode_indices)
    
    # ensure valid parameter type, and convert "equivalent" keys into one common unique key
    if param_type in {"A", "a"}:
        proper_length = 2
        pt = "A"
    elif param_type in {"G", "g"}:
        proper_length = 2
        pt = "G"
    elif param_type in {"B", "b"}:
        proper_length = 4
        pt = "B"
    elif param_type in {"C", "c"}:
        proper_length = 4
        pt = "C"
    elif param_type in {"U", "u"}:
        proper_length = 4
        pt = "U"
    elif param_type in {"O", "o", "do", "dO", "dOmega", "domega"}:
        proper_length = 4
        pt = "O"  
    else:
        sys.exit("Parameter type not recognized")
        
    # validate mi for given wavefunction; rearrange mode indices into proper unique ordering to avoid redundant calculations
    assert len(mi) == proper_length
    if wf_string == "LLL" or wf_string == "LLL_jit":
        for i, spec in enumerate(mi):
            assert len(spec) == 1
            assert isinstance(spec[0], int)
        if pt == "U" or pt == "C":
            # first index should be >= second, third index should be >= fourth
            if mi[0][0] < mi[1][0]:
                mi = (mi[1], mi[0], mi[2], mi[3])
            if mi[2][0] < mi[3][0]:
                mi = (mi[0], mi[1], mi[3], mi[2])
        if pt == "A" or pt == "G":
            if mi[0][0] < mi[1][0]:
                mi = (mi[1], mi[0])
        if pt == "O" or pt == "B":
            if mi[0][0] < mi[1][0]:
                mi = (mi[1], mi[0], mi[2], mi[3])
            
    elif wf_string in wf_dict.keys():
        print("WARNING: Not checking for validity of inputs to this wavefunction yet")
    else:
        print("Unknown wavefunction type")
        print("Valid string inputs: " + str(wf_dict.keys()))
        sys.exit()
        
        
    return wf_string, mi, pt

# This function should only be used if you know what you're doing!
# This will delete all of the pre-computed entries for dOmega parameters.
# It is useful if you redefine O in a way that you can still rapidly recompute it from pre-computed B and A integrals
# because it will NOT delete the pre-computed values of any of those integrals
def delete_O():
    param_shelf = open_shelf()
    for key in param_shelf.keys():
        if key.endswith("O')"):
#             print(len(param_shelf[key]))
#             print(param_shelf[key][0])
            del param_shelf[key]
    
    close_shelf(param_shelf)

def compute_param(param_shelf, wf_string, mi, pt, H, tol_rel=10**-1, force_compute=False, epsabs=1.49e-08, epsrel=0.1, verbose=False):
    # ASSUMES CLEAN INPUT (ex. from call in get_param_precomputed() )
    # MAJOR MAJOR HACK
    # At the moment, force H to be a positive real number, so that when it is multiplied by 1j below it becomes pure imaginay
    # In principle I should absolutely allow arbitrary real parts of H with a positive complex component as well. But that will require
    # a future implementation because the present implementation assumes a sorted array of H values, and it is nontrivial to decide how to sort and handle the logic with complex numbers
    H = np.abs(H)
    
    global warning_flagged
    if not warning_flagged:
        print("WARNING: Major hack in compute_param to make H pure imaginary even though original passed value should be positive real. See comments in compute_param() in PrecomputeParams.py for details")
        warning_flagged = True
    
    if pt in {"A", "B", "C"}:
        if verbose:
            print("Computing integral... this can take O(minute)")
        tic = time.time_ns()
        wf = wf_dict[wf_string]
        result = function_dict[pt](wf, mi, H*1j, epsabs=epsabs, epsrel=epsrel)
        toc = time.time_ns()
        if verbose:
            print("Elapsed time computing integral = " + str((toc-tic)/10**9) + " seconds")
        used_interpolation = False
    else:
        result, used_interpolation = function_dict[pt](param_shelf, wf_string, mi, H, tol_rel=tol_rel, force_compute=force_compute, epsabs=epsabs, epsrel=epsrel)
    
    return result, used_interpolation


def compute_O(param_shelf, wf_string, mi, H, tol_rel=10**-1, force_compute=False, epsabs=1.49e-08, epsrel=0.1):
    # This computes UNITLESS O; the physical scale should be factored in elsewhere as desired
    B, ui1 = get_param_precomputed(wf_string, mi, "B", H, tol_rel=tol_rel, force_compute=force_compute, epsabs=epsabs, epsrel=epsrel, param_shelf_given = param_shelf)
    G, ui2 = get_param_precomputed(wf_string, mi[:2], "G", H, tol_rel=tol_rel, force_compute=force_compute, epsabs=epsabs, epsrel=epsrel, param_shelf_given = param_shelf)
    result = -1j*B/G/np.sqrt(1+(mi[0][0]==mi[1][0])) - ((mi[0][0]==mi[2][0])*(mi[1][0]==mi[3][0]) + (mi[0][0]==mi[3][0])*(mi[1][0]==mi[2][0])) / (1 + (mi[0][0]==mi[1][0]))
    used_interpolation = ui1 | ui2
    return result, used_interpolation

    
def compute_U(param_shelf, wf_string, mi, H, tol_rel=10**-1, force_compute=False, epsabs=1.49e-08, epsrel=0.1):
    # This computes UNITLESS U; the physical scale should be factored in elsewhere as desired
    numerator, ui0 = get_param_precomputed(wf_string, mi, "C", H, tol_rel=tol_rel, force_compute=force_compute, epsabs=epsabs, epsrel=epsrel, param_shelf_given = param_shelf)
    G1, ui1        = get_param_precomputed(wf_string, mi[:2], "G", H, tol_rel=tol_rel, force_compute=force_compute, epsabs=epsabs, epsrel=epsrel, param_shelf_given = param_shelf)
    G2, ui2        = get_param_precomputed(wf_string, mi[2:], "G", H, tol_rel=tol_rel, force_compute=force_compute, epsabs=epsabs, epsrel=epsrel, param_shelf_given = param_shelf)
    denominator    = 2 * G1 * G2 * np.sqrt((1+(mi[0][0]==mi[1][0]))*(1+(mi[2][0]==mi[3][0])))
    scale          = (1 + np.abs(H))**(2/3)
    result         = scale*numerator/denominator
    used_interpolation = ui0 | ui1 | ui2
    return result, used_interpolation


def compute_G(param_shelf, wf_string, mi, H, tol_rel=10**-1, force_compute=False, epsabs=1.49e-08, epsrel=0.1):
    A, used_interpolation = get_param_precomputed(wf_string, mi, "A", H, tol_rel=tol_rel, force_compute=force_compute, epsabs=epsabs, epsrel=epsrel, param_shelf_given = param_shelf) 
    result = np.sqrt(A / 2)
    return result, used_interpolation
    

# LINEAR interpolation used for determining whether new values need to be computed
# for nicer interpolation see InterpolateParams.py
def interpolate_precomputed(H_known, results_known, H, tol_rel, epsabs=1.49e-08, epsrel=0.1):
    # Check if there are at least two values in H_known
    #    if not, return that it is not possible
    # this check may be mildly wasteful, but I don't think this will be slow anyway in any realistic situation
    assert not H in H_known
    interpolation_valid = True
    
    if len(H_known) == 0:
        interpolation_valid = False
        result = 0
        insert_location = 0
    elif H > H_known[-1]:
        interpolation_valid = False
        result = 0
        insert_location = len(H_known)
    elif H < H_known[0]:
        interpolation_valid = False
        result = 0
        insert_location = 0
    else:
        # there should be a pair of adjacent elements, the first of which is below H and the second of which is above H
        # find them 
        # if this were slow, we could use a binary search because the list is pre-sorted
        # but I don't think it will be slow, so I'm not going to bother
        i = 1
        while H > H_known[i]:
            i += 1
        insert_location = i # for returning
        
        # now H_known[i] is the first element greater than H
        H0 = H_known[i-1]
        H1 = H_known[i]
        r0 = results_known[i-1]
        r1 = results_known[i]
        frac1 = (H-H0)/(H1-H0)
        
        result = r0*(1-frac1) + r1*frac1
        
        if r0 == 0 or r1 == 0:
            interpolation_valid = False
        elif np.abs(r1/r0) > 1+tol_rel or np.abs(r0/r1) > 1+tol_rel:
            interpolation_valid = False

    return interpolation_valid, result, insert_location
    

def check_validity(known_values):
    assert len(known_values) == 2
    assert len(known_values[0]) == len(known_values[1])
    H_known = known_values[0]
    results_known = known_values[1]
    return H_known, results_known


def get_param_precomputed(wf_string, mode_indices, param_type, H, tol_rel=10**-1, force_compute=False, epsabs=1.49e-08, epsrel=0.1, param_shelf_given=None, verbose=False):
    # If possible, returns parameter within tol_rel by interpolating between pre-computed values. Otherwise,
    # computes the integral directly and adds it to shelf of pre-computed values.
    
    # Initializing boolean flags that will be used to guide the rest of the function logic
    known_flag = False
    compute_flag = force_compute
    
    # first check for valid inputs
    wf_string, mi, pt = validate_and_clean_param_key(wf_string, mode_indices, param_type)
    key = str((wf_string, mi, pt)) # make a valid shelf key
    
    if verbose:
        print(key)
        
    # Open the shelf
    if param_shelf_given == None:
        param_shelf = open_shelf()
    else:
        param_shelf = param_shelf_given
    
    if key in param_shelf.keys():
        known_values = param_shelf[key]
    else:
        if verbose:
            print("WARNING: this combination of wavefunction, mode indices, and parameter type was not previously known and has no pre-computed values.")
        param_shelf[key] = [[],[]]
        known_values = param_shelf[key]
        
    # check validity of known_values
    # known_values should be a list containing two lists (of the same length); the first is the H values and the second is the corresponding pre-computed values
    H_known, results_known = check_validity(known_values)
        
    # check if H is in H_known
    # if not, attempt to interpolate
    if H in H_known:
        known_flag = True
        insert_location = H_known.index(H)
        result = results_known[insert_location]
        interpolation_flag = False
    else:
        interpolation_flag, result, insert_location = interpolate_precomputed(H_known, results_known, H, tol_rel)
        if not interpolation_flag: # interpolation was not possible within specified tolerance
            compute_flag = True
            
    # compute the parameter as needed
    if compute_flag:
        computed_result, used_interpolation = compute_param(param_shelf, wf_string, mi, pt, H, tol_rel=10**-1, epsabs=epsabs, epsrel=epsrel)
        
        if known_flag and not (computed_result == result):
            # if you compute AND the value already existed, make sure they match; otherwise print a SEVERE warning because something has gone terribly wrong
            print("SEVERE WARNING: Forced to re-compute a pre-existing value, but the two do not match!!")
            print("This might be a numerical comparison issue (i.e. some insignificant difference), but might indicate a SEVERE issue.")
            print("Fractional difference between old and new values = " + str((computed_result-result)/result))
            print("")
        elif not used_interpolation:
            # As long as NO STEP of the calculation used interpolation
            # INSERT COMPUTED VALUE AT insert_location (update the shelf itself)
            H_known.insert(insert_location, H)
            results_known.insert(insert_location, computed_result)
            param_shelf[key] = [H_known, results_known]
        
        result = computed_result            
    else:
        if not interpolation_flag:
            # this is merely a (hopefully redundant) check against the logic of my own code
            assert known_flag
        # if we interpolated, then DO NOT store that value in the shelf!
        
    # only close the shelf if it was opened in this function call
    if param_shelf_given == None:        
        close_shelf(param_shelf)
            
    return result, interpolation_flag
        
        
    
# Global variables
wf_dict = {"LLL": Wavefunctions.wf_lll_radial, "LLL_jit": Wavefunctions.wf_lll_radial_jit}
function_dict = {"A": EFTParams.A, "B": EFTParams.B, "C": EFTParams.C, "G": compute_G, "U": compute_U, "O": compute_O}
dir_path = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + "Precomputed" + os.path.sep
shelf_name = "EFTParams"
warning_flagged = False
        
    
        
        
    
    
    
    