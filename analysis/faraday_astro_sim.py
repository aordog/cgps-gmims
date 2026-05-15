#!/usr/bin/env python3

import numpy as np
import argparse
import json
from typing import List, Dict, Tuple, Optional, Union
from RMtools_1D.do_RMsynth_1D import run_rmsynth
from RMtools_1D.do_RMclean_1D import run_rmclean
from RMutils.util_RM import get_rmsf_planes


def make_los(parameters : Dict,  
             freqs      : Optional[Union[List[float], np.ndarray]] = None,
             noise      : float = 0.001,
             output_file: str = None) -> Dict[str, np.ndarray]:
    """

    """

    # Frequency array (read in or generate default)
    if freqs is None:
        freqs = np.arange(300e6,1801e6,1e6)
    else:
        freqs = np.array(freqs)
    lsq = ((3e8)/freqs)**2

    # Validate parameters
    check_parameters(parameters)

    # Make the complex LOS:
    complex_pol = sim_qu(parameters, lsq)

    # Add noise:
    q_noise = np.random.normal(loc=np.zeros_like(freqs), scale=noise*np.ones_like(freqs))
    u_noise = np.random.normal(loc=np.zeros_like(freqs), scale=noise*np.ones_like(freqs))
    complex_pol.real = complex_pol.real+q_noise
    complex_pol.imag = complex_pol.imag+u_noise

    # RM synthesis and CLEAN:
    rmsfplanes_out = get_rmsf_planes(lsq, np.arange(-200,200,0.1))
    rmsynth_out    = run_rmsynth([freqs, complex_pol.real, 
                                         complex_pol.imag,
                                         np.ones_like(complex_pol.real)*noise, 
                                         np.ones_like(complex_pol.imag)*noise],
                                         dPhi_radm2=0.1, phiMax_radm2=200)
    

    # Assemble results into dictionary:
    results = {'pol' :complex_pol,
               'freq':freqs,
               'lsq' :lsq,
               'fd'  :rmsfplanes_out,
               'fdf_dirty' : rmsynth_out}
    
        
    # If output file is specified, save results
    #if output_file:
    #    with open(output_file, 'w') as f:
    #        f.write("# Line of sight calculations\n")
    #        for model_name, data in results.items():
    #            f.write(f"\n# Model: {model_name}\n")
    #            np.savetxt(f, data)
    
    return results


def check_parameters(parameters: Dict) -> None:
    """
    Validate parameter dictionary structure.
    
    Parameters:
    -----------
    parameters : dict
        Dictionary containing model parameters
        
    Raises:
    -------
    ValueError
        If required parameters are missing or lists have inconsistent lengths
    """
    # Required parameters
    required_params = {"pi", "phi", "psi", "sig", "dphi"}
    
    # Check if all required parameters exist
    missing_params = required_params - set(parameters.keys())
    if missing_params:
        raise ValueError(f"Missing required parameters: {missing_params}")
    
    # Get length of first parameter's list
    first_param = list(required_params)[0]
    required_length = len(parameters[first_param])
    
    # Check if it's at least length 1
    if required_length < 1:
        raise ValueError("Parameter lists must contain at least one element")
    
    # Check if all required parameters have the same length
    inconsistent_params = [param for param in required_params 
                         if len(parameters[param]) != required_length]
    if inconsistent_params:
        raise ValueError(f"Parameters {inconsistent_params} have inconsistent length. "
                       f"All parameter lists must have the same length")
    
    # Check for extra parameters
    extra_params = set(parameters.keys()) - required_params
    if extra_params:
        print(f"Warning: Additional parameters found and will be ignored: {extra_params}")

    return



def sim_qu(parameters, lsq):
    """
    Calculate the complex polarisation for the set of input models
    to include along the LOS.
    
    Parameters:
    -----------
    parameters : dict
        Dictionary containing model parameters
    lsq : array
        Array containing wavelength squared values for sampling

    Returns:
    --------
    pol : complex array
        Array of complex values (Stoke Q - real, Stokes U - imag) 

    """
    print('')
    pol = np.zeros_like(lsq)
    
    for i in range(0,len(parameters["pi"])):

        p_set = {key: values[i] for key, values in parameters.items()}
        
        model_type, model_params = get_model_type(p_set)
        print(model_type)
        print(model_params)
        print('')

        screen = np.exp(2j*(p_set['phi']*lsq + np.radians(p_set['psi'])))
        depol  = np.exp(-2 *p_set['sig']**2 * lsq**2)
        slab   = np.sin(p_set['dphi']*lsq)/(p_set['dphi']*lsq)

        if model_type == 'Burn slab':
            pol = pol + p_set['pi']*screen*depol*slab
        else:
            pol = pol + p_set['pi']*screen*depol

    return pol



def get_model_type(param_set):
    """
    Checks type of model based on parameter values and returns the type.
    
    Parameters:
    -----------
    param_set : dict
        Dictionary containing model parameters
        
    Returns:
    -------
    model_type : str
        String describing the model type. Recognized models:
        - Faraday screen, no depolarization
        - Faraday screen with Burn depolarization
        - Burn slab
    model_params : str
        String summarizing the parameter values
    """

    if param_set["dphi"] == 0.0:
        if param_set["sig"] == 0.0:
            model_type   = 'Faraday screen, no depolarization'
            model_params = (f'PI={param_set["pi"]}, '
                            f'FD={param_set["phi"]}, '
                            f'psi0={param_set["psi"]}')
        else:
            model_type = 'Faraday screen with Burn depolarization'
            model_params = (f'PI={param_set["pi"]}, '
                            f'FD={param_set["phi"]}, '
                            f'psi0={param_set["psi"]}, '
                            f'sigma={param_set["sig"]}')
    else:
        model_type = 'Burn slab'
        model_params = (f'PI={param_set["pi"]}, '
                        f'FD={param_set["phi"]}, '
                        f'psi0={param_set["psi"]}, '
                        f'sigma={param_set["sig"]}, '
                        f'dFD={param_set["dphi"]}')
        
    return model_type, model_params



def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Calculate line-of-sight properties for specified models'
    )
    
    parser.add_argument('models', nargs='+',
                       help='One or more model names from models.json')
    parser.add_argument('--model-file', default='models.json',
                       help='Path to JSON file containing model configurations (default: models.json)')
    parser.add_argument('--output', type=str,
                       help='Output file path (optional)')
    
    args = parser.parse_args()
    
    try:
        # Call the function with command line arguments
        results = make_los(
            args.models,
            model_file=args.model_file,
            output_file=args.output
        )
        
        # Print results to console
        for model_name, data in results.items():
            print(f"\nResults for model '{model_name}':")
            print(data)
            
    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
