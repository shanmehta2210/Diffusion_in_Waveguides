"""Python analysis for the thermally steered OPA project.

Modules, in the order the pipeline uses them:
    constants   physical constants and fixed design choices
    comsol_io   reading COMSOL CSV exports
    thermal     spreading kernels K(x; x_h) and the heater-to-phase matrix G
    farfield    element pattern, far field across the array, beam metrics

Lengths are in micrometres throughout, to match the COMSOL exports.
"""
