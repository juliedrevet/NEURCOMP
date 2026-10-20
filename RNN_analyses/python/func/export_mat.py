# export Python results to .mat files, to plot them in MATLAB
import re
import numpy as np
from scipy.io import savemat


def _field(key):
    # valid MATLAB name (letters, digits and _, starting with a letter, at most 63 characters)
    name = re.sub(r'\W+', '_', str(key)).strip('_')
    if not name or not name[0].isalpha():
        name = 'x_' + name
    return name[:63]


def to_matlab(x, key=None):
    """
    Convert a Python object for scipy.io.savemat:
      dict -> struct (field names made valid MATLAB names), list of dicts -> struct array,
      list of strings -> cell array, None -> [], fields named 'idx' (indices) -> 1-based as in MATLAB,
      integers -> double
    """
    if isinstance(x, dict):
        return {_field(k): to_matlab(v, k) for k, v in x.items()}
    if x is None:
        return np.zeros((0, 0))
    if isinstance(x, (list, tuple)):
        if len(x) > 0 and all(isinstance(e, dict) for e in x):  # struct array
            keys = list(x[0].keys())
            s = np.zeros((1, len(x)), dtype=[(_field(k), object) for k in keys])
            for i, e in enumerate(x):
                for k in keys:
                    s[0, i][_field(k)] = to_matlab(e[k], k)
            return s
        if len(x) > 0 and all(isinstance(e, str) for e in x):  # cell array of strings
            return np.array(x, dtype=object)
        x = np.asarray(x)
    if key == 'idx':
        x = np.asarray(x) + 1
    if not isinstance(x, (bool, np.bool_)) and np.issubdtype(np.asarray(x).dtype, np.integer):
        x = np.asarray(x, dtype=float)  # integers saved as double (MATLAB default)
    return x


def export_mat(fname, **variables):
    """
    Save variables in a .mat file (one MATLAB variable per keyword argument), e.g.
    export_mat('results.mat', epochs=epochs, acc=acc, stat=stat). 1-D arrays are saved as row vectors.
    """
    savemat(fname, {_field(k): to_matlab(v, k) for k, v in variables.items()}, oned_as='row', do_compression=True)
    print(f'saved in {fname}')

