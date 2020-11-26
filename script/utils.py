##############   MULTIPROCESSING    ############
import multiprocessing
import time
from tqdm import tqdm, trange

def apply_ntimes(func, n, args, verbose=True, timeout=None):
    """
    Applies `n` times the function `func` on `args` (useful if, eg, `func` is partly random).
    Parameters
    ----------
    func : function
        func must be pickable, see https://docs.python.org/2/library/pickle.html#what-can-be-pickled-and-unpickled .
    n : int
    args : any
    timeout : int or float
        If given, the computation is cancelled if it hasn't returned a result before `timeout` seconds.
    Returns
    -------
    type
        Result of the computation of func(iter).
    """
    pool = multiprocessing.Pool()

    multiple_results = [pool.apply_async(func, args) for _ in range(n)]

    pool.close()
    
    return [res.get(timeout) for res in tqdm(multiple_results, desc='# castor.parallel.apply_ntimes', disable = True)]