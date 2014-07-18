import numpy as np

class interp1d(object):
    """
    Mimics scipy.interpolate interp1d but is pickleable for multiprocessing module
    """
    def __init__(self,x,y,copy=True,bounds_error=True,fill_value=np.nan,assume_sorted=False):
        self.copy = copy
        self.bounds_error = bounds_error
        self.fill_value = fill_value
        x = np.array(x,copy=self.copy)
        y = np.array(y,copy=self.copy)
        if not assume_sorted:
            ind = np.argsort(x)
            x = x[ind]
            y = y[ind] #np.take(y,ind,axis=-1)
        if not issubclass(y.dtype.type,np.inexact):
            y = y.astype(np.float_)
        self.x = x
        self.y = y
        self._y = y
    
    def _evaluate(self,x_new):
        x_new = np.asarray(x_new)
        out_of_bounds = self._check_bounds(x_new)
        y_new = self._call(x_new)
        if len(y_new) > 0:
            y_new[out_of_bounds] = self.fill_value
        return y_new

    def _call(self,x_new):
        x_new_indices = np.searchsorted(self.x,x_new)
        x_new_indices = x_new_indices.clip(1, len(self.x)-1).astype(int)
        lo = x_new_indices - 1
        hi = x_new_indices
        x_lo = self.x[lo]
        x_hi = self.x[hi]
        y_lo = self._y[lo]
        y_hi = self._y[hi]
        slope = (y_hi-y_lo) / (x_hi-x_lo)#[:,None]
        return slope*(x_new-x_lo) + y_lo
        #return slope*(x_new-x_lo)[:,None] + y_lo
    
    def _check_bounds(self,x_new):
        below_bounds = x_new < self.x[0]
        above_bounds = x_new > self.x[-1]
        if self.bounds_error and below_bounds.any():
            raise ValueError("A value in x_new is below the interpolation "
                "range.")
        if self.bounds_error and above_bounds.any():
            raise ValueError("A value in x_new is above the interpolation "
                "range.")
        out_of_bounds = np.logical_or(below_bounds, above_bounds)
        return out_of_bounds
    
    def __call__(self,x):
        x = np.asarray(x)
        if not np.issubdtype(x.dtype,np.inexact):
            x = x.astype(float)
        x = x.ravel() #x_shape = x.shape
        return self._evaluate(x)
