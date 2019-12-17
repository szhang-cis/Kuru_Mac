import numpy as np
from warnings import warn
from .Numeric import tovoigt #, tovoigt3

__all__ = ['unique2d','in2d','itemfreq','Voigt']

#-------------------------------------------------------------------------#
# UTILITY FUNCTIONS

def unique2d(arr,axis=1,consider_sort=False,order=True,return_index=False,return_inverse=False):
    """Get unique values along an axis for a 2D array.
        see: http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        input:
            arr:
                2D array
            axis:
                Axis along which to take unique values, for instance unique
                rows (axis=1) or unique columns (axis=0). The axis ordering
                should not be confused with the usual numpy style axis argument,
                as here since unique values of a 1D type array is finally computed,
                numpy style axis hence becomes meaningless. Hence, in this context
                axis=0 implies finding unique rows of an array (lumping the rows) and
                axis=1 implies finding unique columns of an array (lumping the columns)
            consider_sort:
                Does permutation of the values in row/column matter. Two rows/columns
                can have the same elements but with different arrangements. If consider_sort
                is True then those rows/columns would be considered equal
            order:
                Similar to 1D unique in numpy, wherein the unique values are always sorted,
                if order is True, unique2d will also sort the values
            return_index:
                Similar to numpy unique. If order is True the indices would be sorted
            return_inverse:
                Similar to numpy unique
        returns:
            2D array of unique values
            If return_index is True also returns indices
            If return_inverse is True also returns the inverse array
            """

    if arr.ndim == 1:
        warn("1D is array is passed to 2D routine. Use the numpy.unique instead")
        arr = arr[:,None]

    if axis == 0:
        arr = np.copy(arr.T,order='C')

    if consider_sort is True:
        a = np.sort(arr,axis=1)
    else:
        a = arr
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))

    if return_inverse is False:
        _, idx = np.unique(b, return_index=True)
    else:
        _, idx, inv = np.unique(b, return_index=True, return_inverse=True)

    if order is True:
        idx.sort()

    if return_index == False and return_inverse == False:
        return arr[idx]
    elif return_index == True and return_inverse == False:
        return arr[idx], idx
    elif return_index == False and return_inverse == True:
        return arr[idx], inv
    elif return_index == True and return_inverse == True:
        return arr[idx], idx, inv
    else:
        return arr[idx]

def in2d(arr1, arr2, axis=1, consider_sort=False):
    """Generalisation of numpy.in1d to 2D arrays
        NOTE: arr_1 and arr_2 should have the same dtype
        input:
            arr1:
                2D array
            arr2:
                2D array
            axis:
                Axis along which to np.in1d values, for instance along
                rows (axis=1) or along columns (axis=0). The axis ordering
                should not be confused with the usual numpy style axis argument,
                as here since np.in1d values of a 1D type array is finally computed,
                numpy style axis hence becomes meaningless. Hence, in this context
                axis=0 implies finding intersection rows of an array (lumping the rows) and
                axis=1 implies finding intersection columns of an array (lumping the columns)
            consider_sort:
                Does permutation of the values in row/column matter. Two rows/columns
                can have the same elements but with different arrangements. If consider_sort
                is True then those rows/columns would be considered equal
        returns:
            1D boolean array of the same length as arr1
            """

    assert arr1.dtype == arr2.dtype

    if axis == 0:
        arr1 = np.copy(arr1.T,order='C')
        arr2 = np.copy(arr2.T,order='C')

    if consider_sort is True:
        sorter_arr1 = np.argsort(arr1)
        arr1 = arr1[np.arange(arr1.shape[0])[:,None],sorter_arr1]
        sorter_arr2 = np.argsort(arr2)
        arr2 = arr2[np.arange(arr2.shape[0])[:,None],sorter_arr2]

    arr1_view = np.ascontiguousarray(arr1).view(np.dtype((np.void, arr1.dtype.itemsize * arr1.shape[1])))
    arr2_view = np.ascontiguousarray(arr2).view(np.dtype((np.void, arr2.dtype.itemsize * arr2.shape[1])))
    intersected = np.in1d(arr1_view, arr2_view)
    return intersected.view(np.bool).reshape(-1)

def itemfreq(arr=None,un_arr=None,inv_arr=None,decimals=None):

    if (arr is None) and (un_arr is None):
        raise ValueError('No input array to work with. Either the array or its unique with inverse should be provided')
    if un_arr is None:
        if decimals is not None:
            un_arr, inv_arr = np.unique(np.round(arr,decimals=decimals),return_inverse=True)
        else:
            un_arr, inv_arr = np.unique(arr,return_inverse=True)

    if arr is not None:
        if len(arr.shape) > 1:
            dtype = type(arr[0,0])
        else:
            dtype = type(arr[0])
    else:
        if len(un_arr.shape) > 1:
            dtype = type(un_arr[0,0])
        else:
            dtype = type(un_arr[0])

    unf_arr = np.zeros((un_arr.shape[0],2),dtype=dtype)
    unf_arr[:,0] = un_arr
    unf_arr[:,1] = np.bincount(inv_arr)

    return unf_arr

def Voigt(A,sym=1):
    """Given a 4th order tensor A, puts it in 6x6 format
        Given a 3rd order tensor A, puts it in 3x6 format (symmetric wrt the first two indices)
        Given a 2nd order tensor A, puts it in 1x6 format
        sym returns the symmetrised tensor (only for 3rd and 4th order). Switched on by default
        """

    if sym==0:
        return __voigt_unsym__(A)


    if A.ndim==4:
         VoigtA = tovoigt(A)
    elif A.ndim==3:
        VoigtA = tovoigt3(A)
    elif A.ndim==2:
        VoigtA = SecondTensor2Vector(A)
        # PURE PYTHON VERSION
        # e=A
        # if e.shape[0]==3:
        #     VoigtA = 0.5*np.array([
        #         [2.*e[0,0,0],2.*e[0,0,1],2.*e[0,0,2]],
        #         [2.*e[1,1,0],2.*e[1,1,1],2.*e[1,1,2]],
        #         [2.*e[2,2,0],2.*e[2,2,1],2.*e[2,2,2]],
        #         [e[0,1,0]+e[1,0,0],e[0,1,1]+e[1,0,1],e[0,1,2]+e[1,0,2]],
        #         [e[0,2,0]+e[2,0,0],e[0,2,1]+e[2,0,1],e[0,2,2]+e[2,0,2]],
        #         [e[1,2,0]+e[2,1,0],e[1,2,1]+e[2,1,1],e[1,2,2]+e[2,1,2]],
        #         ])
        # elif e.shape[0]==2:
        #     VoigtA = 0.5*np.array([
        #         [2.*e[0,0,0],2.*e[0,0,1]],
        #         [2.*e[1,1,0],2.*e[1,1,1]],
        #         [e[0,1,0]+e[1,0,0],e[0,1,1]+e[1,0,1]]
        #         ])

    return VoigtA
