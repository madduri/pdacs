
import numpy as np
import ctypes as ct
import os

#Define where the library is and load it
_path = os.path.dirname('__file__')
libpygio = ct.CDLL(os.path.abspath('dtk/lib/libpygio.so'))
#we need to define the return type ("restype") and
#the argument types
libpygio.get_elem_num.restype=ct.c_int64
libpygio.get_elem_num.argtypes=[ct.c_char_p]

libpygio.get_variable_type.restype=ct.c_int
libpygio.get_variable_type.argtypes=[ct.c_char_p,ct.c_char_p]

libpygio.read_gio_int32.restype=None
libpygio.read_gio_int32.argtypes=[ct.c_char_p,ct.c_char_p,ct.POINTER(ct.c_int)]

libpygio.read_gio_int64.restype=None
libpygio.read_gio_int64.argtypes=[ct.c_char_p,ct.c_char_p,ct.POINTER(ct.c_int64)]

libpygio.read_gio_float.restype=None
libpygio.read_gio_float.argtypes=[ct.c_char_p,ct.c_char_p,ct.POINTER(ct.c_float)]

libpygio.read_gio_double.restype=None
libpygio.read_gio_double.argtypes=[ct.c_char_p,ct.c_char_p,ct.POINTER(ct.c_double)]

libpygio.inspect_gio.restype=None
libpygio.inspect_gio.argtypes=[ct.c_char_p]

def gio_read(file_name,var_name):
    var_size = libpygio.get_elem_num(file_name)
    var_type = libpygio.get_variable_type(file_name,var_name)
    if(var_type==10):
        print("Variable not found")
        return
    elif(var_type==9):
        print("variable type not known (not int32/int64/float/double)")
    elif(var_type==0):
        #float
        result = np.ndarray(var_size,dtype=np.float32)
        libpygio.read_gio_float(file_name,var_name,result.ctypes.data_as(ct.POINTER(ct.c_float)))
        return result
    elif(var_type==1):
        #double
        result = np.ndarray(var_size,dtype=np.float64)
        libpygio.read_gio_double(file_name,var_name,result.ctypes.data_as(ct.POINTER(ct.c_double)))
        return result
    elif(var_type==2):
        #int32
        result = np.ndarray(var_size,dtype=np.int32)
        libpygio.read_gio_int32(file_name,var_name,result.ctypes.data_as(ct.POINTER(ct.c_int32)))
        return result
    elif(var_type==3):
        #int64
        result = np.ndarray(var_size,dtype=np.int64)
        libpygio.read_gio_int64(file_name,var_name,result.ctypes.data_as(ct.POINTER(ct.c_int64)))
        return result        
        
def gio_inspect(file_name):
    libpygio.inspect_gio(file_name)

