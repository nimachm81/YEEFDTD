

__all__ = ["Write1DArray"]


import numpy as np
import struct


def Write1DArray(arr, filename, format="d"):
    ## f: float     d: double
    assert format in ["f", "d"]
    n = len(arr)
    file = open(filename, 'wb')
    
    for i in range(len(arr)):
        file.write(struct.pack(format, arr[i]))
        
    file.close()
    
    


