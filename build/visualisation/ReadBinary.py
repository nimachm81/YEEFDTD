

__all__ = ["GetArrayInfo", "GetArrays"]

import struct
import numpy as np

    
def GetArrayInfo(fileName):
    file = open(fileName, mode='rb')
    
    preambleSize = 1*4 + 1*8 + 3*8
    
    fileContent = file.read(preambleSize)

    ind_st = 0
    size = 4
    typeCode = struct.unpack("i", fileContent[ind_st: ind_st + size])[0]
    ind_st += size    
    size = 8
    typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
    ind_st += size
    size = 8*3
    shape = struct.unpack("QQQ", fileContent[ind_st: ind_st + size])
    size_total = 1*4 + 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*typeSize
    
    file.seek(0,2)
    fileSize = file.tell()
    
    assert fileSize % size_total == 0
    numOfArrays = fileSize//size_total
    
    file.close()

    return {"typeCode":typeCode , "typeSize":typeSize, "shape":shape, "numOfArrays":numOfArrays, \
            "preambleSize":preambleSize}
    

def GetArrays(fileName, indStart=0, indEnd=None):
    
    arrayInfo = GetArrayInfo(fileName)
    numOfArrays = arrayInfo["numOfArrays"]
    preambleSize = arrayInfo["preambleSize"]
    shape = arrayInfo["shape"]
    typeSize = arrayInfo["typeSize"]
    typeCode = arrayInfo["typeCode"]
    
    arraySizeInBytes = shape[0]*shape[1]*shape[2]*typeSize
    arraySizePlusPreambleInBytes = preambleSize + arraySizeInBytes
    
    file = open(fileName, mode='rb')
    
    if indStart < 0:
        indStart += numOfArrays

    file.seek(indStart*arraySizePlusPreambleInBytes)
    
    if indEnd == None:
        indEnd = numOfArrays
    if indEnd < 0:
        indEnd += numOfArrays
        
    numOfArraysToRead = indEnd - indStart
    
    fileContent = file.read(numOfArraysToRead*arraySizePlusPreambleInBytes)
    
    A = None
    if typeCode==1:
        A = np.zeros((numOfArraysToRead, shape[0], shape[1], shape[2]), dtype=np.float32)
    elif typeCode==2:
        A = np.zeros((numOfArraysToRead, shape[0], shape[1], shape[2]), dtype=np.float64)
    elif typeCode==3:
        A = np.zeros((numOfArraysToRead, shape[0], shape[1], shape[2]), dtype=np.complex64)
    elif typeCode==4:
        A = np.zeros((numOfArraysToRead, shape[0], shape[1], shape[2]), dtype=np.complex128)

    ind_st = 0
    for n in range(numOfArraysToRead):
        ind_st += preambleSize
        for i in range(shape[0]):
            for j in range(shape[1]):
                data_ij = None      ## last axis data
                size = shape[2]*typeSize
                if(typeCode == 1):
                    data_ij = struct.unpack("f"*shape[2], fileContent[ind_st: ind_st + size])
                elif(typeCode == 2):
                    data_ij = struct.unpack("d"*shape[2], fileContent[ind_st: ind_st + size])
                elif(typeCode == 3):
                    data_ij = struct.unpack("f"*(2*shape[2]), fileContent[ind_st: ind_st + size])
                    data_ij = [data_ij[2*i]+1j*data_ij[2*i+1] for i in range(shape[2])]
                elif(typeCode == 4):
                    data_ij = struct.unpack("d"*(2*shape[2]), fileContent[ind_st: ind_st + size])
                    data_ij = [data_ij[2*i]+1j*data_ij[2*i+1] for i in range(shape[2])]
                ind_st += size
                A[n, i, j, :] = data_ij
    file.close()
    return A

