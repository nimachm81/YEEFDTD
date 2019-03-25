

__all__ = ["GetArrayInfo", "GetArrays", "ReadParamsFile", "GetDiscreteScalarDataArray", \
           "GetDiscreteVectorDataArray"]

import struct
import numpy as np

    
def GetArrayInfo(fileName):
    """ Reads the parameters of an array saved through MultiDimArray
    """
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
    """ Reads an array saved through MultiDimArray
    """
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
    else:
        assert False

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
    
    
    
def ReadParamsFile(filename_params):
    """ Reads a parameter file saved using UtilityFunctions::WriteParamToFile
    """
    paramfile = open(filename_params, mode='rb')
    paramsfileContent = paramfile.read()
    paramfile.close()
    ind_st = 0
    params = {}
    while ind_st < len(paramsfileContent):
        dsize = 0
        dtype = (struct.unpack("c", paramsfileContent[ind_st:ind_st+1])[0]).decode("utf-8") 
        dcode = None
        if dtype=='f':
            dsize = 4
            dcode = 'f'
        if dtype=='u':
            dsize = 8
            dcode = 'Q'

        #print('dtype: ', dtype)
        assert dsize > 0 and dcode != None
        ind_data = ind_st + 1
        data_i = struct.unpack(dcode, paramsfileContent[ind_data: ind_data + dsize])[0]
        ind_st += dsize + 1
        
        nameLenSize = 8
        nameLen = struct.unpack("Q", paramsfileContent[ind_st:ind_st+nameLenSize])[0]
        #print('nameLen: ', nameLen)
        ind_st += nameLenSize
        
        name = (struct.unpack("{}s".format(nameLen), paramsfileContent[ind_st:ind_st+nameLen])[0]).decode('utf-8')
        #print('name: ', name)
        ind_st += nameLen
        
        params[name] = data_i
        
    return params


def GetDiscreteScalarDataArray(fileName):
    """ Reads a discrete array saved through DiscreteScalarDataView
    """
    file = open(fileName, mode='rb')
    fileContent = file.read()
   
    preambleSize = 1*4 + 1*8 + 1*8
 
    ind_st = 0
    size = 4
    typeCode = struct.unpack("i", fileContent[ind_st: ind_st + size])[0]
    ind_st += size    
    size = 8
    typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
    ind_st += size
    size = 8
    arraySize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
    size_total = preambleSize + arraySize*typeSize

    ind = 0;
    preambleSize = 1*4 + 1*8
    arrays = []
    while ind < len(fileContent):
        ind += preambleSize
        size = 8
        arraySize = struct.unpack("Q", fileContent[ind: ind + size])[0]
        ind += size
        
        size = arraySize * typeSize
        if(typeCode == 1):
            array = struct.unpack("f"*arraySize, fileContent[ind: ind + size])
        elif(typeCode == 2):
            array = struct.unpack("d"*arraySize, fileContent[ind: ind + size])
        elif(typeCode == 3):
            array = struct.unpack("f"*2*arraySize, fileContent[ind: ind + size])
            a_r = [array[2*i] for i in range(len(array)//2) ]
            a_i = [array[2*i + 1] for i in range(len(array)//2) ]
            array = np.array(a_r) + 1j*np.array(a_i)
        elif(typeCode == 4):
            array = struct.unpack("d"*2*arraySize, fileContent[ind: ind + size])
            a_r = [array[2*i] for i in range(len(array)//2) ]
            a_i = [array[2*i + 1] for i in range(len(array)//2) ]
            array = np.array(a_r) + 1j*np.array(a_i)
        else:
            assert False
        
        ind += size
        arrays.append(np.array(array))
        
    
    file.close()

    return arrays
    
def GetDiscreteVectorDataArray(fileName):
    """ Reads a discrete array saved through DiscreteVectorDataView
    """
    file = open(fileName, mode='rb')
    fileContent = file.read()
   
    preambleSize = 1*4 + 1*8 + 1*8
 
    ind_st = 0
    size = 4
    typeCode = struct.unpack("i", fileContent[ind_st: ind_st + size])[0]
    ind_st += size    
    size = 8
    typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
    ind_st += size
    size = 8
    arraySize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
    size_total = preambleSize + arraySize*typeSize

    ind = 0;
    preambleSize = 1*4 + 1*8
    arrays = []
    while ind < len(fileContent):
        ind += preambleSize
        size = 8
        arraySize = 3 * struct.unpack("Q", fileContent[ind: ind + size])[0]
        ind += size
        
        size = arraySize * typeSize
        if(typeCode == 1):
            array = np.array(struct.unpack("f"*arraySize, fileContent[ind: ind + size]))
        elif(typeCode == 2):
            array = np.array(struct.unpack("d"*arraySize, fileContent[ind: ind + size]))
        elif(typeCode == 3):
            array = struct.unpack("f"*2*arraySize, fileContent[ind: ind + size])
            a_r = [array[2*i] for i in range(len(array)//2) ]
            a_i = [array[2*i + 1] for i in range(len(array)//2) ]
            array = np.array(a_r) + 1j*np.array(a_i)
        elif(typeCode == 4):
            array = struct.unpack("d"*2*arraySize, fileContent[ind: ind + size])
            a_r = [array[2*i] for i in range(len(array)//2) ]
            a_i = [array[2*i + 1] for i in range(len(array)//2) ]
            array = np.array(a_r) + 1j*np.array(a_i)
        else:
            assert False
        
        ind += size
        arrays.append(np.array(array.reshape((3, arraySize//3), order='F')))
        
    
    file.close()

    return arrays
    
