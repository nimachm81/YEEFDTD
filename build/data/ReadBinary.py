
import struct
from matplotlib import pyplot as plt
import time

file = open("E-x.data", mode='rb')
fileContent = file.read()

def GetNumOfArrays(fileContent):
    ind_st = 0
    size = 8
    typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])
    ind_st += size
    size = 8*3
    shape = struct.unpack("QQQ", fileContent[ind_st: ind_st + size])
    size_total = 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*8
    assert len(fileContent) % size_total == 0
    n_total = len(fileContent)/size_total
    print("number of arrays: ", n_total)
    print("Type size : ", typeSize)
    print("Shape : ", shape)
    print("size_total : ", size_total)
    print("total bytes : ", len(fileContent))
    return int(n_total)

def PrintArrays(fileContent):
    n_total = GetNumOfArrays(fileContent);
    ind_st = 0
    for n in range(n_total):
        size = 8
        typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])
        ind_st += size
        size = 8*3
        shape = struct.unpack("QQQ", fileContent[ind_st: ind_st + size])
        ind_st += size
        size_total = 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*8
        for i in range(shape[0]):
            for j in range(shape[1]):
                size = shape[2]*8
                data_ij = struct.unpack("d"*shape[2], fileContent[ind_st: ind_st + size])
                ind_st += size
                print(data_ij)

def PlotArrays(fileContent):
    plt.ion()
    n_total = GetNumOfArrays(fileContent);
    ind_st = 0
    for n in range(n_total):
        size = 8
        typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])
        ind_st += size
        size = 8*3
        shape = struct.unpack("QQQ", fileContent[ind_st: ind_st + size])
        ind_st += size
        size_total = 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*8
        for i in range(shape[0]):
            for j in range(shape[1]):
                size = shape[2]*8
                data_ij = struct.unpack("d"*shape[2], fileContent[ind_st: ind_st + size])
                ind_st += size
                if n % 100 == 0:
                    plt.plot(data_ij)
        #plt.show()
        if n % 100 == 0:
            plt.pause(0.05)
            plt.clf()
    plt.show()        


PlotArrays(fileContent)



