

__all__ = ["PrintArrays", "PlotArrays1D", "PlotArrays2D"]

import struct
import numpy as np
from matplotlib import pyplot as plt
import time



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

def PrintArrays(fileName):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    n_total = GetNumOfArrays(fileContent);
    ind_st = 0
    for n in range(n_total):
        print(n, ": ")
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
        print()


def PlotArrays1D(fileName, readEvery = 1, delay = 0.05):
    file = open(fileName, mode='rb')
    fileContent = file.read()
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
                if n % readEvery == 0:
                    plt.plot(data_ij)
        #plt.show()
        if n % readEvery == 0:
            plt.pause(delay)
            plt.clf()
    plt.show()        


def PlotArrays2D(fileName, readEvery = 1):
    file = open(fileName, mode='rb')
    fileContent = file.read()
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
            E = np.zeros((shape[1], shape[2]))
            for j in range(shape[1]):
                size = shape[2]*8
                data_ij = struct.unpack("d"*shape[2], fileContent[ind_st: ind_st + size])
                E[j, :] = data_ij
                ind_st += size
            if n % readEvery == 0:
                print(np.max(np.abs(E)))
                plt.imshow(E, cmap="rainbow")
                plt.colorbar()
        #plt.show()
        if n % readEvery == 0:
            plt.pause(0.05)
            plt.clf()
    plt.show()        







