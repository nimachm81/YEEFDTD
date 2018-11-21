

__all__ = ["PrintArrays", "PlotArrays1D", "PlotArrays2D", "SaveFigure2D", "SaveAnimation2D"]

import struct
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import matplotlib
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


def PlotArrays2D(fileName, readEvery = 1, delay = 0.05):
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
                plt.imshow(E.T, cmap="rainbow", origin='lower', aspect='auto')
                plt.colorbar()
        #plt.show()
        if n % readEvery == 0:
            plt.pause(delay)
            plt.clf()
    plt.show()        



def SaveFigure2D(fileName, n=0):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    plt.ion()
    n_total = GetNumOfArrays(fileContent);

    plt.clf()
    ind_st_arr = 0
    size = 8
    typeSize = struct.unpack("Q", fileContent[ind_st_arr: ind_st_arr + size])
    ind_st_arr = size
    size = 8*3
    shape = struct.unpack("QQQ", fileContent[ind_st_arr: ind_st_arr + size])
    ind_st_arr += size
    size_total = 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*8
    
    ind_st = n*size_total + ind_st_arr
    
    fig = None
    for i in range(shape[0]):
        E = np.zeros((shape[1], shape[2]))
        for j in range(shape[1]):
            size = shape[2]*8
            data_ij = struct.unpack("d"*shape[2], fileContent[ind_st: ind_st + size])
            E[j, :] = data_ij
            ind_st += size

        print(np.max(np.abs(E)))
        fig = plt.imshow(E.T, cmap="rainbow", origin='lower', aspect='auto')
        plt.colorbar()
    plt.savefig(fileName+".png", bbox_inches='tight', pad_inches=0.5)



def SaveAnimation2D(fileName, readEvery=1):
    
    file = open(fileName, mode='rb')
    fileContent = file.read()
    plt.ion()
    n_total = GetNumOfArrays(fileContent);

    fig = plt.gcf()
    def animData(n):
        plt.clf()
        ind_st_arr = 0
        size = 8
        typeSize = struct.unpack("Q", fileContent[ind_st_arr: ind_st_arr + size])
        ind_st_arr = size
        size = 8*3
        shape = struct.unpack("QQQ", fileContent[ind_st_arr: ind_st_arr + size])
        ind_st_arr += size
        size_total = 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*8
        
        ind_st = n*readEvery*size_total + ind_st_arr
        
        for i in range(shape[0]):
            E = np.zeros((shape[1], shape[2]))
            for j in range(shape[1]):
                size = shape[2]*8
                data_ij = struct.unpack("d"*shape[2], fileContent[ind_st: ind_st + size])
                E[j, :] = data_ij
                ind_st += size

            print(np.max(np.abs(E)))
            fig = plt.imshow(E.T, cmap="rainbow", origin='lower', aspect='auto')
            plt.colorbar()
        return fig
    
    anim = animation.FuncAnimation(fig, animData, frames=n_total//readEvery, interval=1, repeat=False)
    anim.save(fileName+".mp4", writer="ffmpeg", fps=15, dpi=200)
        





