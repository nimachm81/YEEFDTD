

__all__ = ["GetNumberOfArrays", 
           "PrintArrays", "PlotArrays1D", "PlotArrays2D", "SaveFigure2D", "SaveAnimation2D",
           "GetSTArrayFrom2DArrayTimeSamples"]

import struct
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import matplotlib
import time



def GetNumberOfArrays(filename, vbose=False):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    n_total = GetNumOfArrays(fileContent, vbose);
    return n_toral

def GetNumOfArrays(fileContent, vbose=False):
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
    assert len(fileContent) % size_total == 0
    n_total = len(fileContent)/size_total
    if vbose:
        print("number of arrays: ", n_total)
        print("Type code : ", typeCode)
        print("Type size : ", typeSize)
        print("Shape : ", shape)
        print("size_total : ", size_total)
        print("total bytes : ", len(fileContent))
    return int(n_total)
    
def GetTypeCodeTypeSizeShapeNumOfArraysAndPreambleSize(fileContent):
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
    assert len(fileContent) % size_total == 0
    n_total = len(fileContent)//size_total
    preambleSize = 1*4 + 1*8 + 3*8
    return typeCode, typeSize, shape, n_total, preambleSize
    

def PrintArrays(fileName):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    n_total = GetNumOfArrays(fileContent);
    ind_st = 0
    for n in range(n_total):
        print(n, ": ")
        size = 4
        typeCode = struct.unpack("i", fileContent[ind_st: ind_st + size])[0]
        ind_st += size    
        size = 8
        typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
        ind_st += size
        size = 8*3
        shape = struct.unpack("QQQ", fileContent[ind_st: ind_st + size])
        ind_st += size
        size_total = 1*4 + 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*typeSize
        for i in range(shape[0]):
            for j in range(shape[1]):
                data_ij = None
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
                print(data_ij)
        print()


def PlotArrays1D(fileName, readEvery = 1, delay = 0.05):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    plt.ion()
    n_total = GetNumOfArrays(fileContent);
    ind_st = 0
    for n in range(n_total):
        size = 4
        typeCode = struct.unpack("i", fileContent[ind_st: ind_st + size])[0]
        ind_st += size    
        size = 8
        typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
        ind_st += size
        size = 8*3
        shape = struct.unpack("QQQ", fileContent[ind_st: ind_st + size])
        ind_st += size
        size_total = 1*4 + 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*typeSize
        for i in range(shape[0]):
            for j in range(shape[1]):
                data_ij = None
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
                if n % readEvery == 0:
                    plt.plot(np.real(data_ij))
        #plt.show()
        if n % readEvery == 0:
            plt.pause(delay)
            plt.clf()
    plt.show()        


def PlotArrays2D(fileName, readEvery = 1, delay = 0.05):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    plt.ion()
    n_total = GetNumOfArrays(fileContent)
    ind_st = 0
    for n in range(n_total):
        size = 4
        typeCode = struct.unpack("i", fileContent[ind_st: ind_st + size])[0]
        ind_st += size    
        size = 8
        typeSize = struct.unpack("Q", fileContent[ind_st: ind_st + size])[0]
        ind_st += size
        size = 8*3
        shape = struct.unpack("QQQ", fileContent[ind_st: ind_st + size])
        ind_st += size
        size_total = 1*4 + 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*typeSize
        for i in range(shape[0]):
            E = np.zeros((shape[1], shape[2]), dtype=complex)
            for j in range(shape[1]):
                data_ij = None
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
                E[j, :] = data_ij
            if n % readEvery == 0:
                print(np.max(np.abs(E)))
                plt.imshow(np.real(E), cmap="rainbow", origin='lower', aspect='auto')
                plt.colorbar()
        #plt.show()
        if n % readEvery == 0:
            plt.pause(delay)
            plt.clf()
    plt.show()        




def GetSTArrayFrom2DArrayTimeSamples(fileName):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    n_total = GetNumOfArrays(fileContent)
    typeCode, typeSize, shape, n_total, preambleSize = GetTypeCodeTypeSizeShapeNumOfArraysAndPreambleSize(fileContent);
    size_total = preambleSize + shape[0]*shape[1]*shape[2]*typeSize
    E = np.zeros((n_total, shape[0], shape[1], shape[2]), dtype=complex)

    ind_st = 0
    for n in range(n_total):
        ind_st += preambleSize
        for i in range(shape[0]):
            for j in range(shape[1]):
                data_ij = None
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
                E[n, i, j, :] = data_ij
    return E

def GetArrays(fileName, indStart, indEnd):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    n_total = GetNumOfArrays(fileContent)
    typeCode, typeSize, shape, n_total, preambleSize = GetTypeCodeTypeSizeShapeNumOfArraysAndPreambleSize(fileContent);
    size_total = preambleSize + shape[0]*shape[1]*shape[2]*typeSize
    E = None
    if typeCode==1:
        E = np.zeros((n_total, shape[0], shape[1], shape[2]), dtype=float)
    elif typeCode==2:
        E = np.zeros((n_total, shape[0], shape[1], shape[2]), dtype=float)
    elif typeCode==3 or typeCode==4:
        E = np.zeros((n_total, shape[0], shape[1], shape[2]), dtype=complex)

    ind_st = 0
    for n in range(n_total):
        ind_st += preambleSize
        for i in range(shape[0]):
            for j in range(shape[1]):
                data_ij = None
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
                E[n, i, j, :] = data_ij
    return E



###------   Below to be fixed


def SaveFigure2D(fileName, n=0):
    file = open(fileName, mode='rb')
    fileContent = file.read()
    plt.ion()
    n_total = GetNumOfArrays(fileContent)

    plt.clf()
    ind_st_arr = 0
    size = 4
    typeCode = struct.unpack("i", fileContent[ind_st_arr: ind_st_arr + size])[0]
    ind_st_arr += size    
    size = 8
    typeSize = struct.unpack("Q", fileContent[ind_st_arr: ind_st_arr + size])[0]
    ind_st_arr = size
    size = 8*3
    shape = struct.unpack("QQQ", fileContent[ind_st_arr: ind_st_arr + size])
    ind_st_arr += size
    size_total = 1*4 + 1*8 + 3*8 + shape[0]*shape[1]*shape[2]*typeSize
    
    ind_st = n*size_total + ind_st_arr
    
    fig = None
    for i in range(shape[0]):
        E = np.zeros((shape[1], shape[2]), dtype=complex)
        for j in range(shape[1]):
            data_ij = None
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
            E[j, :] = data_ij
            ind_st += size

        print(np.max(np.abs(E)))
        fig = plt.imshow(np.real(E).T, cmap="rainbow", origin='lower', aspect='auto')
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
        





