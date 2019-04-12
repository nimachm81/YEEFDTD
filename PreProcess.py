




def FindBlocks(str, begin, end):
    blocks = []
    blocksLimits = []
    while True:
        block_ind0 = str.find(begin)
        if block_ind0 >= 0:            
            block_ind1 = str.find(end, block_ind0 + len(begin))
            str_block = str[block_ind0 + len(begin): block_ind1]
            
            blocks.append(str_block)
            blockEndPrevious = 0
            if(len(blocksLimits) > 0):
                blockEndPrevious = blocksLimits[-1][1]
            blocksLimits.append([blockEndPrevious + block_ind0, blockEndPrevious + block_ind1 + len(end)])
            
            str = str[block_ind1 + len(end): ]
        else:
            break
    return blocks, blocksLimits
 


def ApplyFileSubstitution(filename_in, filename_out):
    file_in = open(filename_in, "r")
    str = file_in.read()
    file_in.close()
    blocks, blocksLimits = FindBlocks(str, r"\begin{substitution}", r"\end{substitution}")

    ind_last = 0
    str_out = ""
    for i in range(len(blocks)):
        str_out += str[ind_last: blocksLimits[i][0]]
        
        fromBlock, _ = FindBlocks(blocks[i], r"\begin{from}", r"\end{from}")
        fromList = FindBlocks(fromBlock[0], r"\begin{list}", r"\end{list}")[0][0].split(" ")
        toBlock, toBlockLimits = FindBlocks(blocks[i], r"\begin{to}", r"\end{to}")
        toList = FindBlocks(toBlock[0], r"\begin{list}", r"\end{list}")[0]
        toList = [a.split(' ') for a in toList] 
        
        blockIndStart = toBlockLimits[-1][1]
        
        for j in range(len(toList)):
            block_replaced = blocks[i][blockIndStart:]
            for k in range(len(fromList)):
                block_replaced = block_replaced.replace(fromList[k], toList[j][k])
                
            str_out += block_replaced
        ind_last = blocksLimits[i][1]
        
    str_out += str[ind_last:]
    print(str_out)
    
    file_out = open(filename_out, "w")
    file_out.write(str_out)
    file_out.close()
    

ApplyFileSubstitution("include/MultiDimArray.hpp.sub", "include/MultiDimArray.hpp")


        
