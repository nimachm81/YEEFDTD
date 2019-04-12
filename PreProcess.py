

def FindSubstitutionBlock(str):
    blocks = []
    block_ind0 = str.find(r"\begin{substitution}")
    if ind_start >= 0:
        params_ind0 = str.find(r"{", block_ind0)
        assert params_ind0 >= 0
        params_ind1 = str.find(r"}", params_ind0)
        assert params_ind1 >= 0
        
        str_params = str[params_ind0 + 1: params_ind1]
        paramsFrom_ind0 = str_params.find(r"\subFrom")
        assert paramsFrom_ind0 > 0
        paramsFrom_ind1 = str_params.find(r"\subTo")
        
        block_ind1 = str.find(r"\end{substitution}", block_ind0)
