
#ifndef FDTD_FDINSTRUCTIONCODE_H_
#define FDTD_FDINSTRUCTIONCODE_H_


enum class FDInstructionCode {
    A_plusequal_sum_b_C  // Update array A as:     A += Sum{b_n*C_n}
                // A[i][j][k] += b0*C[i+i0][j+j0][k+k0] + b1*C[i+i1][j+j1][k+k1] + ... for ijk in the range
                // ind_start...ind_end
                // parameters : tuple{pair(ind_start, ind_end), A, ind_xyz_A, vector<b>, vector<C>, vector<ind_xyz_C>,
                //              vector<[i0,j0,k0]>} with types
                // std::tuple<std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,      // summation range
                //            std::string,                    // A array names
                //            int,                            // 0:A_x, 1:A_y, 2:A_z
                //            std::vector<FPNumber>,        // b scalars
                //            std::vector<std::string>,       // C arrays names
                //            std::vector<int>,               // C arrays components : 0:C_x, 1:C_y, 2:C_z
                //            std::vector<int, 3>>>    // [i0,j0,k0], ...
                //

    ,A_equal_sum_b_C    // Sets array A as:     A = Sum{b_n*C_n}
                // same parameters as A_plusequal_sum_b_C

    ,A_plusequal_sum_bB_C    // Update array A as:       A += Sum{b_n*B_n*C_n}
                // same parameters as in A_plusequal_sum_b_C except B is a grid array with the same shape as A and C

    ,A_equal_sum_bB_C    // Sets array A as:     A = Sum{b_n*B_n*C_n}
                // same parameters as in A_plusequal_sum_b_C except B is a grid array with the same shape as A and C

    ,A_equal_func_r_t   // update array A as a function of space and time        A = f(r, t)
                // The function takes the array A and time t as parameters and calculates the parameter r (position)
                // based on the location of the origin of the array eith respect to the corners of the Yee grid, and
                // then sets the elements of A based on their positions and time.
                // parameters : tuple{A, f, t}
                // std::tuple<std::string     // name of GridArrayManipulator (it has a pointer to A) and a method to
                //           >                // calculate time
    ,A_plusequal_func_r_t   // update array A as a function of space and time    A += f(r, t)
                // same parameters as in A_equal_func_r_t
    ,A_multequal_func_r_t   // update array A as a function of space and time    A *= f(r, t)
                // same parameters as in A_equal_func_r_t

    ,A_plusequal_sum_b_C_neighbor    // Sets array A from a neighbor grid as:     A += Sum{b_n*C_n}
                // same parameters as A_plusequal_sum_b_C, except it takes an extra argument as a refrence to the
                // neighbor grid
                // A belongs to this,  C belongs to neighbor
    ,A_equal_sum_b_C_neighbor    // Sets array A from a neighbor grid as:     A = Sum{b_n*C_n}
                // same parameters as A_neighbor_plusequal_sum_b_C
    ,A_to_buffer        // writes array A to a buffer
                // parameters : tuple{pair(ind_start, ind_end), A, ind_xyz_A, buffer, ind_start_buffer
    ,A_from_buffer      // reads A from buffer
                // same parameters as A_to_buffer

    ,A_plusequal_sum_b_C_shifted // same as A_plusequal_sum_b_C except if origins of A and C are not the same, A will
                // be shifted by originC - originA

};



#endif // FDTD_FDINSTRUCTIONCODE_H_
