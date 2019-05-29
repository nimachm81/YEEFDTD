
__all__ = ["GridBlock", "GridCollectionner"]

import json
import re

def MultiWordReplace(text, wordDic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """
    rc = re.compile('|'.join(map(re.escape, wordDic)))
    def translate(match):
        return wordDic[match.group(0)]
    return rc.sub(translate, text)

class GridBlock:
    def __init__(self, name, blockLevel, blockPosition):
        """ blockLevel: 0, 1, 2
            blockPosition: "c" (center), "r" (right), "l" (left), "u" (up), "d" (down), "f" (front), "b" (back)
        """
        self.name = name
        self.blockLevel = blockLevel
        self.blockPosition = blockPosition
        self.connections = {}
        self.sources = []
        self.views = []
    
    def SetCorners(self, r0, r1):
        self.r0 = r0
        self.r1 = r1
    
    def SetNumOfCells(self, nx, ny, nz):
        self.nx = int(nx)
        self.ny = int(ny)
        self.nz = int(nz)
        
    def SetCellDimentions(self, dx, dy, dz):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        
    def SetTimeStep(self, dt):
        self.dt = dt
        
    def AddConnection(self, face, gridBlock):
        """ face: "r" (right), "l" (left), "u" (up), "d" (down), "f" (front), "b" (back)
        """
        self.connections[face] = gridBlock
        
    def AddSource(self, sourceParams):
        self.sources.append(sourceParams)

    def AddView(self, viewParams):
        self.views.append(viewParams)
        
    def SetupSources(self, grid):
        ## grid : json grid 
        for sourceIndex in range(len(self.sources)):
            sourceParams = self.sources[sourceIndex]
            if sourceParams["type"] == "GaussianPointSource":
                json_file = open('sources/GaussianInTime/J_point_3D.json')
                file_content = json_file.read()
                json_file.close()
                
                ind_x, ind_y, ind_z = int(sourceParams["ind_x"]), int(sourceParams["ind_y"]), int(sourceParams["ind_z"])
                polarization = '"' + sourceParams["polarization"] + '"'
                amplitude = sourceParams["amplitude"]
                t_center = sourceParams["t_center"]
                t_decay = sourceParams["t_decay"]
                modulationFrequency = sourceParams["modulationFrequency"]
                modulationPhase = sourceParams["modulationPhase"]
                timeOffsetFraction = sourceParams["timeOffsetFraction"]
                dV = self.dx*self.dy*self.dz
                dA = None
                if sourceParams["polarization"] == 'x':
                    dA = self.dy*self.dz
                elif sourceParams["polarization"] == 'y':
                    dA = self.dx*self.dz
                elif sourceParams["polarization"] == 'z':
                    dA = self.dy*self.dx
                else:
                    assert False
                    
                
                replaceDic = {'"_indxJ_"': str(ind_x), '"_indyJ_"': str(ind_y), '"_indzJ_"': str(ind_z),
                              '"_indxJ_p1_"': str(ind_x + 1), '"_indyJ_p1_"': str(ind_y + 1), '"_indzJ_p1_"': str(ind_z + 1),
                              '"_j_polarization_"':polarization, '"_j_amplitude_"': str(amplitude),
                              '"_j_t_center_"': str(t_center), '"_j_t_decay_"': str(t_decay),
                              '"_j_mod_freq_"': str(modulationFrequency), '"_j_mod_phase_"': str(modulationPhase),
                              '"_j_time_offset_frac_"': str(timeOffsetFraction),
                              '"_m_dt_dA_"': str(-self.dt/dA),
                              '"_m_dt_dV_"': str(-self.dt/dV),
                              '"J"': '"J{}"'.format(sourceIndex),
                              '"Jupdater"': '"Jupdater{}"'.format(sourceIndex),
                              '"J_update"': '"J_update{}"'.format(sourceIndex),
                              '"E_me_J"': '"E_me_J{}"'.format(sourceIndex)
                              }
                file_content = MultiWordReplace(file_content, replaceDic)
                currentData = json.loads(file_content)
                
                grid["partialGridArrays"].extend(currentData["partialGridArrays"])
                grid["girdArrayManipulators"].extend(currentData["girdArrayManipulators"])
                grid["updateInstructions"].extend(currentData["updateInstructions"])
                for j_updateSequence in currentData["updateSequences"]:
                    seqName = j_updateSequence["name"]
                    for updateSequence in grid["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(j_updateSequence["sequence"])
        
        
    def SetupViews(self, grid):
        ## grid : json grid 
        for viewIndex in range(len(self.views)):
            viewParams = self.views[viewIndex]
            if viewParams["type"] == "partial":
                json_file = open('views/partial.json')
                file_content = json_file.read()
                json_file.close()
                
                plane = viewParams["plane"]             ## "x", "y", "z"
                at = viewParams["at"]                   ## 1.2  ---> x = 1.2
                direction = viewParams["direction"]     ## "x", "y", "z"
                
                indx0, indy0, indz0 = [None]*3
                indx1, indy1, indz1 = [None]*3
                
                if plane == "x":
                    indx0 = int((at - self.r0[0])/self.dx)
                    indx1 = indx0 + 1
                    indy0 = 0
                    indz0 = 0
                    indy1 = self.ny
                    indz1 = self.nz
                elif plane == "y":
                    indy0 = int((at - self.r0[1])/self.dy)
                    indy1 = indy0 + 1
                    indx0 = 0
                    indz0 = 0
                    indx1 = self.nx
                    indz1 = self.nz
                elif plane == "z":
                    indz0 = int((at - self.r0[2])/self.dz)
                    indz1 = indz0 + 1
                    indy0 = 0
                    indx0 = 0
                    indy1 = self.ny
                    indx1 = self.nx
                    
                if indx0 >= 0 and indx0 <= self.nx and indy0 >= 0 and indy0 <= self.ny and indz0 >= 0 and indz0 <= self.nz:
                    saveRate = viewParams['saveRate']
                    arrayName = viewParams['arrayName']
                    fileName = viewParams['fileName']
                    
                    replaceDic = {'"_indxSave_0_"': str(indx0), '"_indySave_0_"': str(indy0), '"_indzSave_0_"': str(indz0),
                                  '"_indxSave_1_"': str(indx1), '"_indySave_1_"': str(indy1), '"_indzSave_1_"': str(indz1),
                                  '"_save_rate_"': str(saveRate),
                                  '"_array_"': '"' + arrayName + '"', 
                                  '"_direction_"': '"' + direction + '"', 
                                  '"_file_name_"': '"' + fileName + '"'
                                  }
                    file_content = MultiWordReplace(file_content, replaceDic)
                    viewData = json.loads(file_content)
                    
                    grid["gridViews"].extend(viewData["gridViews"])
            
        
    def SetupGrid(self):
        ##---------------------------------------------- grid_m ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        if self.blockPosition == "c":
            json_file = open('layer0/grid_m/grid_m.json')
            file_content = json_file.read()
            json_file.close()
            
            nx, ny, nz = self.nx, self.ny, self.nz
            r0, r1 = self.r0, self.r1
            dt, dx, dy, dz = self.dt, self.dx, self.dy, self.dz
            gm_name = '"' + self.name + '"'
            replaceDic = {'"_nx_"': str(nx), '"_ny_"': str(ny), '"_nz_"': str(nz),
                          '"_nx_p1_"': str(nx + 1), '"_ny_p1_"': str(ny + 1), '"_nz_p1_"': str(nz + 1),
                          '"_x0_"': str(r0[0]), '"_y0_"': str(r0[1]), '"_z0_"': str(r0[2]),
                          '"_x1_"': str(r1[0]), '"_y1_"': str(r1[1]), '"_z1_"': str(r1[2]),
                          '"_dt_dx_"': str(dt/dx), '"_dt_dy_"': str(dt/dy), '"_dt_dz_"': str(dt/dz),
                          '"_m_dt_dx_"': str(-dt/dx), '"_m_dt_dy_"': str(-dt/dy), '"_m_dt_dz_"': str(-dt/dz),
                          '"_dt_"': str(dt),
                          '"grid_m"': gm_name
                          }
            file_content = MultiWordReplace(file_content, replaceDic)
            grids = json.loads(file_content)
            grid_m = grids[self.name]
            
            self.SetupSources(grid_m)
            self.SetupViews(grid_m)
            
            connections = self.connections
            
            #------------------- r ------------------
            #----------------------------------------
            if "r" in connections:
                json_file = open('layer0/grid_m/connections/grid_r.json')
                file_content = json_file.read()
                json_file.close()
                
                grBlock = connections["r"]
                gr_nx, gr_ny = grBlock.nx, grBlock.ny
                gr_name = '"' + grBlock.name + '"'
                gr_replaceDic = {'"_gr_nx_m1_"': str(gr_nx - 1), '"_gr_ny_m1_"': str(gr_ny - 1), 
                                 '"_nx_m1_"': str(nx - 1), '"_ny_m1_"': str(ny - 1), 
                                 '"_nx_m2_"': str(nx - 2), '"_ny_m2_"': str(ny - 2), 
                                 '"grid_r"': gr_name
                                 }
                
                file_content = MultiWordReplace(file_content, replaceDic)
                file_content = MultiWordReplace(file_content, gr_replaceDic)
                grid_r_conn = json.loads(file_content)
                
                if "u" in connections or "d" in connections:
                    guBlock = connections["u"]
                    gdBlock = connections["d"]
                    gu_nx, gu_ny, gu_nz = guBlock.nx, guBlock.ny, guBlock.nz
                    gd_ny = gdBlock.ny
                    gu_name = guBlock.name
                    gd_name = gdBlock.name
                    gl_nz = connections["l"].nz
                    gr_ud_replaceDic = {
                        '"_gd_ny_m1_"': str(gd_ny - 1),
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                        '"grid_u"': '"' + gu_name + '"',
                        '"grid_d"': '"' + gd_name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gr_ud_replaceDic)
                    grid_r_conn = json.loads(file_content)
                
                if "f" in connections or "b" in connections:
                    gfBlock = connections["f"]
                    gbBlock = connections["b"]
                    gf_nx, gf_ny, gf_nz = gfBlock.nx, gfBlock.ny, gfBlock.nz
                    gb_nx, gb_ny, gb_nz = gbBlock.nx, gbBlock.ny, gbBlock.nz
                    gd_ny = connections["d"].ny
                    gr_fb_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gb_nx_m1_"': str(gb_nx - 1),
                        '"grid_f"': '"' + gfBlock.name + '"',
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gr_fb_replaceDic)
                    grid_r_conn = json.loads(file_content)
                    
                
                grid_m["updateInstructions"].extend(grid_r_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_m["updateInstructions"].extend(grid_r_conn["updateInstructions"]["right_only"])
                elif "f" not in connections and "b" not in connections:
                    grid_m["updateInstructions"].extend(grid_r_conn["updateInstructions"]["right_up_down"])
                else:
                    grid_m["updateInstructions"].extend(grid_r_conn["updateInstructions"]["right_up_down_front_back"])
                
                for gr_updateSequence in grid_r_conn["updateSequences"]:
                    seqName = gr_updateSequence["name"]
                    for updateSequence in grid_m["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gr_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(gr_updateSequence["sequence"]["right_only"])
                            elif "f" not in connections and "b" not in connections:
                                updateSequence["sequence"].extend(gr_updateSequence["sequence"]["right_up_down"])
                            else:
                                updateSequence["sequence"].extend(gr_updateSequence["sequence"]["right_up_down_front_back"])

            #------------------- l ------------------
            #----------------------------------------
            if "l" in connections:
                json_file = open('layer0/grid_m/connections/grid_l.json')
                file_content = json_file.read()
                json_file.close()
                
                glBlock = connections["l"]
                gl_nx, gl_ny, gl_nz = glBlock.nx, glBlock.ny, glBlock.nz
                gl_name = '"' + glBlock.name + '"'
                gl_replaceDic = {'"_gl_nx_m1_"': str(gl_nx - 1), '"_gl_ny_m1_"': str(gl_ny - 1), 
                                 '"_gl_nz_"': str(gl_nz),
                                 '"_nx_m1_"': str(nx - 1), '"_ny_m1_"': str(ny - 1), 
                                 '"_nx_m2_"': str(nx - 2), '"_ny_m2_"': str(ny - 2), 
                                 '"grid_l"': gl_name
                                 }
                
                file_content = MultiWordReplace(file_content, replaceDic)
                file_content = MultiWordReplace(file_content, gl_replaceDic)
                grid_l_conn = json.loads(file_content)

                if "u" in connections or "d" in connections:
                    guBlock = connections["u"]
                    gdBlock = connections["d"]
                    gu_nx, gu_ny, gu_nz = guBlock.nx, guBlock.ny, guBlock.nz
                    gu_name = guBlock.name
                    gd_name = gdBlock.name
                    gl_ud_replaceDic = {
                        '"_gd_ny_m1_"': str(gd_ny - 1),
                        '"_gd_ny_p_ny2_"': str(gd_ny + int(ny/2)),
                        '"grid_u"': '"' + gu_name + '"',
                        '"grid_d"': '"' + gd_name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gl_ud_replaceDic)
                    grid_l_conn = json.loads(file_content)

                if "f" in connections or "b" in connections:
                    gfBlock = connections["f"]
                    gbBlock = connections["b"]
                    gf_nx, gf_ny, gf_nz = gfBlock.nx, gfBlock.ny, gfBlock.nz
                    gb_nx, gb_ny, gb_nz = gbBlock.nx, gbBlock.ny, gbBlock.nz
                    gd_ny = connections["d"].ny
                    gl_fb_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gb_nx_m1_"': str(gb_nx - 1),
                        '"grid_f"': '"' + gfBlock.name + '"',
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gl_fb_replaceDic)
                    grid_l_conn = json.loads(file_content)
                
                grid_m["updateInstructions"].extend(grid_l_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_m["updateInstructions"].extend(grid_l_conn["updateInstructions"]["left_only"])
                elif "f" not in connections and "b" not in connections:
                    grid_m["updateInstructions"].extend(grid_l_conn["updateInstructions"]["left_up_down"])
                else:
                    grid_m["updateInstructions"].extend(grid_l_conn["updateInstructions"]["left_up_down_front_back"])
                
                for gl_updateSequence in grid_l_conn["updateSequences"]:
                    seqName = gl_updateSequence["name"]
                    for updateSequence in grid_m["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gl_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(gl_updateSequence["sequence"]["left_only"])
                            elif "f" not in connections and "b" not in connections:
                                updateSequence["sequence"].extend(gl_updateSequence["sequence"]["left_up_down"])
                            else:
                                updateSequence["sequence"].extend(gl_updateSequence["sequence"]["left_up_down_front_back"])

            #------------------- u ------------------
            #----------------------------------------
            if "u" in connections:
                json_file = open('layer0/grid_m/connections/grid_u.json')
                file_content = json_file.read()
                json_file.close()
                
                guBlock = connections["u"]
                gu_nx, gu_ny = guBlock.nx, guBlock.ny
                gl_nz = connections["l"].nz
                gu_replaceDic = {'"_nx_m1_"': str(nx - 1), '"_nx_m2_"': str(nx - 2),
                    '"_nz_m1_"': str(nz - 1), '"_nz_m2_"': str(nz - 2), 
                    '"_gl_nz_"': str(gl_nz), '"_gl_nz_m1_"': str(gl_nz - 1), '"_gl_nz_p1_"': str(gl_nz + 1), 
                    '"_gu_nx_m1_"': str(gu_nx - 1),
                    '"grid_u"': '"' + guBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, replaceDic)
                file_content = MultiWordReplace(file_content, gu_replaceDic)
                grid_u_conn = json.loads(file_content)
                
                if "f" in connections or "b" in connections:
                    gfBlock = connections["f"]
                    gbBlock = connections["b"]
                    gb_nx = gbBlock.nx
                    gd_ny = connections["d"].ny
                    gu_fb_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gd_ny_p_ny2_"': str(gd_ny + int(ny/2)),
                        '"_gb_nx_m1_"': str(gb_nx - 1),
                        '"grid_f"': '"' + gfBlock.name + '"',
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gu_fb_replaceDic)
                    grid_u_conn = json.loads(file_content)
                

                grid_m["updateInstructions"].extend(grid_u_conn["updateInstructions"]["general"])
                if "f" not in connections and "b" not in connections:
                    grid_m["updateInstructions"].extend(grid_u_conn["updateInstructions"]["up_only"])
                else:
                    grid_m["updateInstructions"].extend(grid_u_conn["updateInstructions"]["up_front_back"])
                
                                            
                for gu_updateSequence in grid_u_conn["updateSequences"]:
                    seqName = gu_updateSequence["name"]
                    for updateSequence in grid_m["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gu_updateSequence["sequence"]["general"])
                            if "f" not in connections and "b" not in connections:
                                updateSequence["sequence"].extend(gu_updateSequence["sequence"]["up_only"])
                            else:
                                updateSequence["sequence"].extend(gu_updateSequence["sequence"]["up_front_back"])

            #------------------- d ------------------
            #----------------------------------------
            if "d" in connections:
                json_file = open('layer0/grid_m/connections/grid_d.json')
                file_content = json_file.read()
                json_file.close()
                
                gdBlock = connections["d"]
                gd_nx, gd_ny = gdBlock.nx, gdBlock.ny
                gl_nz = connections["l"].nz
                gd_replaceDic = {
                    '"_nx_m1_"': str(nx - 1), '"_nx_m2_"': str(nx - 2),
                    '"_nz_m1_"': str(nz - 1),
                    '"_gl_nz_"': str(gl_nz), '"_gl_nz_m1_"': str(gl_nz - 1), '"_gl_nz_p1_"': str(gl_nz + 1),
                    '"_gd_ny_"': str(gd_ny), '"_gd_nx_m1_"': str(gd_nx - 1),
                    '"grid_d"': '"' + gdBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, replaceDic)
                file_content = MultiWordReplace(file_content, gd_replaceDic)
                grid_d_conn = json.loads(file_content)
                
                if "f" in connections or "b" in connections:
                    gfBlock = connections["f"]
                    gbBlock = connections["b"]
                    gb_nx = gbBlock.nx
                    gd_ny = connections["d"].ny
                    gd_fb_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gb_nx_m1_"': str(gb_nx - 1),
                        '"grid_f"': '"' + gfBlock.name + '"',
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gd_fb_replaceDic)
                    grid_d_conn = json.loads(file_content)
                

                grid_m["updateInstructions"].extend(grid_d_conn["updateInstructions"]["general"])
                if "f" not in connections and "b" not in connections:
                    grid_m["updateInstructions"].extend(grid_d_conn["updateInstructions"]["down_only"])
                else:
                    grid_m["updateInstructions"].extend(grid_d_conn["updateInstructions"]["down_front_back"])
                
                                            
                for gd_updateSequence in grid_d_conn["updateSequences"]:
                    seqName = gd_updateSequence["name"]
                    for updateSequence in grid_m["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gd_updateSequence["sequence"]["general"])
                            if "f" not in connections and "b" not in connections:
                                updateSequence["sequence"].extend(gd_updateSequence["sequence"]["down_only"])
                            else:
                                updateSequence["sequence"].extend(gd_updateSequence["sequence"]["down_front_back"])

            #------------------- f ------------------
            #----------------------------------------
            if "f" in connections:
                json_file = open('layer0/grid_m/connections/grid_f.json')
                file_content = json_file.read()
                json_file.close()
                
                gfBlock = connections["f"]
                gd_ny = connections["d"].ny
                gl_nz = connections["l"].nz
                gf_replaceDic = {'"_nx_m1_"': str(nx - 1), '"_nx_m2_"': str(nx - 2),
                    '"_ny_m1_"': str(ny - 1),
                    '"_nz_m1_"': str(nz - 1), '"_nz_m2_"': str(nz - 2), 
                    '"_gl_nz_"': str(gl_nz), '"_gl_nz_m1_"': str(gl_nz - 1), '"_gl_nz_p1_"': str(gl_nz + 1), 
                    '"_gd_ny_"': str(gd_ny), '"_gd_ny_p1_"': str(gd_ny + 1), '"_gd_ny_m1_"': str(gd_ny - 1),
                    '"grid_f"': '"' + gfBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, replaceDic)
                file_content = MultiWordReplace(file_content, gf_replaceDic)
                grid_f_conn = json.loads(file_content)
                
                grid_m["updateInstructions"].extend(grid_f_conn["updateInstructions"])                
                                            
                for gf_updateSequence in grid_f_conn["updateSequences"]:
                    seqName = gf_updateSequence["name"]
                    for updateSequence in grid_m["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gf_updateSequence["sequence"])

            #------------------- b ------------------
            #----------------------------------------
            if "b" in connections:
                json_file = open('layer0/grid_m/connections/grid_b.json')
                file_content = json_file.read()
                json_file.close()
                
                gbBlock = connections["b"]
                gb_nx = gbBlock.nx
                gd_ny = connections["d"].ny
                gl_nz = connections["l"].nz
                gb_replaceDic = {'"_nx_m1_"': str(nx - 1), '"_nx_m2_"': str(nx - 2),
                    '"_ny_m1_"': str(ny - 1),
                    '"_nz_m1_"': str(nz - 1), '"_nz_m2_"': str(nz - 2), 
                    '"_gl_nz_"': str(gl_nz), '"_gl_nz_m1_"': str(gl_nz - 1), '"_gl_nz_p1_"': str(gl_nz + 1), 
                    '"_gd_ny_"': str(gd_ny), '"_gd_ny_p1_"': str(gd_ny + 1), '"_gd_ny_m1_"': str(gd_ny - 1),
                    '"_gb_nx_"': str(gb_nx),
                    '"grid_b"': '"' + gbBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, replaceDic)
                file_content = MultiWordReplace(file_content, gb_replaceDic)
                grid_b_conn = json.loads(file_content)

                grid_m["updateInstructions"].extend(grid_b_conn["updateInstructions"])                
                                            
                for gb_updateSequence in grid_b_conn["updateSequences"]:
                    seqName = gb_updateSequence["name"]
                    for updateSequence in grid_m["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gb_updateSequence["sequence"])

            return grids

        ##---------------------------------------------- grid_r ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "r":
            if self.blockLevel == 1:
                json_file = open('layer1/grid_r/grid_r.json')
                file_content = json_file.read()
                json_file.close()
                
                gr_nx, gr_ny, gr_nz = self.nx, self.ny, self.nz
                gr_r0, gr_r1 = self.r0, self.r1
                gr_dt, gr_dx, gr_dy, gr_dz = self.dt, self.dx, self.dy, self.dz 
                gr_name = '"' + self.name + '"'
                assert "l" in self.connections
                assert self.connections["l"].blockLevel == self.blockLevel - 1
                gm_name = '"' + self.connections["l"].name + '"'
                nz = self.connections["l"].nz
                replaceDic = {'"_gr_nx_"': str(gr_nx), '"_gr_ny_"': str(gr_ny), '"_gr_nz_"': str(gr_nz),
                    '"_gr_nx_p1_"': str(gr_nx + 1), '"_gr_ny_p1_"': str(gr_ny + 1), '"_gr_nz_p1_"': str(gr_nz + 1),
                    '"_gr_x0_"': str(gr_r0[0]), '"_gr_y0_"': str(gr_r0[1]), '"_gr_z0_"': str(gr_r0[2]),
                    '"_gr_x1_"': str(gr_r1[0]), '"_gr_y1_"': str(gr_r1[1]), '"_gr_z1_"': str(gr_r1[2]),
                    '"_gr_dt_dx_"': str(gr_dt/gr_dx), '"_gr_dt_dy_"': str(gr_dt/gr_dy), '"_gr_dt_dz_"': str(gr_dt/gr_dz),
                    '"_gr_m_dt_dx_"': str(-gr_dt/gr_dx), '"_gr_m_dt_dy_"': str(-gr_dt/gr_dy), '"_gr_m_dt_dz_"': str(-gr_dt/gr_dz),
                    '"_gr_dt_dx_2_"': str(gr_dt/gr_dx/2.0), '"_gr_dt_dy_2_"': str(gr_dt/gr_dy/2.0), '"_gr_dt_dz_2_"': str(gr_dt/gr_dz/2.0),
                    '"_gr_m_dt_dx_2_"': str(-gr_dt/gr_dx/2.0), '"_gr_m_dt_dy_2_"': str(-gr_dt/gr_dy/2.0), '"_gr_m_dt_dz_2_"': str(-gr_dt/gr_dz/2.0),
                    '"_nz_m1_"': str(nz - 1), '"_nz_m2_"': str(nz - 2),
                    '"_gr_dt_"': str(gr_dt),
                    '"grid_r"': gr_name,
                    '"grid_m"': gm_name
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_r = grids[self.name]
                
                self.SetupSources(grid_r)
                self.SetupViews(grid_r)
                    
                connections = self.connections
                
                #------------------- u ------------------
                #----------------------------------------
                if "u" in connections:
                    json_file = open('layer1/grid_r/connections/grid_u.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    guBlock = connections["u"]
                    gl_nz = self.nz
                    gu_replaceDic = {
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                        '"grid_u"': '"' + guBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gu_replaceDic)
                    grid_u_conn = json.loads(file_content)
                    
                    grid_r["updateInstructions"].extend(grid_u_conn["updateInstructions"])
                    
                    for gu_updateSequence in grid_u_conn["updateSequences"]:
                        seqName = gu_updateSequence["name"]
                        for updateSequence in grid_r["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gu_updateSequence["sequence"])

                #------------------- d ------------------
                #----------------------------------------
                if "d" in connections:
                    json_file = open('layer1/grid_r/connections/grid_d.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gdBlock = connections["d"]
                    gd_ny = gdBlock.ny
                    gl_nz = self.nz
                    gd_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                        '"grid_d"': '"' + gdBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gd_replaceDic)
                    grid_d_conn = json.loads(file_content)
                    
                    grid_r["updateInstructions"].extend(grid_d_conn["updateInstructions"])
                    
                    for gd_updateSequence in grid_d_conn["updateSequences"]:
                        seqName = gd_updateSequence["name"]
                        for updateSequence in grid_r["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gd_updateSequence["sequence"])
                    

                #------------------- f ------------------
                #----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_r/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gfBlock = connections["f"]
                    gl_nz = self.nz
                    gd_ny = connections["d"].ny
                    gf_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                        '"grid_f"': '"' + gfBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gf_replaceDic)
                    grid_f_conn = json.loads(file_content)
                    
                    grid_r["updateInstructions"].extend(grid_f_conn["updateInstructions"])
                    
                    for gf_updateSequence in grid_f_conn["updateSequences"]:
                        seqName = gf_updateSequence["name"]
                        for updateSequence in grid_r["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gf_updateSequence["sequence"])

                #------------------- b ------------------
                #----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_r/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gbBlock = connections["b"]
                    gb_nx = gbBlock.nx
                    gl_nz = self.nz
                    gd_ny = connections["d"].ny
                    gb_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gb_nx_"': str(gb_nx),
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gb_replaceDic)
                    grid_b_conn = json.loads(file_content)
                    
                    grid_r["updateInstructions"].extend(grid_b_conn["updateInstructions"])
                    
                    for gb_updateSequence in grid_b_conn["updateSequences"]:
                        seqName = gb_updateSequence["name"]
                        for updateSequence in grid_r["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gb_updateSequence["sequence"])

                #------------------- rr ------------------
                #----------------------------------------
                if "r" in connections:
                    json_file = open('layer1/grid_r/connections/grid_rr.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    grrBlock = connections["r"]
                    gb_nx = connections["b"].nx
                    gd_ny = connections["d"].ny
                    grr_replaceDic = {
                        '"_gr_nx_m1_"': str(gr_nx - 1), '"_gr_ny_m1_"': str(gr_ny - 1),
                        '"_gd_ny2_"': str(int(gd_ny/2)), '"_gd_ny2_m1_"': str(int(gd_ny/2) - 1), '"_gd_ny2_p1_"': str(int(gd_ny/2) + 1),
                        '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                        '"grid_rr"': '"' + grrBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, grr_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                    
                    grid_r["updateInstructions"].extend(grid_rr_conn["updateInstructions"])
                    
                    for grr_updateSequence in grid_rr_conn["updateSequences"]:
                        seqName = grr_updateSequence["name"]
                        for updateSequence in grid_r["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"])

                ##---- return grid
                return grids

            ##---------------------------------------------- grid_rr ------------------------------------------------
            ##-------------------------------------------------------------------------------------------------------
            elif self.blockLevel >= 2:
                json_file = open('layer2/grid_rr/grid_rr.json')
                file_content = json_file.read()
                json_file.close()
                
                grr_nx, grr_ny, grr_nz = self.nx, self.ny, self.nz
                grr_r0, grr_r1 = self.r0, self.r1
                grr_dt, grr_dx, grr_dy, grr_dz = self.dt, self.dx, self.dy, self.dz 
                assert "lc" in self.connections
                assert "lu" in self.connections
                assert "ld" in self.connections
                assert "lf" in self.connections
                assert "lb" in self.connections
                assert self.connections["lc"].blockLevel == self.blockLevel - 1
                assert self.connections["lu"].blockLevel == self.blockLevel - 1
                assert self.connections["ld"].blockLevel == self.blockLevel - 1
                assert self.connections["lf"].blockLevel == self.blockLevel - 1
                assert self.connections["lb"].blockLevel == self.blockLevel - 1
                gr_nx, gr_ny, gr_nz = self.connections["lc"].nx, self.connections["lc"].ny, self.connections["lc"].nz
                gu_nz = self.connections["lu"].nz
                gd_ny, gd_nz = self.connections["ld"].ny, self.connections["ld"].nz
                gf_nz = self.connections["lf"].nz
                gb_nx, gb_nz = self.connections["lb"].nx, self.connections["lb"].nz
                replaceDic = {'"_grr_nx_"': str(grr_nx), '"_grr_ny_"': str(grr_ny), '"_grr_nz_"': str(grr_nz),
                    '"_grr_nx_p1_"': str(grr_nx + 1), '"_grr_ny_p1_"': str(grr_ny + 1), '"_grr_nz_p1_"': str(grr_nz + 1),
                    '"_grr_x0_"': str(grr_r0[0]), '"_grr_y0_"': str(grr_r0[1]), '"_grr_z0_"': str(grr_r0[2]),
                    '"_grr_x1_"': str(grr_r1[0]), '"_grr_y1_"': str(grr_r1[1]), '"_grr_z1_"': str(grr_r1[2]),
                    '"_grr_dt_dx_"': str(grr_dt/grr_dx), '"_grr_dt_dy_"': str(grr_dt/grr_dy), '"_grr_dt_dz_"': str(grr_dt/grr_dz),
                    '"_grr_m_dt_dx_"': str(-grr_dt/grr_dx), '"_grr_m_dt_dy_"': str(-grr_dt/grr_dy), '"_grr_m_dt_dz_"': str(-grr_dt/grr_dz),
                    '"_grr_dt_dx_2_"': str(grr_dt/grr_dx/2.0), '"_grr_dt_dy_2_"': str(grr_dt/grr_dy/2.0), '"_grr_dt_dz_2_"': str(grr_dt/grr_dz/2.0),
                    '"_grr_m_dt_dx_2_"': str(-grr_dt/grr_dx/2.0), '"_grr_m_dt_dy_2_"': str(-grr_dt/grr_dy/2.0), '"_grr_m_dt_dz_2_"': str(-grr_dt/grr_dz/2.0),
                    '"_gr_nz_m1_"': str(gr_nz - 1), '"_gr_nz_m2_"': str(gr_nz - 2),
                    '"_gu_nz_m1_"': str(gu_nz - 1), '"_gu_nz_m2_"': str(gu_nz - 2),
                    '"_gd_nz_m1_"': str(gd_nz - 1), '"_gd_nz_m2_"': str(gd_nz - 2),
                    '"_gf_nz_m1_"': str(gf_nz - 1), '"_gf_nz_m2_"': str(gf_nz - 2),
                    '"_gb_nz_m1_"': str(gb_nz - 1), '"_gb_nz_m2_"': str(gb_nz - 2),
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                    '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)),
                    '"_gd_ny2_"': str(int(gd_ny/2)),
                    '"_gd_ny2_p_gr_ny2_"': str(int(gd_ny/2) + int(gr_ny/2)),
                    '"_grr_dt_"': str(grr_dt),
                    '"grid_rr"': '"' + self.name + '"',
                    '"grid_r"': '"' + self.connections['lc'].name + '"',
                    '"grid_u"': '"' + self.connections['lu'].name + '"',
                    '"grid_d"': '"' + self.connections['ld'].name + '"',
                    '"grid_f"': '"' + self.connections['lf'].name + '"',
                    '"grid_b"': '"' + self.connections['lb'].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_rr = grids[self.name]
                
                self.SetupSources(grid_rr)
                self.SetupViews(grid_rr)
                    
                connections = self.connections

                ##---- grid_rr
                return grids
        
        ##---------------------------------------------- grid_l ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "l":
            if self.blockLevel == 1:
                json_file = open('layer1/grid_l/grid_l.json')
                file_content = json_file.read()
                json_file.close()
                
                gl_nx, gl_ny, gl_nz = self.nx, self.ny, self.nz
                gl_r0, gl_r1 = self.r0, self.r1
                gl_dt, gl_dx, gl_dy, gl_dz = self.dt, self.dx, self.dy, self.dz 
                gl_name = '"' + self.name + '"'
                assert "r" in self.connections
                assert self.connections["r"].blockLevel == self.blockLevel - 1
                gm_name = '"' + self.connections["r"].name + '"'
                replaceDic = {'"_gl_nx_"': str(gl_nx), '"_gl_ny_"': str(gl_ny), '"_gl_nz_"': str(gl_nz),
                    '"_gl_nx_p1_"': str(gl_nx + 1), '"_gl_ny_p1_"': str(gl_ny + 1), '"_gl_nz_p1_"': str(gl_nz + 1),
                    '"_gl_nz_m1_"': str(gl_nz - 1),
                    '"_gl_x0_"': str(gl_r0[0]), '"_gl_y0_"': str(gl_r0[1]), '"_gl_z0_"': str(gl_r0[2]),
                    '"_gl_x1_"': str(gl_r1[0]), '"_gl_y1_"': str(gl_r1[1]), '"_gl_z1_"': str(gl_r1[2]),
                    '"_gl_dt_dx_"': str(gl_dt/gl_dx), '"_gl_dt_dy_"': str(gl_dt/gl_dy), '"_gl_dt_dz_"': str(gl_dt/gl_dz),
                    '"_gl_m_dt_dx_"': str(-gl_dt/gl_dx), '"_gl_m_dt_dy_"': str(-gl_dt/gl_dy), '"_gl_m_dt_dz_"': str(-gl_dt/gl_dz),
                    '"_gl_dt_dx_2_"': str(gl_dt/gl_dx/2.0), '"_gl_dt_dy_2_"': str(gl_dt/gl_dy/2.0), '"_gl_dt_dz_2_"': str(gl_dt/gl_dz/2.0),
                    '"_gl_m_dt_dx_2_"': str(-gl_dt/gl_dx/2.0), '"_gl_m_dt_dy_2_"': str(-gl_dt/gl_dy/2.0), '"_gl_m_dt_dz_2_"': str(-gl_dt/gl_dz/2.0),
                    '"_gl_dt_"': str(gl_dt),
                    '"grid_l"': gl_name,
                    '"grid_m"': gm_name
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_l = grids[self.name]
                
                self.SetupSources(grid_l)
                self.SetupViews(grid_l)

                connections = self.connections
                
                #------------------- u ------------------
                #----------------------------------------
                if "u" in connections:
                    json_file = open('layer1/grid_l/connections/grid_u.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    guBlock = connections["u"]
                    gl_nz = self.nz
                    gu_replaceDic = {
                        '"grid_u"': '"' + guBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gu_replaceDic)
                    grid_u_conn = json.loads(file_content)
                    
                    grid_l["updateInstructions"].extend(grid_u_conn["updateInstructions"])
                    
                    for gu_updateSequence in grid_u_conn["updateSequences"]:
                        seqName = gu_updateSequence["name"]
                        for updateSequence in grid_l["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gu_updateSequence["sequence"])

                #------------------- d ------------------
                #----------------------------------------
                if "d" in connections:
                    json_file = open('layer1/grid_l/connections/grid_d.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gdBlock = connections["d"]
                    gd_ny = gdBlock.ny
                    gl_nz = self.nz
                    gd_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"grid_d"': '"' + gdBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gd_replaceDic)
                    grid_d_conn = json.loads(file_content)
                    
                    grid_l["updateInstructions"].extend(grid_d_conn["updateInstructions"])
                    
                    for gd_updateSequence in grid_d_conn["updateSequences"]:
                        seqName = gd_updateSequence["name"]
                        for updateSequence in grid_l["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gd_updateSequence["sequence"])
                    
                #------------------- f ------------------
                #----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_l/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gfBlock = connections["f"]
                    gd_ny = connections["d"].ny
                    gf_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"grid_f"': '"' + gfBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gf_replaceDic)
                    grid_f_conn = json.loads(file_content)
                    
                    grid_l["updateInstructions"].extend(grid_f_conn["updateInstructions"])
                    
                    for gf_updateSequence in grid_f_conn["updateSequences"]:
                        seqName = gf_updateSequence["name"]
                        for updateSequence in grid_l["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gf_updateSequence["sequence"])

                #------------------- b ------------------
                #----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_l/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gbBlock = connections["b"]
                    gb_nx = gbBlock.nx
                    gd_ny = connections["d"].ny
                    gb_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gb_nx_"': str(gb_nx),
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gb_replaceDic)
                    grid_b_conn = json.loads(file_content)
                    
                    grid_l["updateInstructions"].extend(grid_b_conn["updateInstructions"])
                    
                    for gb_updateSequence in grid_b_conn["updateSequences"]:
                        seqName = gb_updateSequence["name"]
                        for updateSequence in grid_l["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gb_updateSequence["sequence"])

                #------------------- ll ------------------
                #----------------------------------------
                if "l" in connections:
                    json_file = open('layer1/grid_l/connections/grid_ll.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gllBlock = connections["l"]
                    gll_nz = gllBlock.nz
                    gb_nx = connections["b"].nx
                    gd_ny = connections["d"].ny
                    gll_replaceDic = {
                        '"_gl_nx_m1_"': str(gl_nx - 1), '"_gl_ny_m1_"': str(gl_ny - 1),
                        '"_gll_nz_"': str(gll_nz),
                        '"_gd_ny2_"': str(int(gd_ny/2)), '"_gd_ny2_m1_"': str(int(gd_ny/2) - 1), '"_gd_ny2_p1_"': str(int(gd_ny/2) + 1),
                        '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                        '"grid_ll"': '"' + gllBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gll_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                    
                    grid_l["updateInstructions"].extend(grid_ll_conn["updateInstructions"])
                    
                    for gll_updateSequence in grid_ll_conn["updateSequences"]:
                        seqName = gll_updateSequence["name"]
                        for updateSequence in grid_l["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"])

                ##--- grid_l --
                return grids
        
            ##---------------------------------------------- grid_ll ------------------------------------------------
            ##------------------------------------------------------------------------------------------------------
            if self.blockLevel >= 2:
                json_file = open('layer2/grid_ll/grid_ll.json')
                file_content = json_file.read()
                json_file.close()
                
                gll_nx, gll_ny, gll_nz = self.nx, self.ny, self.nz
                gll_r0, gll_r1 = self.r0, self.r1
                gll_dt, gll_dx, gll_dy, gll_dz = self.dt, self.dx, self.dy, self.dz 
                assert "rc" in self.connections
                assert "ru" in self.connections
                assert "rd" in self.connections
                assert "rf" in self.connections
                assert "rb" in self.connections
                assert self.connections["rc"].blockLevel == self.blockLevel - 1
                assert self.connections["ru"].blockLevel == self.blockLevel - 1
                assert self.connections["rd"].blockLevel == self.blockLevel - 1
                assert self.connections["rf"].blockLevel == self.blockLevel - 1
                assert self.connections["rb"].blockLevel == self.blockLevel - 1
                gr_nx = self.connections["rc"].nx   # gr_nx = gl_nx
                gr_ny = self.connections["rc"].ny
                gb_nx = self.connections["rb"].nx
                gd_ny = self.connections["rd"].ny
                replaceDic = {'"_gll_nx_"': str(gll_nx), '"_gll_ny_"': str(gll_ny), '"_gll_nz_"': str(gll_nz),
                    '"_gll_nx_p1_"': str(gll_nx + 1), '"_gll_ny_p1_"': str(gll_ny + 1), '"_gll_nz_p1_"': str(gll_nz + 1),
                    '"_gll_nz_m1_"': str(gll_nz - 1),
                    '"_gll_x0_"': str(gll_r0[0]), '"_gll_y0_"': str(gll_r0[1]), '"_gll_z0_"': str(gll_r0[2]),
                    '"_gll_x1_"': str(gll_r1[0]), '"_gll_y1_"': str(gll_r1[1]), '"_gll_z1_"': str(gll_r1[2]),
                    '"_gll_dt_dx_"': str(gll_dt/gll_dx), '"_gll_dt_dy_"': str(gll_dt/gll_dy), '"_gll_dt_dz_"': str(gll_dt/gll_dz),
                    '"_gll_m_dt_dx_"': str(-gll_dt/gll_dx), '"_gll_m_dt_dy_"': str(-gll_dt/gll_dy), '"_gll_m_dt_dz_"': str(-gll_dt/gll_dz),
                    '"_gll_dt_dx_2_"': str(gll_dt/gll_dx/2.0), '"_gll_dt_dy_2_"': str(gll_dt/gll_dy/2.0), '"_gll_dt_dz_2_"': str(gll_dt/gll_dz/2.0),
                    '"_gll_m_dt_dx_2_"': str(-gll_dt/gll_dx/2.0), '"_gll_m_dt_dy_2_"': str(-gll_dt/gll_dy/2.0), '"_gll_m_dt_dz_2_"': str(-gll_dt/gll_dz/2.0),
                    '"_gll_dt_"': str(gll_dt),
                    '"_gb_nx2_"': str(int(gb_nx/2)),
                    '"_gd_ny2_"': str(int(gd_ny/2)),
                    '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)),
                    '"_gd_ny2_p_gr_ny2_"': str(int(gd_ny/2) + int(gr_ny/2)),
                    '"_gd_ny2_p_gr_ny2_p1_"': str(int(gd_ny/2) + int(gr_ny/2) + 1),
                    '"grid_ll"': '"' + self.name + '"',
                    '"grid_l"': '"' + self.connections['rc'].name + '"',
                    '"grid_u"': '"' + self.connections['ru'].name + '"',
                    '"grid_d"': '"' + self.connections['rd'].name + '"',
                    '"grid_f"': '"' + self.connections['rf'].name + '"',
                    '"grid_b"': '"' + self.connections['rb'].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_ll = grids[self.name]
                
                self.SetupSources(grid_ll)
                self.SetupViews(grid_ll)

                connections = self.connections

                ##--- grid_ll --
                return grids

        ##---------------------------------------------- grid_u ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "u":
            if self.blockLevel == 1:
                json_file = open('layer1/grid_u/grid_u.json')
                file_content = json_file.read()
                json_file.close()
                
                gu_nx, gu_ny, gu_nz = self.nx, self.ny, self.nz
                gu_r0, gu_r1 = self.r0, self.r1
                gu_dt, gu_dx, gu_dy, gu_dz = self.dt, self.dx, self.dy, self.dz 
                assert "dc" in self.connections
                assert "dl" in self.connections
                assert "dr" in self.connections
                assert self.connections["dc"].blockLevel == self.blockLevel - 1
                assert self.connections["dl"].blockLevel == self.blockLevel
                assert self.connections["dr"].blockLevel == self.blockLevel
                gr_ny = self.connections["dr"].ny
                gl_ny, gl_nz = self.connections["dl"].ny, self.connections["dl"].nz
                ny, nz = self.connections["dc"].ny, self.connections["dc"].nz
        
                replaceDic = {'"_gu_nx_"': str(gu_nx), '"_gu_ny_"': str(gu_ny), '"_gu_nz_"': str(gu_nz),
                    '"_gu_nx_p1_"': str(gu_nx + 1), '"_gu_ny_p1_"': str(gu_ny + 1), '"_gu_nz_p1_"': str(gu_nz + 1),
                    '"_gu_nz_m1_"': str(gu_nz - 1),
                    '"_gu_x0_"': str(gu_r0[0]), '"_gu_y0_"': str(gu_r0[1]), '"_gu_z0_"': str(gu_r0[2]),
                    '"_gu_x1_"': str(gu_r1[0]), '"_gu_y1_"': str(gu_r1[1]), '"_gu_z1_"': str(gu_r1[2]),
                    '"_gu_dt_dx_"': str(gu_dt/gu_dx), '"_gu_dt_dy_"': str(gu_dt/gu_dy), '"_gu_dt_dz_"': str(gu_dt/gu_dz),
                    '"_gu_m_dt_dx_"': str(-gu_dt/gu_dx), '"_gu_m_dt_dy_"': str(-gu_dt/gu_dy), '"_gu_m_dt_dz_"': str(-gu_dt/gu_dz),
                    '"_gu_dt_dx_2_"': str(gu_dt/gu_dx/2.0), '"_gu_dt_dy_2_"': str(gu_dt/gu_dy/2.0), '"_gu_dt_dz_2_"': str(gu_dt/gu_dz/2.0),
                    '"_gu_m_dt_dx_2_"': str(-gu_dt/gu_dx/2.0), '"_gu_m_dt_dy_2_"': str(-gu_dt/gu_dy/2.0), '"_gu_m_dt_dz_2_"': str(-gu_dt/gu_dz/2.0),
                    '"_gu_dt_"': str(gu_dt),
                    '"_ny_m1_"': str(ny - 1), '"_ny_m2_"': str(ny - 2),
                    '"_gl_nz_"': str(gl_nz),
                    '"_gl_ny_m1_"': str(gl_ny - 1), '"_gl_nz_p1_"': str(gl_nz + 1), 
                    '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                    '"_gr_ny_m1_"': str(gr_ny - 1), 
                    '"grid_u"': '"' + self.name + '"',
                    '"grid_r"': '"' + self.connections["dr"].name + '"',
                    '"grid_l"': '"' + self.connections["dl"].name + '"',
                    '"grid_m"': '"' + self.connections["dc"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_u = grids[self.name]
                
                self.SetupSources(grid_u)
                self.SetupViews(grid_u)
                
                connections = self.connections

                #------------------- f ------------------
                #----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_u/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gfBlock = connections["f"]
                    gd_ny = self.ny
                    gf_replaceDic = {
                        '"_gd_ny_p_ny2_"': str(gd_ny + int(ny/2)),
                        '"grid_f"': '"' + gfBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gf_replaceDic)
                    grid_f_conn = json.loads(file_content)
                    
                    grid_u["updateInstructions"].extend(grid_f_conn["updateInstructions"])
                    
                    for gf_updateSequence in grid_f_conn["updateSequences"]:
                        seqName = gf_updateSequence["name"]
                        for updateSequence in grid_u["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gf_updateSequence["sequence"])

                #------------------- b ------------------
                #----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_u/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gbBlock = connections["b"]
                    gb_nx = gbBlock.nx
                    gd_ny = self.ny
                    gb_replaceDic = {
                        '"_gb_nx_"': str(gb_nx),
                        '"_gd_ny_p_ny2_"': str(gd_ny + int(ny/2)),
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gb_replaceDic)
                    grid_b_conn = json.loads(file_content)
                    
                    grid_u["updateInstructions"].extend(grid_b_conn["updateInstructions"])
                    
                    for gb_updateSequence in grid_b_conn["updateSequences"]:
                        seqName = gb_updateSequence["name"]
                        for updateSequence in grid_u["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gb_updateSequence["sequence"])

                #------------------- rr ------------------
                #----------------------------------------
                if "r" in connections:
                    json_file = open('layer1/grid_u/connections/grid_rr.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    grrBlock = connections["r"]
                    grr_ny = grrBlock.ny
                    gb_nx = connections["b"].nx
                    gd_ny = self.ny
                    gr_ny = connections["dr"].ny
                    grr_replaceDic = {
                        '"_gu_nx_m1_"': str(gu_nx - 1), '"_gu_ny_m1_"': str(gu_ny - 1), '"_gu_ny_m2_"': str(gu_ny - 2),
                        '"_grr_ny_m1_"': str(grr_ny - 1), 
                        '"_gd_ny2_p_gr_ny2_"': str(int(gd_ny/2) + int(gr_ny/2)), 
                        '"_gd_ny2_p_gr_ny2_m1_"': str(int(gd_ny/2) + int(gr_ny/2) - 1), 
                        '"_gd_ny2_p_gr_ny2_p1_"': str(int(gd_ny/2) + int(gr_ny/2) + 1), 
                        '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                        '"grid_rr"': '"' + grrBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, grr_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                    
                    grid_u["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                    if "u" not in connections:
                        grid_u["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                    else:
                        grid_u["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up"])
                    
                    
                    for grr_updateSequence in grid_rr_conn["updateSequences"]:
                        seqName = grr_updateSequence["name"]
                        for updateSequence in grid_u["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                                if "u" not in connections:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                                else:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up"])                                

                #------------------- ll ------------------
                #----------------------------------------
                if "l" in connections:
                    json_file = open('layer1/grid_u/connections/grid_ll.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gllBlock = connections["l"]
                    gll_ny, gll_nz = gllBlock.ny, gllBlock.nz
                    gb_nx = connections["b"].nx
                    gd_ny = self.ny
                    gr_ny = connections["dr"].ny
                    gll_replaceDic = {
                        '"_gu_nx_m1_"': str(gu_nx - 1), '"_gu_ny_m1_"': str(gu_ny - 1), '"_gu_ny_m2_"': str(gu_ny - 2),
                        '"_gll_ny_m1_"': str(gll_ny - 1), 
                        '"_gll_nz_"': str(gll_nz),
                        '"_gd_ny2_p_gr_ny2_"': str(int(gd_ny/2) + int(gr_ny/2)), 
                        '"_gd_ny2_p_gr_ny2_m1_"': str(int(gd_ny/2) + int(gr_ny/2) - 1), 
                        '"_gd_ny2_p_gr_ny2_p1_"': str(int(gd_ny/2) + int(gr_ny/2) + 1), 
                        '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                        '"grid_ll"': '"' + gllBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gll_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                    
                    grid_u["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                    if "u" not in connections:
                        grid_u["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                    else:
                        grid_u["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up"])
                    
                    for gll_updateSequence in grid_ll_conn["updateSequences"]:
                        seqName = gll_updateSequence["name"]
                        for updateSequence in grid_u["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                                if "u" not in connections:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                                else:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up"])                                
                
                ##------         
                return grids
        
        ##---------------------------------------------- grid_d ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "d":
            if self.blockLevel == 1:
                json_file = open('layer1/grid_d/grid_d.json')
                file_content = json_file.read()
                json_file.close()
                
                gd_nx, gd_ny, gd_nz = self.nx, self.ny, self.nz
                gd_r0, gd_r1 = self.r0, self.r1
                gd_dt, gd_dx, gd_dy, gd_dz = self.dt, self.dx, self.dy, self.dz 
                assert "uc" in self.connections
                assert "ul" in self.connections
                assert "ur" in self.connections
                assert self.connections["uc"].blockLevel == self.blockLevel - 1
                assert self.connections["ul"].blockLevel == self.blockLevel
                assert self.connections["ur"].blockLevel == self.blockLevel
                
                gl_nz = self.connections["ul"].nz
                nz = self.connections["uc"].nz
        
                replaceDic = {'"_gd_nx_"': str(gd_nx), '"_gd_ny_"': str(gd_ny), '"_gd_nz_"': str(gd_nz),
                    '"_gd_nx_p1_"': str(gd_nx + 1), '"_gd_ny_p1_"': str(gd_ny + 1), '"_gd_nz_p1_"': str(gd_nz + 1),
                    '"_gd_ny_m1_"': str(gd_ny - 1), '"_gd_nz_m1_"': str(gd_nz - 1),
                    '"_gd_x0_"': str(gd_r0[0]), '"_gd_y0_"': str(gd_r0[1]), '"_gd_z0_"': str(gd_r0[2]),
                    '"_gd_x1_"': str(gd_r1[0]), '"_gd_y1_"': str(gd_r1[1]), '"_gd_z1_"': str(gd_r1[2]),
                    '"_gd_dt_dx_"': str(gd_dt/gd_dx), '"_gd_dt_dy_"': str(gd_dt/gd_dy), '"_gd_dt_dz_"': str(gd_dt/gd_dz),
                    '"_gd_m_dt_dx_"': str(-gd_dt/gd_dx), '"_gd_m_dt_dy_"': str(-gd_dt/gd_dy), '"_gd_m_dt_dz_"': str(-gd_dt/gd_dz),
                    '"_gd_dt_dx_2_"': str(gd_dt/gd_dx/2.0), '"_gd_dt_dy_2_"': str(gd_dt/gd_dy/2.0), '"_gd_dt_dz_2_"': str(gd_dt/gd_dz/2.0),
                    '"_gd_m_dt_dx_2_"': str(-gd_dt/gd_dx/2.0), '"_gd_m_dt_dy_2_"': str(-gd_dt/gd_dy/2.0), '"_gd_m_dt_dz_2_"': str(-gd_dt/gd_dz/2.0),
                    '"_gd_dt_"': str(gd_dt),
                    '"_gl_nz_"': str(gl_nz), '"_gl_nz_p1_"': str(gl_nz + 1), '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                    '"grid_d"': '"' + self.name + '"',
                    '"grid_r"': '"' + self.connections["ur"].name + '"',
                    '"grid_l"': '"' + self.connections["ul"].name + '"',
                    '"grid_m"': '"' + self.connections["uc"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_d = grids[self.name]
                
                self.SetupSources(grid_d)
                self.SetupViews(grid_d)

                connections = self.connections

                #------------------- f ------------------
                #----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_d/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gfBlock = connections["f"]
                    gf_replaceDic = {
                        '"grid_f"': '"' + gfBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gf_replaceDic)
                    grid_f_conn = json.loads(file_content)
                    
                    grid_d["updateInstructions"].extend(grid_f_conn["updateInstructions"])
                    
                    for gf_updateSequence in grid_f_conn["updateSequences"]:
                        seqName = gf_updateSequence["name"]
                        for updateSequence in grid_d["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gf_updateSequence["sequence"])

                #------------------- b ------------------
                #----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_d/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gbBlock = connections["b"]
                    gb_nx = gbBlock.nx
                    gb_replaceDic = {
                        '"_gb_nx_"': str(gb_nx),
                        '"grid_b"': '"' + gbBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gb_replaceDic)
                    grid_b_conn = json.loads(file_content)
                    
                    grid_d["updateInstructions"].extend(grid_b_conn["updateInstructions"])
                    
                    for gb_updateSequence in grid_b_conn["updateSequences"]:
                        seqName = gb_updateSequence["name"]
                        for updateSequence in grid_d["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gb_updateSequence["sequence"])
        
                #------------------- rr ------------------
                #----------------------------------------
                if "r" in connections:
                    json_file = open('layer1/grid_d/connections/grid_rr.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    grrBlock = connections["r"]
                    gb_nx = connections["b"].nx
                    grr_replaceDic = {
                        '"_gd_nx_m1_"': str(gd_nx - 1),
                        '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                        '"grid_rr"': '"' + grrBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, grr_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                    
                    grid_d["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                    if "d" not in connections:
                        grid_d["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                    else:
                        grid_d["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_down"])
                    
                    for grr_updateSequence in grid_rr_conn["updateSequences"]:
                        seqName = grr_updateSequence["name"]
                        for updateSequence in grid_d["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                                if "d" not in connections:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                                else:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_down"])                                
                                
                                
                #------------------- ll ------------------
                #----------------------------------------
                if "l" in connections:
                    json_file = open('layer1/grid_d/connections/grid_ll.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gllBlock = connections["l"]
                    gll_ny, gll_nz = gllBlock.ny, gllBlock.nz
                    gb_nx = connections["b"].nx
                    gll_replaceDic = {
                        '"_gd_nx_m1_"': str(gd_nx - 1),
                        '"_gll_ny_m1_"': str(gll_ny - 1), 
                        '"_gll_nz_"': str(gll_nz),
                        '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                        '"grid_ll"': '"' + gllBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gll_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                    
                    grid_d["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                    if "d" not in connections:
                        grid_d["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                    else:
                        grid_d["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_down"])
                    
                    for gll_updateSequence in grid_ll_conn["updateSequences"]:
                        seqName = gll_updateSequence["name"]
                        for updateSequence in grid_d["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                                if "d" not in connections:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                                else:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_down"])                                
                
                ##------         
                return grids
        
        ##---------------------------------------------- grid_f ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "f":
            if self.blockLevel == 1:
                json_file = open('layer1/grid_f/grid_f.json')
                file_content = json_file.read()
                json_file.close()
                
                gf_nx, gf_ny, gf_nz = self.nx, self.ny, self.nz
                gf_r0, gf_r1 = self.r0, self.r1
                gf_dt, gf_dx, gf_dy, gf_dz = self.dt, self.dx, self.dy, self.dz 
                assert "bc" in self.connections
                assert "bl" in self.connections
                assert "br" in self.connections
                assert "bu" in self.connections
                assert "bd" in self.connections
                assert self.connections["bc"].blockLevel == self.blockLevel - 1
                assert self.connections["bl"].blockLevel == self.blockLevel
                assert self.connections["br"].blockLevel == self.blockLevel
                assert self.connections["bu"].blockLevel == self.blockLevel
                assert self.connections["bd"].blockLevel == self.blockLevel

                gmBlock = self.connections["bc"]
                nx, ny, nz = gmBlock.nx, gmBlock.ny, gmBlock.nz 
                gl_nx, gl_nz = self.connections["bl"].nx, self.connections["bl"].nz
                gd_nx, gd_ny = self.connections["bd"].nx, self.connections["bd"].ny
                gr_nx = self.connections["br"].nx
                gu_nx = self.connections["bu"].nx
        
                replaceDic = {'"_gf_nx_"': str(gf_nx), '"_gf_ny_"': str(gf_ny), '"_gf_nz_"': str(gf_nz),
                    '"_gf_nx_p1_"': str(gf_nx + 1), '"_gf_ny_p1_"': str(gf_ny + 1), '"_gf_nz_p1_"': str(gf_nz + 1),
                    '"_gf_nz_m1_"': str(gf_nz - 1),
                    '"_gf_x0_"': str(gf_r0[0]), '"_gf_y0_"': str(gf_r0[1]), '"_gf_z0_"': str(gf_r0[2]),
                    '"_gf_x1_"': str(gf_r1[0]), '"_gf_y1_"': str(gf_r1[1]), '"_gf_z1_"': str(gf_r1[2]),
                    '"_gf_dt_dx_"': str(gf_dt/gf_dx), '"_gf_dt_dy_"': str(gf_dt/gf_dy), '"_gf_dt_dz_"': str(gf_dt/gf_dz),
                    '"_gf_m_dt_dx_"': str(-gf_dt/gf_dx), '"_gf_m_dt_dy_"': str(-gf_dt/gf_dy), '"_gf_m_dt_dz_"': str(-gf_dt/gf_dz),
                    '"_gf_dt_dx_2_"': str(gf_dt/gf_dx/2.0), '"_gf_dt_dy_2_"': str(gf_dt/gf_dy/2.0), '"_gf_dt_dz_2_"': str(gf_dt/gf_dz/2.0),
                    '"_gf_m_dt_dx_2_"': str(-gf_dt/gf_dx/2.0), '"_gf_m_dt_dy_2_"': str(-gf_dt/gf_dy/2.0), '"_gf_m_dt_dz_2_"': str(-gf_dt/gf_dz/2.0),
                    '"_gf_dt_"': str(gf_dt),
                    '"_gl_nz_"': str(gl_nz), '"_gl_nz_p1_"': str(gl_nz + 1), '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                    '"_gd_ny_"': str(gd_ny), '"_gd_ny_p1_"': str(gd_ny + 1), '"_gd_ny_p_ny2_"': str(gd_ny + int(ny/2)),
                    '"_nx_m1_"': str(nx - 1), '"_nx_m2_"': str(nx - 2),
                    '"_gl_nx_m1_"': str(gl_nx - 1), 
                    '"_gr_nx_m1_"': str(gr_nx - 1), 
                    '"_gu_nx_m1_"': str(gu_nx - 1),
                    '"_gd_nx_m1_"': str(gd_nx - 1),
                    '"grid_f"': '"' + self.name + '"',
                    '"grid_u"': '"' + self.connections["bu"].name + '"',
                    '"grid_d"': '"' + self.connections["bd"].name + '"',
                    '"grid_r"': '"' + self.connections["br"].name + '"',
                    '"grid_l"': '"' + self.connections["bl"].name + '"',
                    '"grid_m"': '"' + self.connections["bc"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_f = grids[self.name]
                
                self.SetupSources(grid_f)
                self.SetupViews(grid_f)
        
                connections = self.connections

                #------------------- rr ------------------
                #----------------------------------------
                if "r" in connections:
                    json_file = open('layer1/grid_f/connections/grid_rr.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    grrBlock = connections["r"]
                    grr_nx, grr_ny = grrBlock.nx, grrBlock.ny
                    gb_nx = self.nx
                    gr_nx = connections["br"].nx
                    grr_replaceDic = {
                        '"_gf_nx_m1_"': str(gf_nx - 1), '"_gf_nx_m2_"': str(gf_nx - 2),
                        '"_gf_ny_m1_"': str(gf_ny - 1), '"_gf_ny_m2_"': str(gf_ny - 2),
                        '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)),
                        '"_gb_nx2_p_gr_nx2_p1_"': str(int(gb_nx/2) + int(gr_nx/2) + 1),
                        '"_gb_nx2_p_gr_nx2_m1_"': str(int(gb_nx/2) + int(gr_nx/2) - 1),
                        '"_grr_nx_m1_"': str(grr_nx - 1),
                        '"_grr_ny_m1_"': str(grr_ny - 1),
                        '"grid_rr"': '"' + grrBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, grr_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                    
                    grid_f["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                    if "u" not in connections and "d" not in connections:
                        grid_f["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                    elif "f" not in connections:
                        grid_f["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down"])
                    else:
                        grid_f["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down_front"])
                    
                    for grr_updateSequence in grid_rr_conn["updateSequences"]:
                        seqName = grr_updateSequence["name"]
                        for updateSequence in grid_f["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                                if "u" not in connections and "d" not in connections:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                                elif "f" not in connections:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down"])  
                                else:                              
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down_front"])  
                                
                                
                #------------------- ll ------------------
                #----------------------------------------
                if "l" in connections:
                    json_file = open('layer1/grid_f/connections/grid_ll.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gllBlock = connections["l"]
                    gll_nx, gll_ny, gll_nz = gllBlock.nx, gllBlock.ny, gllBlock.nz
                    gb_nx = self.nx
                    gll_replaceDic = {
                        '"_gf_nx_m1_"': str(gf_nx - 1), '"_gf_nx_m2_"': str(gf_nx - 2),
                        '"_gf_ny_m1_"': str(gf_ny - 1), '"_gf_ny_m2_"': str(gf_ny - 2),
                        '"_gll_nx_m1_"': str(gll_nx - 1), 
                        '"_gll_ny_m1_"': str(gll_ny - 1), 
                        '"_gll_nz_"': str(gll_nz),
                        '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)),
                        '"_gb_nx2_p_gr_nx2_p1_"': str(int(gb_nx/2) + int(gr_nx/2) + 1),
                        '"_gb_nx2_p_gr_nx2_m1_"': str(int(gb_nx/2) + int(gr_nx/2) - 1),
                        '"grid_ll"': '"' + gllBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gll_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                    
                    grid_f["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                    if "u" not in connections and "d" not in connections:
                        grid_f["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                    elif "f" not in connections:
                        grid_f["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down"])
                    else:
                        grid_f["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down_front"])
                    
                    for gll_updateSequence in grid_ll_conn["updateSequences"]:
                        seqName = gll_updateSequence["name"]
                        for updateSequence in grid_f["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                                if "u" not in connections and "d" not in connections:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                                elif "f" not in connections:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down"])
                                else:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down_front"])
                
                ##----
                return grids

        ##---------------------------------------------- grid_b ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "b":
            if self.blockLevel == 1:
                json_file = open('layer1/grid_b/grid_b.json')
                file_content = json_file.read()
                json_file.close()
                
                gb_nx, gb_ny, gb_nz = self.nx, self.ny, self.nz
                gb_r0, gb_r1 = self.r0, self.r1
                gb_dt, gb_dx, gb_dy, gb_dz = self.dt, self.dx, self.dy, self.dz 
                assert "fc" in self.connections
                assert "fl" in self.connections
                assert "fr" in self.connections
                assert "fu" in self.connections
                assert "fd" in self.connections
                assert self.connections["fc"].blockLevel == self.blockLevel - 1
                assert self.connections["fl"].blockLevel == self.blockLevel
                assert self.connections["fr"].blockLevel == self.blockLevel
                assert self.connections["fu"].blockLevel == self.blockLevel
                assert self.connections["fd"].blockLevel == self.blockLevel
                
                ny, nz = self.connections["fc"].ny, self.connections["fc"].nz 
                gl_nz = self.connections["fl"].nz
                gd_ny = self.connections["fd"].ny
        
                replaceDic = {'"_gb_nx_"': str(gb_nx), '"_gb_ny_"': str(gb_ny), '"_gb_nz_"': str(gb_nz),
                    '"_gb_nx_p1_"': str(gb_nx + 1), '"_gb_ny_p1_"': str(gb_ny + 1), '"_gb_nz_p1_"': str(gb_nz + 1),
                    '"_gb_nx_m1_"': str(gb_nx - 1), '"_gb_nz_m1_"': str(gb_nz - 1),
                    '"_gb_x0_"': str(gb_r0[0]), '"_gb_y0_"': str(gb_r0[1]), '"_gb_z0_"': str(gb_r0[2]),
                    '"_gb_x1_"': str(gb_r1[0]), '"_gb_y1_"': str(gb_r1[1]), '"_gb_z1_"': str(gb_r1[2]),
                    '"_gb_dt_dx_"': str(gb_dt/gb_dx), '"_gb_dt_dy_"': str(gb_dt/gb_dy), '"_gb_dt_dz_"': str(gb_dt/gb_dz),
                    '"_gb_m_dt_dx_"': str(-gb_dt/gb_dx), '"_gb_m_dt_dy_"': str(-gb_dt/gb_dy), '"_gb_m_dt_dz_"': str(-gb_dt/gb_dz),
                    '"_gb_dt_dx_2_"': str(gb_dt/gb_dx/2.0), '"_gb_dt_dy_2_"': str(gb_dt/gb_dy/2.0), '"_gb_dt_dz_2_"': str(gb_dt/gb_dz/2.0),
                    '"_gb_m_dt_dx_2_"': str(-gb_dt/gb_dx/2.0), '"_gb_m_dt_dy_2_"': str(-gb_dt/gb_dy/2.0), '"_gb_m_dt_dz_2_"': str(-gb_dt/gb_dz/2.0),
                    '"_gb_dt_"': str(gb_dt),
                    '"_gl_nz_"': str(gl_nz), '"_gl_nz_p1_"': str(gl_nz + 1), '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                    '"_gd_ny_"': str(gd_ny), '"_gd_ny_p1_"': str(gd_ny + 1), '"_gd_ny_p_ny2_"': str(gd_ny + int(ny/2)),
                    '"grid_b"': '"' + self.name + '"',
                    '"grid_u"': '"' + self.connections["fu"].name + '"',
                    '"grid_d"': '"' + self.connections["fd"].name + '"',
                    '"grid_r"': '"' + self.connections["fr"].name + '"',
                    '"grid_l"': '"' + self.connections["fl"].name + '"',
                    '"grid_m"': '"' + self.connections["fc"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_b = grids[self.name]
                
                self.SetupSources(grid_b)
                self.SetupViews(grid_b)
        
                connections = self.connections

                #------------------- rr ------------------
                #----------------------------------------
                if "r" in connections:
                    json_file = open('layer1/grid_b/connections/grid_rr.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    grrBlock = connections["r"]
                    grr_ny = grrBlock.ny
                    grr_replaceDic = {
                        '"_gb_ny_m1_"': str(gb_ny - 1), '"_gb_ny_m2_"': str(gb_ny - 2),
                        '"_grr_ny_m1_"': str(grr_ny - 1),
                        '"grid_rr"': '"' + grrBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, grr_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                    
                    grid_b["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                    if "u" not in connections and "d" not in connections:
                        grid_b["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                    elif "b" not in connections:
                        grid_b["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down"])
                    else:
                        grid_b["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down_back"])
                    
                    for grr_updateSequence in grid_rr_conn["updateSequences"]:
                        seqName = grr_updateSequence["name"]
                        for updateSequence in grid_b["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                                if "u" not in connections and "d" not in connections:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                                elif "b" not in connections:
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down"])  
                                else:                              
                                    updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down_back"])  
                                
                                
                #------------------- ll ------------------
                #-----------------------------------------
                if "l" in connections:
                    json_file = open('layer1/grid_b/connections/grid_ll.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gllBlock = connections["l"]
                    gll_ny, gll_nz = gllBlock.ny, gllBlock.nz
                    gll_replaceDic = {
                        '"_gb_ny_m1_"': str(gb_ny - 1), '"_gb_ny_m2_"': str(gb_ny - 2),
                        '"_gll_ny_m1_"': str(gll_ny - 1), 
                        '"_gll_nz_"': str(gll_nz),
                        '"grid_ll"': '"' + gllBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, replaceDic)
                    file_content = MultiWordReplace(file_content, gll_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                    
                    grid_b["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                    if "u" not in connections and "d" not in connections:
                        grid_b["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                    elif "b" not in connections:
                        grid_b["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down"])
                    else:
                        grid_b["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down_back"])
                    
                    for gll_updateSequence in grid_ll_conn["updateSequences"]:
                        seqName = gll_updateSequence["name"]
                        for updateSequence in grid_b["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                                if "u" not in connections and "d" not in connections:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                                elif "b" not in connections:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down"])
                                else:
                                    updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down_back"])
                
                ##----
                return grids

        
class GridCollectionner:
    def __init__(self):
        self.grids = []
        self.gridsByName = {}
        self.gridsByLevel = {}
        self.maxLevel = None
        self.numOfTimeSteps = None
        self.gridsArranged = False
        
    def SetNumOfCoarseTimeSteps(self, nt):
        self.numOfTimeSteps = nt
        
    def AddGrid(self, grid):
        self.grids.append(grid)
        self.gridsArranged = False
        
    def ArrangeGridsBasedOnLevel(self):
        if self.gridsArranged:
            return
        self.gridsByName = {}
        self.gridsByLevel = {}
        self.maxLevel = None
        for grid in self.grids:
            level = grid.blockLevel
            name = grid.name
            assert name not in self.gridsByName
            self.gridsByName[name] = grid
            
            if self.maxLevel is None:
                self.maxLevel = level
            elif self.maxLevel < level:
                self.maxLevel = level
            
            if level in self.gridsByLevel:
                self.gridsByLevel[level].append(name)            
            else:
                self.gridsByLevel[level] = [name]
        
        assert len(self.gridsByLevel) == self.maxLevel + 1
        
        sortDic = {"c":0, "r":1, "l":2, "u":3, "d":4, "f":5, "b":6}
        def sortFun(blockName):
            return sortDic[self.gridsByName[blockName].blockPosition]
        for i in range(self.maxLevel + 1):
            self.gridsByLevel[i].sort(key=sortFun, reverse=True)
        
        self.gridsArranged = True


    def GenerateGridCollection(self):
        self.ArrangeGridsBasedOnLevel()
        
        json_file = open('GridCollection.json')
        file_content = json_file.read()
        json_file.close()
        
        assert self.numOfTimeSteps != None
        replaceDic = {'"_nt_coarse_"': str(self.numOfTimeSteps)}
        file_content = MultiWordReplace(file_content, replaceDic)
        
        gridCollection = json.loads(file_content)
                
        for i in range(self.maxLevel + 1):
            gridNames_i = self.gridsByLevel[i]
            
            for gridName in gridNames_i:
                gridCollection["simulationParameters"]["grids"].update(self.gridsByName[gridName].SetupGrid())

        initSequence = gridCollection["simulationParameters"]["runSequence"][0]["sequence"]
        for i in range(self.maxLevel + 1):
            gridNames_i = self.gridsByLevel[i]
            
            for gridName in gridNames_i:
                initSequence.append([gridName, "reset_time"])

        timeSteppingSequence = gridCollection["simulationParameters"]["runSequence"][1]["sequence"]
        self.SetupEHSteps(self.maxLevel, timeSteppingSequence)
            
        return gridCollection

    def SetupEHSteps(self, level, timeSteppingSequence):
        if level < 0:
            return
        self.SetupHSteps(level, timeSteppingSequence)
        self.SetupEHSteps(level - 1, timeSteppingSequence)
        self.SetupESteps(level, timeSteppingSequence)
        self.SetupEHSteps(level - 1, timeSteppingSequence)
        
    def SetupHSteps(self, level, timeSteppingSequence):
        gridNames_lev = self.gridsByLevel[level]
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "update_H"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "update_H_out"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "update_E_edge_dt2"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "set_E_Edge_out"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "interpolate_E_edge_IP1_02"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "interpolate_E_edge_IP1_1"])


    def SetupESteps(self, level, timeSteppingSequence):
        gridNames_lev = self.gridsByLevel[level]
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "update_E_inside"])
            timeSteppingSequence.append([gridName, "update_E_edge_dt2"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "set_E_Edge_out"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "interpolate_E_edge_IP1_02"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "interpolate_E_edge_IP1_1"])
        for gridName in gridNames_lev:
            timeSteppingSequence.append([gridName, "update_time"])



        
