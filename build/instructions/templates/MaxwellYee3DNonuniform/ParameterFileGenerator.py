
__all__ = ["GridBlock", "GridCollectionner"]

import numpy as np
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
    def __init__(self, name, blockLevel, blockPosition, gridType = "normal"):
        """ blockLevel: 0, 1, 2
            blockPosition: "c" (center), "r" (right), "l" (left), "u" (up), "d" (down), "f" (front), "b" (back)
        """
        self.name = name
        self.blockLevel = blockLevel
        self.blockPosition = blockPosition
        self.connections = {}
        self.sources = []
        self.views = []
        self.materials = []
        self.geometries = []
        self.gridType = gridType    # "normal" or "pml"
        self.sig_e, self.sig_h = 1.0, 1.0
    
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

    def AddMaterial(self, materialParams):
        self.materials.append(materialParams)

    def AddGeometry(self, geomParams):
        self.geometries.append(geomParams)

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
                polarization = sourceParams["polarization"]
                amplitude = sourceParams["amplitude"]
                t_center = sourceParams["t_center"]
                t_decay = sourceParams["t_decay"]
                modulationFrequency = sourceParams["modulationFrequency"]
                modulationPhase = sourceParams["modulationPhase"]
                timeOffsetFraction = sourceParams["timeOffsetFraction"]
                dV = self.dx*self.dy*self.dz
                dA = None
                if polarization == 'x':
                    dA = self.dy*self.dz
                elif polarization == 'y':
                    dA = self.dx*self.dz
                elif polarization == 'z':
                    dA = self.dy*self.dx
                else:
                    assert False
                    
                
                replaceDic = {'"_indxJ_"': str(ind_x), '"_indyJ_"': str(ind_y), '"_indzJ_"': str(ind_z),
                              '"_indxJ_p1_"': str(ind_x + 1), '"_indyJ_p1_"': str(ind_y + 1), '"_indzJ_p1_"': str(ind_z + 1),
                              '"_j_polarization_"': '"' + polarization + '"', 
                              '"_j_amplitude_"': str(amplitude),
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

            if sourceParams["type"] == "GaussianLineSource_y":
                json_file = open('sources/GaussianInTime/J_line_y_3D.json')
                file_content = json_file.read()
                json_file.close()
                
                j_r = sourceParams["j_r"]
                height = sourceParams["height"]
                polarization = sourceParams["polarization"]
                amplitude = sourceParams["amplitude"]
                t_center = sourceParams["t_center"]
                t_decay = sourceParams["t_decay"]
                modulationFrequency = sourceParams["modulationFrequency"]
                modulationPhase = sourceParams["modulationPhase"]
                timeOffsetFraction = sourceParams["timeOffsetFraction"]
                dA = self.dx*self.dz
                
                dr = np.array([self.dx, self.dy, self.dz])
                j_inds = ((j_r - self.r0)/dr).astype("int")
                ny = int(height/self.dy)
                ind_y_0 = int(j_inds[1] - ny/2)
                ind_y_1 = ind_y_0 + ny
                
                if ind_y_0 < 0:
                    ind_y_0 = 0
                if ind_y_1 > self.ny:
                    ind_y_1 = self.ny
                if ny > self.ny:
                    ny = self.ny
                
                ind_x, ind_z = j_inds[0], j_inds[2]

                replaceDic = {'"_indxJ_"': str(ind_x), '"_indyJ_"': str(ind_y_0), '"_indzJ_"': str(ind_z),
                              '"_indxJ_p1_"': str(ind_x + 1), '"_indyJ_p_ny_"': str(ind_y_1), '"_indzJ_p1_"': str(ind_z + 1),
                              '"_ny_"': str(ny),
                              '"_j_polarization_"': '"' + polarization + '"', 
                              '"_j_amplitude_"': str(amplitude),
                              '"_j_t_center_"': str(t_center), '"_j_t_decay_"': str(t_decay),
                              '"_j_mod_freq_"': str(modulationFrequency), '"_j_mod_phase_"': str(modulationPhase),
                              '"_j_time_offset_frac_"': str(timeOffsetFraction),
                              '"_m_dt_dA_"': str(-self.dt/dA),
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

            if sourceParams["type"] == "GaussianSheetSource_z":
                json_file = open('sources/GaussianInTime/J_sheet_z_3D.json')
                file_content = json_file.read()
                json_file.close()
                
                j_r = sourceParams["j_r"]
                y_width = sourceParams["y_width"]
                x_width = sourceParams["x_width"]
                polarization = sourceParams["polarization"]
                amplitude = sourceParams["amplitude"]
                t_center = sourceParams["t_center"]
                t_decay = sourceParams["t_decay"]
                modulationFrequency = sourceParams["modulationFrequency"]
                modulationPhase = sourceParams["modulationPhase"]
                timeOffsetFraction = sourceParams["timeOffsetFraction"]
                dw = self.dz    ## thickness
                
                dr = np.array([self.dx, self.dy, self.dz])
                j_inds = ((j_r - self.r0)/dr).astype("int")
                nx = int(x_width/self.dx)
                ny = int(y_width/self.dy)
                ind_x_0 = int(j_inds[0] - nx/2 - 1)
                ind_x_1 = ind_x_0 + nx + 1
                ind_y_0 = int(j_inds[1] - ny/2)
                ind_y_1 = ind_y_0 + ny
                
                if ind_x_0 < 0:
                    ind_x_0 = 0
                if ind_x_1 > self.nx + 1:
                    ind_x_1 = self.nx + 1
                if nx > self.nx:
                    nx = self.nx

                if ind_y_0 < 0:
                    ind_y_0 = 0
                if ind_y_1 > self.ny:
                    ind_y_1 = self.ny
                if ny > self.ny:
                    ny = self.ny
                
                ind_x, ind_z = j_inds[0], j_inds[2]

                replaceDic = {'"_indxJ_"': str(ind_x_0), '"_indyJ_"': str(ind_y_0), '"_indzJ_"': str(ind_z),
                              '"_indxJ_p_nx_p1_"': str(ind_x_1), '"_indyJ_p_ny_"': str(ind_y_1), '"_indzJ_p1_"': str(ind_z + 1),
                              '"_nx_"': str(nx), '"_ny_"': str(ny),
                              '"_j_polarization_"': '"' + polarization + '"', 
                              '"_j_amplitude_"': str(amplitude),
                              '"_j_t_center_"': str(t_center), '"_j_t_decay_"': str(t_decay),
                              '"_j_mod_freq_"': str(modulationFrequency), '"_j_mod_phase_"': str(modulationPhase),
                              '"_j_time_offset_frac_"': str(timeOffsetFraction),
                              '"_m_dt_dw_"': str(-self.dt/dw),
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


            elif sourceParams["type"] == "PureScatteredRectPlaneWave":
                json_file = open('sources/PureScatteredIncidentPlaneWave/RectPlaneWave.json')
                file_content = json_file.read()
                json_file.close()

                ## get partial grid start and end indices         
                r0, r1 = sourceParams["boundingBox"]    ##  bounding box around the metal
                dr = np.array([self.dx, self.dy, self.dz])
                inds_0 = np.round((r0 - self.r0) / dr).astype(int) - 1
                inds_1 = np.round((r1 - self.r0) / dr).astype(int) + 1
                self.truncateIndsToRange(inds_0)
                self.truncateIndsToRange(inds_1)
                if np.prod(inds_1 - inds_0) == 0:
                    continue
                inds_0 -= (inds_0 % 2)      # align to even grids
                    
                polarization = sourceParams["polarization"]
                amplitude = sourceParams["amplitude"]
                propagationDirection = sourceParams["propagationDirection"]
                propagationVelocity = sourceParams["propagationVelocity"]
                t_center = sourceParams["t_center"]
                rectWidth = sourceParams["rectWidth"]
                rectEdgeWidth = sourceParams["rectEdgeWidth"]
                modulationFrequency = sourceParams["modulationFrequency"]
                modulationPhase = sourceParams["modulationPhase"]
                
                if polarization != 'y':
                    assert False
                                
                n = (inds_1 - inds_0).astype(int)
                replaceDic = {
                    '"_indx0_"': str(inds_0[0]), '"_indy0_"': str(inds_0[1]), '"_indz0_"': str(inds_0[2]), 
                    '"_nx_"': str(n[0]), '"_ny_"': str(n[1]), '"_nz_"': str(n[2]),
                    '"_polarization_"': '"' + polarization + '"',
                    '"_planewave_direction_x_"': str(propagationDirection[0]),
                    '"_planewave_direction_y_"': str(propagationDirection[1]),
                    '"_planewave_direction_z_"': str(propagationDirection[2]),
                    '"_planewave_velocity_"': str(propagationVelocity),
                    '"_Einc_amp_"': str(amplitude),
                    '"_planewave_t_center_"': str(t_center),
                    '"_rect_width_"': str(rectWidth),
                    '"_rect_edge_width_"': str(rectEdgeWidth),
                    '"_planewave_freq_"': str(modulationFrequency),
                    '"_planewave_phase_"': str(modulationPhase),
                    '"_m_dt_"': str(-self.dt),
                    "__ext__": "_" + str(sourceParams["geometryName"])
                    }

                file_content = MultiWordReplace(file_content, replaceDic)
                sourceData = json.loads(file_content)
                
                grid["partialGridArrays"].extend(sourceData["partialGridArrays"])
                grid["girdArrayManipulators"].extend(sourceData["girdArrayManipulators"])
                grid["updateInstructions"].extend(sourceData["updateInstructions"])
                for s_updateSequence in sourceData["updateSequences"]:
                    seqName = s_updateSequence["name"]
                    for updateSequence in grid["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(s_updateSequence["sequence"])
                
                
    def truncateIndsToRange(self, inds):
        #print("inds: ", inds)
        n = np.array([self.nx, self.ny, self.nz])
        for i in range(3):
            if inds[i] < 0:
                inds[i] = 0
            elif inds[i] > n[i]:
                inds[i] = n[i]
        #print("truncated inds: ", inds)
    
    
    def SetupMaterials(self, grid):
        for materialIndex in range(len(self.materials)):
            materialParams = self.materials[materialIndex]
            if materialParams["type"] == "DrudeMetal_PureScattered":
                file_content = None
                if "wireMeshAlong" in materialParams and materialParams["wireMeshAlong"] == "y":
                    json_file = open('materials/metal/Drude_y_pureScattered.json')
                    file_content = json_file.read()
                    json_file.close()                
                else:
                    json_file = open('materials/metal/Drude_pureScattered.json')
                    file_content = json_file.read()
                    json_file.close()
    
                ## get partial grid start and end indices         
                r0, r1 = materialParams["boundingBox"]    ##  bounding box around the metal
                dr = np.array([self.dx, self.dy, self.dz])
                inds_0 = np.round((r0 - self.r0) / dr).astype(int)
                inds_1 = np.round((r1 - self.r0) / dr).astype(int)
                self.truncateIndsToRange(inds_0)
                self.truncateIndsToRange(inds_1)
                if np.prod(inds_1 - inds_0) == 0:
                    continue
                inds_0 -= (inds_0 % 2)      # align to even grids
                    
                geometryName = materialParams["geometryName"]
                wp = materialParams["plasmaFrequency"]
                gamma = materialParams["scatteringRate"]
                wp_sq = wp**2
                
                if "wireMeshAlong" in materialParams and materialParams["wireMeshAlong"] == "y":
                    geom = None
                    for geomParams in self.geometries:
                        if geomParams["geometryName"] == geometryName:
                            geom = geomParams
                            break
                    assert geom != None
                    if geom["type"] == "cylinder": 
                        #print(geom)
                        if  geom["alignEven"] == "yes":
                            wp_sq *= 4
                            print("adjusting cylinder amp")
                
                                
                n = (inds_1 - inds_0).astype(int) 
                replaceDic = {
                    '"_indx0_"': str(inds_0[0]), '"_indy0_"': str(inds_0[1]), '"_indz0_"': str(inds_0[2]), 
                    '"_indx1_"': str(inds_1[0]), '"_indy1_"': str(inds_1[1]), '"_indz1_"': str(inds_1[2]), 
                    '"_indx0_p1_"': str(inds_0[0] + 1), '"_indy0_p1_"': str(inds_0[1] + 1), '"_indz0_p1_"': str(inds_0[2] + 1), 
                    '"_indx1_p1_"': str(inds_1[0] + 1), '"_indy1_p1_"': str(inds_1[1] + 1), '"_indz1_p1_"': str(inds_1[2] + 1), 
                    '"_nx_"': str(n[0]), '"_ny_"': str(n[1]), '"_nz_"': str(n[2]),
                    '"_nx_p1_"': str(n[0] + 1), '"_ny_p1_"': str(n[1] + 1), '"_nz_p1_"': str(n[2] + 1),
                    '"_geometry_name_"': '"' + geometryName + '"',
                    '"_wp_sq_"': str(wp_sq),
                    '"_gamma_"': str(gamma),
                    '"_m_dt_"': str(-self.dt),
                    '"_dt_"': str(self.dt),
                    '"_dt_2_"': str(self.dt/2.0),
                    '"_m_dt_2_"': str(-self.dt/2.0),
                    "__ext__": "_" + geometryName
                    }
    
                file_content = MultiWordReplace(file_content, replaceDic)
                matData = json.loads(file_content)
                
                grid["partialGridArrays"].extend(matData["partialGridArrays"])
                grid["girdArrayManipulators"].extend(matData["girdArrayManipulators"])
                grid["updateInstructions"].extend(matData["updateInstructions"])
                for m_updateSequence in matData["updateSequences"]:
                    seqName = m_updateSequence["name"]
                    for updateSequence in grid["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(m_updateSequence["sequence"])
            elif materialParams["type"] == "pec_PureScattered":
                json_file = open('materials/metal/pec_pureScattered.json')
                file_content = json_file.read()
                json_file.close()
    
                ## get partial grid start and end indices         
                r0, r1 = materialParams["boundingBox"]    ##  bounding box around the metal
                dr = np.array([self.dx, self.dy, self.dz])
                inds_0 = np.round((r0 - self.r0) / dr).astype(int) - 1
                inds_1 = np.round((r1 - self.r0) / dr).astype(int) + 1
                self.truncateIndsToRange(inds_0)
                self.truncateIndsToRange(inds_1)
                if np.prod(inds_1 - inds_0) == 0:
                    continue
                inds_0 -= (inds_0 % 2)      # align to even grids
                    
                geometryName = materialParams["geometryName"]                                
                n = (inds_1 - inds_0).astype(int) 
                replaceDic = {
                    '"_indx0_"': str(inds_0[0]), '"_indy0_"': str(inds_0[1]), '"_indz0_"': str(inds_0[2]), 
                    '"_indx1_"': str(inds_1[0]), '"_indy1_"': str(inds_1[1]), '"_indz1_"': str(inds_1[2]), 
                    '"_indx0_p1_"': str(inds_0[0] + 1), '"_indy0_p1_"': str(inds_0[1] + 1), '"_indz0_p1_"': str(inds_0[2] + 1), 
                    '"_indx1_p1_"': str(inds_1[0] + 1), '"_indy1_p1_"': str(inds_1[1] + 1), '"_indz1_p1_"': str(inds_1[2] + 1), 
                    '"_nx_"': str(n[0]), '"_ny_"': str(n[1]), '"_nz_"': str(n[2]),
                    '"_nx_p1_"': str(n[0] + 1), '"_ny_p1_"': str(n[1] + 1), '"_nz_p1_"': str(n[2] + 1),
                    '"_geometry_name_"': '"' + geometryName + '"',
                    "__ext__": "_" + str(geometryName)
                    }
    
                file_content = MultiWordReplace(file_content, replaceDic)
                matData = json.loads(file_content)
                
                """
                ## merge partial grids
                for partialGAList in matData["partialGridArrays"]:
                    partialGAname = partialGAList["name"]
                    exists = False
                    for ind in range(len(grid["partialGridArrays"])):
                        g_partialGA = grid["partialGridArrays"][ind]
                        if g_partialGA["name"] == partialGAname:
                            exists = True
                            # merge
                            partialGA_merged = self.MergeParialGridArrays(g_partialGA, partialGAList)
                            grid["partialGridArrays"][ind] = partialGA_merged
                            break
                    if not exists:
                        grid["partialGridArrays"].append(partialGAList)
                """
                
                grid["partialGridArrays"].extend(matData["partialGridArrays"])        
                grid["girdArrayManipulators"].extend(matData["girdArrayManipulators"])
                grid["updateInstructions"].extend(matData["updateInstructions"])
                for m_updateSequence in matData["updateSequences"]:
                    seqName = m_updateSequence["name"]
                    for updateSequence in grid["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(m_updateSequence["sequence"])

            elif materialParams["type"] == "pec":
                json_file = open('materials/metal/pec.json')
                file_content = json_file.read()
                json_file.close()
    
                ## get partial grid start and end indices         
                r0, r1 = materialParams["boundingBox"]    ##  bounding box around the metal
                dr = np.array([self.dx, self.dy, self.dz])
                inds_0 = np.round((r0 - self.r0) / dr).astype(int) - 1
                inds_1 = np.round((r1 - self.r0) / dr).astype(int) + 1
                self.truncateIndsToRange(inds_0)
                self.truncateIndsToRange(inds_1)
                if np.prod(inds_1 - inds_0) == 0:
                    continue
                inds_0 -= (inds_0 % 2)      # align to even grids
                    
                geometryName = materialParams["geometryName"]                                
                n = (inds_1 - inds_0).astype(int) 
                replaceDic = {
                    '"_indx0_"': str(inds_0[0]), '"_indy0_"': str(inds_0[1]), '"_indz0_"': str(inds_0[2]), 
                    '"_indx1_"': str(inds_1[0]), '"_indy1_"': str(inds_1[1]), '"_indz1_"': str(inds_1[2]), 
                    '"_indx0_p1_"': str(inds_0[0] + 1), '"_indy0_p1_"': str(inds_0[1] + 1), '"_indz0_p1_"': str(inds_0[2] + 1), 
                    '"_indx1_p1_"': str(inds_1[0] + 1), '"_indy1_p1_"': str(inds_1[1] + 1), '"_indz1_p1_"': str(inds_1[2] + 1), 
                    '"_nx_"': str(n[0]), '"_ny_"': str(n[1]), '"_nz_"': str(n[2]),
                    '"_nx_p1_"': str(n[0] + 1), '"_ny_p1_"': str(n[1] + 1), '"_nz_p1_"': str(n[2] + 1),
                    '"_geometry_name_"': '"' + geometryName + '"',
                    "__ext__": "_" + str(geometryName)
                    }
    
                file_content = MultiWordReplace(file_content, replaceDic)
                matData = json.loads(file_content)
                
                grid["partialGridArrays"].extend(matData["partialGridArrays"])        
                grid["girdArrayManipulators"].extend(matData["girdArrayManipulators"])
                grid["updateInstructions"].extend(matData["updateInstructions"])
                for m_updateSequence in matData["updateSequences"]:
                    seqName = m_updateSequence["name"]
                    for updateSequence in grid["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(m_updateSequence["sequence"])
            else:
                assert False
 
    def MergeParialGridArrays(self, partialGA1, partialGA2):
        assert partialGA1["name"] == partialGA2["name"]
        assert partialGA1["type"] == partialGA2["type"]
        indStart = np.minimum(np.array(partialGA1["indStart"]), np.array(partialGA2["indStart"]))
        indEnd = np.maximum(np.array(partialGA1["indStart"]) + np.array(partialGA1["nCells"]), 
                            np.array(partialGA2["indStart"]) + np.array(partialGA2["nCells"]))
        nCells = indEnd - indStart
        
        partialGA_merged = {"name":partialGA1["name"] , "type":partialGA1["type"], 
                            "indStart":indStart.tolist(), "nCells":nCells.tolist()}
       
        return partialGA_merged
        
     
    def SetupGeometries(self, grid):
        for geomIndex in range(len(self.geometries)):
            geomParams = self.geometries[geomIndex]
            if geomParams["type"] == "hyperboloid":
                json_file = open('geometries/hyperboloid.json')
                file_content = json_file.read()
                json_file.close()
    
                ## get partial grid start and end indices         
                r0, r1 = geomParams["boundingBox"]    ##  bounding box around the metal
                dr = np.array([self.dx, self.dy, self.dz])
                inds_0 = np.round((r0 - self.r0) / dr).astype(int)
                inds_1 = np.round((r1 - self.r0) / dr).astype(int)
                self.truncateIndsToRange(inds_0)
                self.truncateIndsToRange(inds_1)
                if np.prod(inds_1 - inds_0) == 0:
                    continue
                inds_0 -= (inds_0 % 2)      # align to even grids
                    
                geometryName = geomParams["geometryName"]
                coneAngle = geomParams["coneAngle"]
                tipRadius = geomParams["tipRadius"]
                height = geomParams["height"]
                apexPosition = geomParams["apexPosition"]

                replaceDic = {
                    '"_name_"': '"' + geometryName + '"',
                    '"_coneAngle_"': str(coneAngle),
                    '"_tipRadius_"': str(tipRadius),
                    '"_height_"': str(height),
                    '"_apex_x_"': str(apexPosition[0]),
                    '"_apex_y_"': str(apexPosition[1]),
                    '"_apex_z_"': str(apexPosition[2])
                    }
    
                file_content = MultiWordReplace(file_content, replaceDic)
                geomData = json.loads(file_content)
                
                #print(geomData)
           
                if "geometries" in grid:
                    grid["geometries"].extend(geomData["geometries"])
                else:
                    grid["geometries"] = geomData["geometries"]
                    
            elif geomParams["type"] == "cone":
                json_file = open('geometries/cone.json')
                file_content = json_file.read()
                json_file.close()
    
                ## get partial grid start and end indices         
                r0, r1 = geomParams["boundingBox"]    ##  bounding box around the metal
                dr = np.array([self.dx, self.dy, self.dz])
                inds_0 = np.round((r0 - self.r0) / dr).astype(int)
                inds_1 = np.round((r1 - self.r0) / dr).astype(int)
                self.truncateIndsToRange(inds_0)
                self.truncateIndsToRange(inds_1)
                if np.prod(inds_1 - inds_0) == 0:
                    continue
                inds_0 -= (inds_0 % 2)      # align to even grids
                    
                geometryName = geomParams["geometryName"]
                coneAngle = geomParams["coneAngle"]
                tipRadius = geomParams["tipRadius"]
                height = geomParams["height"]
                apexPosition = geomParams["apexPosition"]

                replaceDic = {
                    '"_name_"': '"' + geometryName + '"',
                    '"_coneAngle_"': str(coneAngle),
                    '"_tipRadius_"': str(tipRadius),
                    '"_height_"': str(height),
                    '"_apex_x_"': str(apexPosition[0]),
                    '"_apex_y_"': str(apexPosition[1]),
                    '"_apex_z_"': str(apexPosition[2])
                    }
    
                file_content = MultiWordReplace(file_content, replaceDic)
                geomData = json.loads(file_content)
                
                #print(geomData)
           
                if "geometries" in grid:
                    grid["geometries"].extend(geomData["geometries"])
                else:
                    grid["geometries"] = geomData["geometries"]
            
            elif geomParams["type"] == "cylinder":
                json_file = open('geometries/cylinder.json')
                file_content = json_file.read()
                json_file.close()
    
                ## get partial grid start and end indices         
                r0, r1 = geomParams["boundingBox"]    ##  bounding box around the metal
                dr = np.array([self.dx, self.dy, self.dz])
                inds_0 = np.round((r0 - self.r0) / dr).astype(int)
                inds_1 = np.round((r1 - self.r0) / dr).astype(int)
                self.truncateIndsToRange(inds_0)
                self.truncateIndsToRange(inds_1)
                if np.prod(inds_1 - inds_0) == 0:
                    continue
                inds_0 -= (inds_0 % 2)      # align to even grids
                    
                geometryName = geomParams["geometryName"]
                radius = geomParams["radius"]
                height = geomParams["height"]
                topCenter = geomParams["topCenter"]
                alignEven = geomParams["alignEven"]     ## yes/no

                replaceDic = {
                    '"_name_"': '"' + geometryName + '"',
                    '"_radius_"': str(radius),
                    '"_height_"': str(height),
                    '"_top_x_"': str(topCenter[0]),
                    '"_top_y_"': str(topCenter[1]),
                    '"_top_z_"': str(topCenter[2]),
                    '"_align_even_"': '"' + alignEven + '"'
                    }
    
                file_content = MultiWordReplace(file_content, replaceDic)
                geomData = json.loads(file_content) 
                
                #print(geomData)
           
                if "geometries" in grid:
                    grid["geometries"].extend(geomData["geometries"])
                else:
                    grid["geometries"] = geomData["geometries"]
            else:
                raise NotImplementedError

        
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
                    bufferSize = 1024*1024*10;
                    if "bufferSize"in viewParams:
                        bufferSize = viewParams["bufferSize"]
                    
                    replaceDic = {'"_indxSave_0_"': str(indx0), '"_indySave_0_"': str(indy0), '"_indzSave_0_"': str(indz0),
                                  '"_indxSave_1_"': str(indx1), '"_indySave_1_"': str(indy1), '"_indzSave_1_"': str(indz1),
                                  '"_save_rate_"': str(saveRate),
                                  '"_buffer_size_"': str(bufferSize),
                                  '"_array_"': '"' + arrayName + '"', 
                                  '"_direction_"': '"' + direction + '"', 
                                  '"_file_name_"': '"' + fileName + '"'
                                  }
                    file_content = MultiWordReplace(file_content, replaceDic)
                    viewData = json.loads(file_content)
                    
                    grid["gridViews"].extend(viewData["gridViews"])

            elif viewParams["type"] == "entire":
                json_file = open('views/entire.json')
                file_content = json_file.read()
                json_file.close()
            
                saveRate = viewParams['saveRate']
                arrayName = viewParams['arrayName']
                fileName = viewParams['fileName']
                direction = viewParams["direction"]     ## "x", "y", "z"
            
                replaceDic = {'"_save_rate_"': str(saveRate),
                              '"_array_"': '"' + arrayName + '"', 
                              '"_direction_"': '"' + direction + '"', 
                              '"_file_name_"': '"' + fileName + '"'
                              }
                file_content = MultiWordReplace(file_content, replaceDic)
                viewData = json.loads(file_content)
                
                grid["gridViews"].extend(viewData["gridViews"])

    def SetupSourcesMateriasAndViews(self, grid):
        self.SetupSources(grid)
        self.SetupGeometries(grid)
        self.SetupMaterials(grid)
        self.SetupViews(grid)
    

        
    
    def SetupGrid(self):
        if self.gridType == "normal":
            return self.SetupNormalGrid()
        elif self.gridType == "pml":
            return self.SetupPMLGrid()
        
        
    def SetupNormalGrid(self):
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
            
            self.SetupSourcesMateriasAndViews(grid_m)
            
            connections = self.connections
            
            #------------------- r ------------------
            #----------------------------------------
            if "r" in connections:
                self.ConnectToRight(grid_m)
            
            #------------------- l ------------------
            #----------------------------------------
            if "l" in connections:
                self.ConnectToLeft(grid_m)

            #------------------- u ------------------
            #----------------------------------------
            if "u" in connections:
                self.ConnectToUp(grid_m)
                
            #------------------- d ------------------
            #----------------------------------------
            if "d" in connections:
                self.ConnectToDown(grid_m)
                
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

            ##---- grid_m
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
                
                self.SetupSourcesMateriasAndViews(grid_r)
                    
                connections = self.connections
                
                #------------------- u ------------------
                #----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_r)
                    
                #------------------- d ------------------
                #----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_r)
                    
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
                #-----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_r)

                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_r)

                ##---- grid_r
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
                
                self.SetupSourcesMateriasAndViews(grid_rr)
                    
                connections = self.connections

                #------------------- uu ------------------
                #-----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_rr)
                
                #------------------- dd ------------------
                #-----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_rr)
                    
                #------------------- ff ------------------
                #-----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_r/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _nz = connections['lf'].nz
                    _gr_nx, _gr_ny, _gr_nz = self.nx, self.ny, self.nz
                    _gl_nz = self.nz
                    _gd_ny = connections["d"].ny
                    gff_replaceDic = {
                        '"_gr_nx_"': str(_gr_nx),
                        '"_gr_nx_p1_"': str(_gr_nx + 1),
                        '"_gr_ny_"': str(_gr_ny),
                        '"_gr_ny_p1_"': str(_gr_ny + 1),
                        '"_gr_nz_"': str(_gr_nz),
                        '"_gd_ny_"': str(_gd_ny),
                        '"_gl_nz_p_nz2_"': str(_gl_nz + int(_nz/2)),
                        '"grid_f"': '"' + connections["f"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gff_replaceDic)
                    grid_ff_conn = json.loads(file_content)
                    
                    grid_rr["updateInstructions"].extend(grid_ff_conn["updateInstructions"])
                    
                    for gff_updateSequence in grid_ff_conn["updateSequences"]:
                        seqName = gff_updateSequence["name"]
                        for updateSequence in grid_rr["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gff_updateSequence["sequence"])


                #------------------- bb ------------------
                #-----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_r/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _gr_ny, _gr_nz = self.ny, self.nz
                    _gb_nx = connections["b"].nx
                    _gl_nz = self.nz
                    _gd_ny = connections["d"].ny
                    _nz = connections["lb"].nz
                    gbb_replaceDic = {
                        '"_gr_ny_"': str(_gr_ny), '"_gr_ny_p1_"': str(_gr_ny + 1),
                        '"_gr_nz_"': str(_gr_nz),
                        '"_gd_ny_"': str(_gd_ny),
                        '"_gb_nx_"': str(_gb_nx),
                        '"_gl_nz_p_nz2_"': str(_gl_nz + int(_nz/2)),
                        '"grid_b"': '"' + connections["b"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gbb_replaceDic)
                    grid_bb_conn = json.loads(file_content)
                    
                    grid_rr["updateInstructions"].extend(grid_bb_conn["updateInstructions"])
                    
                    for gbb_updateSequence in grid_bb_conn["updateSequences"]:
                        seqName = gbb_updateSequence["name"]
                        for updateSequence in grid_rr["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gbb_updateSequence["sequence"])

                #------------------- rrr ------------------
                #-----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_rr)

                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_rr)


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
                
                self.SetupSourcesMateriasAndViews(grid_l)

                connections = self.connections
                
                #------------------- u ------------------
                #----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_l)
                
                #------------------- d ------------------
                #----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_l)
                                    
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
                #-----------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_l)
                    
                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_l)

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
                
                self.SetupSourcesMateriasAndViews(grid_ll)

                connections = self.connections

                #------------------- uu ------------------
                #-----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_ll)
                
                #------------------- dd ------------------
                #-----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_ll)
                    
                #------------------- ff ------------------
                #-----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_l/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _gl_nx, _gl_ny, _gl_nz = self.nx, self.ny, self.nz
                    _gd_ny = connections["d"].ny
                    gff_replaceDic = {
                        '"_gl_nx_"': str(_gl_nx), '"_gl_nx_p1_"': str(_gl_nx + 1),
                        '"_gl_ny_"': str(_gl_ny), '"_gl_ny_p1_"': str(_gl_ny + 1),
                        '"_gl_nz_"': str(_gl_nz), '"_gl_nz_p1_"': str(_gl_nz + 1),
                        '"_gd_ny_"': str(_gd_ny),
                        '"grid_f"': '"' + connections["f"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gff_replaceDic)
                    grid_ff_conn = json.loads(file_content)
                    
                    grid_ll["updateInstructions"].extend(grid_ff_conn["updateInstructions"])
                    
                    for gff_updateSequence in grid_ff_conn["updateSequences"]:
                        seqName = gff_updateSequence["name"]
                        for updateSequence in grid_ll["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gff_updateSequence["sequence"])

                #------------------- bb ------------------
                #-----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_l/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _gl_nx, _gl_ny, _gl_nz = self.nx, self.ny, self.nz
                    _gb_nx = connections["b"].nx
                    _gd_ny = connections["d"].ny
                    gbb_replaceDic = {
                        '"_gl_ny_"': str(_gl_ny), '"_gl_ny_p1_"': str(_gl_ny + 1),
                        '"_gl_nz_"': str(_gl_nz), '"_gl_nz_p1_"': str(_gl_nz + 1),
                        '"_gd_ny_"': str(_gd_ny),
                        '"_gb_nx_"': str(_gb_nx),
                        '"grid_b"': '"' + connections["b"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gbb_replaceDic)
                    grid_bb_conn = json.loads(file_content)
                    
                    grid_ll["updateInstructions"].extend(grid_bb_conn["updateInstructions"])
                    
                    for gbb_updateSequence in grid_bb_conn["updateSequences"]:
                        seqName = gbb_updateSequence["name"]
                        for updateSequence in grid_ll["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gbb_updateSequence["sequence"])

                #------------------- lll ------------------
                #-----------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_ll)

                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_ll)

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
                
                self.SetupSourcesMateriasAndViews(grid_u)
                
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
                #-----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_u)

                #------------------- ll ------------------
                #-----------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_u)
                                
                #------------------- uu -----------------
                #----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_u)
                
                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_u)

                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_u)

                ##------ grid_u        
                return grids

            ##---------------------------------------------- grid_uu ------------------------------------------------
            ##-------------------------------------------------------------------------------------------------------
            elif self.blockLevel >= 2:
                json_file = open('layer2/grid_uu/grid_uu.json')
                file_content = json_file.read()
                json_file.close()
                
                guu_nx, guu_ny, guu_nz = self.nx, self.ny, self.nz
                guu_r0, guu_r1 = self.r0, self.r1
                guu_dt, guu_dx, guu_dy, guu_dz = self.dt, self.dx, self.dy, self.dz 
                assert "dc" in self.connections
                assert "dl" in self.connections
                assert "dr" in self.connections
                assert "df" in self.connections
                assert "db" in self.connections
                assert self.connections["dc"].blockLevel == self.blockLevel - 1
                assert self.connections["df"].blockLevel == self.blockLevel - 1
                assert self.connections["db"].blockLevel == self.blockLevel - 1
                assert self.connections["dl"].blockLevel == self.blockLevel
                assert self.connections["dr"].blockLevel == self.blockLevel
        
                gll_ny, gll_nz = self.connections["dl"].ny, self.connections["dl"].nz
                grr_ny = self.connections["dr"].ny
                gf_ny, gf_nz = self.connections["df"].ny, self.connections["df"].nz 
                gb_nx, gb_ny = self.connections["db"].nx, self.connections["db"].ny
                gu_nx, gu_ny = self.connections["dc"].nx, self.connections["dc"].ny
                replaceDic = {'"_guu_nx_"': str(guu_nx), '"_guu_ny_"': str(guu_ny), '"_guu_nz_"': str(guu_nz),
                    '"_guu_nx_p1_"': str(guu_nx + 1), '"_guu_ny_p1_"': str(guu_ny + 1), '"_guu_nz_p1_"': str(guu_nz + 1),
                    '"_guu_nz_m1_"': str(guu_nz - 1),
                    '"_guu_x0_"': str(guu_r0[0]), '"_guu_y0_"': str(guu_r0[1]), '"_guu_z0_"': str(guu_r0[2]),
                    '"_guu_x1_"': str(guu_r1[0]), '"_guu_y1_"': str(guu_r1[1]), '"_guu_z1_"': str(guu_r1[2]),
                    '"_guu_dt_dx_"': str(guu_dt/guu_dx), '"_guu_dt_dy_"': str(guu_dt/guu_dy), '"_guu_dt_dz_"': str(guu_dt/guu_dz),
                    '"_guu_m_dt_dx_"': str(-guu_dt/guu_dx), '"_guu_m_dt_dy_"': str(-guu_dt/guu_dy), '"_guu_m_dt_dz_"': str(-guu_dt/guu_dz),
                    '"_guu_dt_dx_2_"': str(guu_dt/guu_dx/2.0), '"_guu_dt_dy_2_"': str(guu_dt/guu_dy/2.0), '"_guu_dt_dz_2_"': str(guu_dt/guu_dz/2.0),
                    '"_guu_m_dt_dx_2_"': str(-guu_dt/guu_dx/2.0), '"_guu_m_dt_dy_2_"': str(-guu_dt/guu_dy/2.0), 
                    '"_guu_m_dt_dz_2_"': str(-guu_dt/guu_dz/2.0),
                    '"_guu_dt_"': str(guu_dt),
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_p1_"': str(gll_nz + 1),
                    '"_gll_ny_m1_"': str(gll_ny - 1),
                    '"_grr_ny_m1_"': str(grr_ny - 1),
                    '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                    '"_gb_nx2_"': str(int(gb_nx/2)),
                    '"_gb_ny_m1_"': str(gb_ny - 1), '"_gb_ny_m2_"': str(gb_ny - 2),
                    '"_gb_nx2_p_gu_nx2_"': str(int(gb_nx/2) + int(gu_nx/2)),
                    '"_gu_ny_m1_"': str(gu_ny - 1), '"_gu_ny_m2_"': str(gu_ny - 2),
                    '"_gf_ny_m1_"': str(gf_ny - 1), '"_gf_ny_m2_"': str(gf_ny - 2),
                    '"grid_uu"': '"' + self.name + '"',
                    '"grid_u"':  '"' + self.connections["dc"].name + '"',
                    '"grid_rr"': '"' + self.connections["dr"].name + '"',
                    '"grid_ll"': '"' + self.connections["dl"].name + '"',
                    '"grid_f"':  '"' + self.connections["df"].name + '"',
                    '"grid_b"':  '"' + self.connections["db"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_uu = grids[self.name]
                
                self.SetupSourcesMateriasAndViews(grid_uu)
                
                connections = self.connections

                #------------------- ff ------------------
                #-----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_u/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _gu_nx, _gu_ny, _gu_nz = self.nx, self.ny, self.nz
                    _gd_ny = self.ny
                    _ny = connections["df"].ny      ## f: to take the entire lower level middle cube's ny
                    gff_replaceDic = {
                        '"_gu_nx_"': str(_gu_nx), '"_gu_nx_p1_"': str(_gu_nx + 1), 
                        '"_gu_ny_"': str(_gu_ny), '"_gu_ny_p1_"': str(_gu_ny + 1), 
                        '"_gu_nz_"': str(_gu_nz), '"_gu_nz_p1_"': str(_gu_nz + 1), 
                        '"_gd_ny_p_ny2_"': str(_gd_ny + int(_ny/2)),
                        '"grid_f"': '"' + connections["f"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gff_replaceDic)
                    grid_ff_conn = json.loads(file_content)
                    
                    grid_uu["updateInstructions"].extend(grid_ff_conn["updateInstructions"])
                    
                    for gff_updateSequence in grid_ff_conn["updateSequences"]:
                        seqName = gff_updateSequence["name"]
                        for updateSequence in grid_uu["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gff_updateSequence["sequence"])

                #------------------- bb ------------------
                #-----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_u/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _gu_nx, _gu_ny, _gu_nz = self.nx, self.ny, self.nz
                    _gb_nx = connections["b"].nx
                    _gd_ny = self.ny
                    _ny = connections["df"].ny      ## f: to take the entire lower level middle cube's ny
                    gbb_replaceDic = {
                        '"_gu_ny_"': str(_gu_ny), '"_gu_ny_p1_"': str(_gu_ny + 1), 
                        '"_gu_nz_"': str(_gu_nz), '"_gu_nz_p1_"': str(_gu_nz + 1), 
                        '"_gb_nx_"': str(_gb_nx),
                        '"_gd_ny_p_ny2_"': str(_gd_ny + int(_ny/2)),
                        '"grid_b"': '"' + connections["b"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gbb_replaceDic)
                    grid_bb_conn = json.loads(file_content)
                    
                    grid_uu["updateInstructions"].extend(grid_bb_conn["updateInstructions"])
                    
                    for gbb_updateSequence in grid_bb_conn["updateSequences"]:
                        seqName = gbb_updateSequence["name"]
                        for updateSequence in grid_uu["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gbb_updateSequence["sequence"])

                #------------------- rrr ------------------
                #-----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_uu)
                    
                #------------------- lll ------------------
                #-----------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_uu)
                    
                #------------------- uuu -----------------
                #-----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_uu)

                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_uu)

                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_uu)

                ##------ grid_uu        
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
                
                self.SetupSourcesMateriasAndViews(grid_d)

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
        
                #------------------- rr -----------------
                #----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_d)
                                
                #------------------- ll ------------------
                #----------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_d)
                                
                #------------------- dd ------------------
                #-----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_d)

                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_d)
                
                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_d)
                
                ##------ grid_d
                return grids

            ##---------------------------------------------- grid_dd ------------------------------------------------
            ##-------------------------------------------------------------------------------------------------------
            if self.blockLevel >= 2:
                json_file = open('layer2/grid_dd/grid_dd.json')
                file_content = json_file.read()
                json_file.close()
                
                gdd_nx, gdd_ny, gdd_nz = self.nx, self.ny, self.nz
                gdd_r0, gdd_r1 = self.r0, self.r1
                gdd_dt, gdd_dx, gdd_dy, gdd_dz = self.dt, self.dx, self.dy, self.dz 
                assert "uc" in self.connections
                assert "ul" in self.connections
                assert "ur" in self.connections
                assert "uf" in self.connections
                assert "ub" in self.connections
                assert self.connections["uc"].blockLevel == self.blockLevel - 1
                assert self.connections["uf"].blockLevel == self.blockLevel - 1
                assert self.connections["ub"].blockLevel == self.blockLevel - 1
                assert self.connections["ul"].blockLevel == self.blockLevel
                assert self.connections["ur"].blockLevel == self.blockLevel
                
                gll_nz = self.connections['ul'].nz
                gf_nz =  self.connections['uf'].nz
                gb_nx =  self.connections['ub'].nx
                gd_nx =  self.connections['uc'].nx
                replaceDic = {'"_gdd_nx_"': str(gdd_nx), '"_gdd_ny_"': str(gdd_ny), '"_gdd_nz_"': str(gdd_nz),
                    '"_gdd_nx_p1_"': str(gdd_nx + 1), '"_gdd_ny_p1_"': str(gdd_ny + 1), '"_gdd_nz_p1_"': str(gdd_nz + 1),
                    '"_gdd_ny_m1_"': str(gdd_ny - 1), '"_gdd_nz_m1_"': str(gdd_nz - 1),
                    '"_gdd_x0_"': str(gdd_r0[0]), '"_gdd_y0_"': str(gdd_r0[1]), '"_gdd_z0_"': str(gdd_r0[2]),
                    '"_gdd_x1_"': str(gdd_r1[0]), '"_gdd_y1_"': str(gdd_r1[1]), '"_gdd_z1_"': str(gdd_r1[2]),
                    '"_gdd_dt_dx_"': str(gdd_dt/gdd_dx), '"_gdd_dt_dy_"': str(gdd_dt/gdd_dy), '"_gdd_dt_dz_"': str(gdd_dt/gdd_dz),
                    '"_gdd_m_dt_dx_"': str(-gdd_dt/gdd_dx), '"_gdd_m_dt_dy_"': str(-gdd_dt/gdd_dy), '"_gdd_m_dt_dz_"': str(-gdd_dt/gdd_dz),
                    '"_gdd_dt_dx_2_"': str(gdd_dt/gdd_dx/2.0), '"_gdd_dt_dy_2_"': str(gdd_dt/gdd_dy/2.0), '"_gdd_dt_dz_2_"': str(gdd_dt/gdd_dz/2.0),
                    '"_gdd_m_dt_dx_2_"': str(-gdd_dt/gdd_dx/2.0), '"_gdd_m_dt_dy_2_"': str(-gdd_dt/gdd_dy/2.0), 
                    '"_gdd_m_dt_dz_2_"': str(-gdd_dt/gdd_dz/2.0),
                    '"_gdd_dt_"': str(gdd_dt),
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_p1_"': str(gll_nz + 1),
                    '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                    '"_gb_nx2_"': str(int(gb_nx/2)),
                    '"_gb_nx2_p_gd_nx2_"': str(int(gb_nx/2) + int(gd_nx/2)),
                    '"grid_dd"': '"' + self.name + '"',
                    '"grid_d"':  '"' + self.connections["uc"].name + '"',
                    '"grid_rr"': '"' + self.connections["ur"].name + '"',
                    '"grid_ll"': '"' + self.connections["ul"].name + '"',
                    '"grid_f"':  '"' + self.connections["uf"].name + '"',
                    '"grid_b"':  '"' + self.connections["ub"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                grid_dd = grids[self.name]
                
                self.SetupSourcesMateriasAndViews(grid_dd)

                connections = self.connections

                #------------------- ff ------------------
                #----------------------------------------
                if "f" in connections:
                    json_file = open('layer1/grid_d/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _gd_nx, _gd_ny, _gd_nz = self.nx, self.ny, self.nz
                    gff_replaceDic = {
                        '"_gd_nx_"': str(_gd_nx), '"_gd_nx_p1_"': str(_gd_nx + 1),
                        '"_gd_ny_"': str(_gd_ny), '"_gd_ny_p1_"': str(_gd_ny + 1),
                        '"_gd_nz_"': str(_gd_nz), '"_gd_nz_p1_"': str(_gd_nz + 1),
                        '"grid_f"': '"' + connections["f"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gff_replaceDic)
                    grid_ff_conn = json.loads(file_content)
                    
                    grid_dd["updateInstructions"].extend(grid_ff_conn["updateInstructions"])
                    
                    for gff_updateSequence in grid_ff_conn["updateSequences"]:
                        seqName = gff_updateSequence["name"]
                        for updateSequence in grid_dd["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gff_updateSequence["sequence"])

                #------------------- bb ------------------
                #----------------------------------------
                if "b" in connections:
                    json_file = open('layer1/grid_d/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    _gd_nx, _gd_ny, _gd_nz = self.nx, self.ny, self.nz
                    _gb_nx = connections["b"].nx
                    gbb_replaceDic = {
                        '"_gd_ny_"': str(_gd_ny), '"_gd_ny_p1_"': str(_gd_ny + 1),
                        '"_gd_nz_"': str(_gd_nz), '"_gd_nz_p1_"': str(_gd_nz + 1),
                        '"_gb_nx_"': str(_gb_nx),
                        '"grid_b"': '"' + connections["b"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gbb_replaceDic)
                    grid_bb_conn = json.loads(file_content)
                    
                    grid_dd["updateInstructions"].extend(grid_bb_conn["updateInstructions"])
                    
                    for gbb_updateSequence in grid_bb_conn["updateSequences"]:
                        seqName = gbb_updateSequence["name"]
                        for updateSequence in grid_dd["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gbb_updateSequence["sequence"])

                #------------------- rrr -----------------
                #-----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_dd)
                    
                #------------------- lll ------------------
                #------------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_dd)

                #------------------- ddd -----------------
                #-----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_dd)

                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_dd)
                    
                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_dd)

                ##------ grid_dd
                return grids
        
        ##---------------------------------------------- grid_f ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "f":
            if self.blockLevel >= 1:
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
                
                self.SetupSourcesMateriasAndViews(grid_f)
        
                connections = self.connections

                #------------------- rr ------------------
                #-----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_f)
                
                                
                #------------------- ll ------------------
                #-----------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_f)
                                
                #------------------- uu ------------------
                #----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_f)
                    
                #------------------- dd ------------------
                #----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_f)
                
                #------------------- ff ------------------
                #-----------------------------------------
                if "f" in connections:
                    json_file = open('layer0/grid_m/connections/grid_f.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    ## comment: grid_m's nx and grid_f's nx are not the same, but nx = gf_nx does the job when it is 
                    ## supposed to select the last elements of grid_f along x
                    _nx, _ny, _nz = self.nx, self.ny, self.nz
                    _gd_ny = connections["d"].ny
                    _gl_nz = connections["l"].nz
                    gff_replaceDic = {
                        '"_nx_"': str(_nx), '"_nx_p1_"': str(_nx + 1),  
                        '"_ny_"': str(_ny), '"_ny_p1_"': str(_ny + 1),  
                        '"_nz_"': str(_nz), '"_nz_p1_"': str(_nz + 1),  
                        '"_nx_m1_"': str(_nx - 1), '"_nx_m2_"': str(_nx - 2),
                        '"_ny_m1_"': str(_ny - 1),
                        '"_nz_m1_"': str(_nz - 1), '"_nz_m2_"': str(_nz - 2), 
                        '"_gl_nz_"': str(_gl_nz), '"_gl_nz_m1_"': str(_gl_nz - 1), '"_gl_nz_p1_"': str(_gl_nz + 1), 
                        '"_gd_ny_"': str(_gd_ny), '"_gd_ny_p1_"': str(_gd_ny + 1), '"_gd_ny_m1_"': str(_gd_ny - 1),
                        '"grid_f"': '"' + connections["f"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gff_replaceDic)
                    grid_ff_conn = json.loads(file_content)
                    
                    grid_f["updateInstructions"].extend(grid_ff_conn["updateInstructions"])                
                                                
                    for gff_updateSequence in grid_ff_conn["updateSequences"]:
                        seqName = gff_updateSequence["name"]
                        for updateSequence in grid_f["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gff_updateSequence["sequence"])

                #------------------- pml-f ----------------
                #------------------------------------------
                if "pml-f" in connections:
                    self.ConnectToFrontPml(grid_f)

                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_f)

                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_f)

                ##---- grid_f
                return grids

        ##---------------------------------------------- grid_b ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "b":
            if self.blockLevel >= 1:
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
                
                self.SetupSourcesMateriasAndViews(grid_b)
        
                connections = self.connections

                #------------------- rr ------------------
                #-----------------------------------------
                if "r" in connections:
                    self.ConnectToRight(grid_b)
                                                
                #------------------- ll ------------------
                #-----------------------------------------
                if "l" in connections:
                    self.ConnectToLeft(grid_b)
                    
                #------------------- uu ------------------
                #-----------------------------------------
                if "u" in connections:
                    self.ConnectToUp(grid_b)
                    
                #------------------- dd ------------------
                #-----------------------------------------
                if "d" in connections:
                    self.ConnectToDown(grid_b)
                

                #------------------- bb ------------------
                #-----------------------------------------
                if "b" in connections:
                    json_file = open('layer0/grid_m/connections/grid_b.json')
                    file_content = json_file.read()
                    json_file.close()

                    ## comment: grid_m's nx and grid_b's nx are not the same
                    _ny, _nz = self.ny, self.nz                    
                    _gb_nx = connections["b"].nx
                    _gd_ny = connections["d"].ny
                    _gl_nz = connections["l"].nz
                    gbb_replaceDic = {
                        '"_ny_"': str(_ny), '"_ny_p1_"': str(_ny + 1),  
                        '"_nz_"': str(_nz), '"_nz_p1_"': str(_nz + 1),  
                        '"_ny_m1_"': str(_ny - 1),
                        '"_nz_m1_"': str(_nz - 1), '"_nz_m2_"': str(_nz - 2), 
                        '"_gl_nz_"': str(_gl_nz), '"_gl_nz_m1_"': str(_gl_nz - 1), '"_gl_nz_p1_"': str(_gl_nz + 1), 
                        '"_gd_ny_"': str(_gd_ny), '"_gd_ny_p1_"': str(_gd_ny + 1), '"_gd_ny_m1_"': str(_gd_ny - 1),
                        '"_gb_nx_"': str(_gb_nx),
                        '"grid_b"': '"' + connections["b"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gbb_replaceDic)
                    grid_bb_conn = json.loads(file_content)

                    grid_b["updateInstructions"].extend(grid_bb_conn["updateInstructions"])                
                                                
                    for gbb_updateSequence in grid_bb_conn["updateSequences"]:
                        seqName = gbb_updateSequence["name"]
                        for updateSequence in grid_b["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gbb_updateSequence["sequence"])


                #------------------- pml-b ------------------
                #-----------------------------------------
                if "pml-b" in connections:
                    self.ConnectToBackPml(grid_b)

                #------------------- pml-r ----------------
                #------------------------------------------
                if "pml-r" in connections:
                    self.ConnectToRightPml(grid_b)

                #------------------- pml-l ----------------
                #------------------------------------------
                if "pml-l" in connections:
                    self.ConnectToLeftPml(grid_b)

                ##---- grid_b
                return grids



    def SetupPMLGrid(self):
        ##---------------------------------------------- pml_f ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        if self.blockPosition == "f":
            if self.blockLevel >= 1:
                json_file = open('pmls/pml_f/pml_f.json')
                file_content = json_file.read()
                json_file.close()
                
                pf_nx, pf_ny, pf_nz = self.nx, self.ny, self.nz
                pf_r0, pf_r1 = self.r0, self.r1
                pf_dt, pf_dx, pf_dy, pf_dz = self.dt, self.dx, self.dy, self.dz 
                assert "b" in self.connections
                assert self.connections["b"].blockLevel == self.blockLevel
        
                pf_sig_ex, pf_sig_hx = self.sig_e, self.sig_h
                replaceDic = {
                    '"_pmlf_nx_"': str(pf_nx), '"_pmlf_ny_"': str(pf_ny), '"_pmlf_nz_"': str(pf_nz),
                    '"_pmlf_nx_p1_"': str(pf_nx + 1), '"_pmlf_ny_p1_"': str(pf_ny + 1), '"_pmlf_nz_p1_"': str(pf_nz + 1),
                    '"_pmlf_nz_m1_"': str(pf_nz - 1),
                    '"_pmlf_x0_"': str(pf_r0[0]), '"_pmlf_y0_"': str(pf_r0[1]), '"_pmlf_z0_"': str(pf_r0[2]),
                    '"_pmlf_x1_"': str(pf_r1[0]), '"_pmlf_y1_"': str(pf_r1[1]), '"_pmlf_z1_"': str(pf_r1[2]),
                    '"_pmlf_dt_dx_"': str(pf_dt/pf_dx), '"_pmlf_dt_dy_"': str(pf_dt/pf_dy), '"_pmlf_dt_dz_"': str(pf_dt/pf_dz),
                    '"_pmlf_m_dt_dx_"': str(-pf_dt/pf_dx), '"_pmlf_m_dt_dy_"': str(-pf_dt/pf_dy), '"_pmlf_m_dt_dz_"': str(-pf_dt/pf_dz),
                    '"_pmlf_dt_"': str(pf_dt),
                    '"_pmlf_dt_sigEx_"': str(pf_dt*pf_sig_ex), '"_pmlf_m_dt_sigEx_"': str(-pf_dt*pf_sig_ex),
                    '"_pmlf_dt_sigHx_"': str(pf_dt*pf_sig_hx), '"_pmlf_m_dt_sigHx_"': str(-pf_dt*pf_sig_hx),
                    '"pml_f"': '"' + self.name + '"',
                    '"grid_f"': '"' + self.connections["b"].name + '"',
                    '"_gf_nx_m1_"': str(self.connections["b"].nx - 1)
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                pml_f = grids[self.name]
                
                ##----- pml_f --
                return grids

        ##---------------------------------------------- pml_b ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "b":
            if self.blockLevel >= 1:
                json_file = open('pmls/pml_b/pml_b.json')
                file_content = json_file.read()
                json_file.close()
                
                pb_nx, pb_ny, pb_nz = self.nx, self.ny, self.nz
                pb_r0, pb_r1 = self.r0, self.r1
                pb_dt, pb_dx, pb_dy, pb_dz = self.dt, self.dx, self.dy, self.dz 
                assert "f" in self.connections
                assert self.connections["f"].blockLevel == self.blockLevel
        
                pb_sig_ex, pb_sig_hx = self.sig_e, self.sig_h
                replaceDic = {
                    '"_pmlb_nx_"': str(pb_nx), '"_pmlb_ny_"': str(pb_ny), '"_pmlb_nz_"': str(pb_nz),
                    '"_pmlb_nx_p1_"': str(pb_nx + 1), '"_pmlb_ny_p1_"': str(pb_ny + 1), '"_pmlb_nz_p1_"': str(pb_nz + 1),
                    '"_pmlb_nx_m1_"': str(pb_nx - 1),
                    '"_pmlb_nz_m1_"': str(pb_nz - 1),
                    '"_pmlb_x0_"': str(pb_r0[0]), '"_pmlb_y0_"': str(pb_r0[1]), '"_pmlb_z0_"': str(pb_r0[2]),
                    '"_pmlb_x1_"': str(pb_r1[0]), '"_pmlb_y1_"': str(pb_r1[1]), '"_pmlb_z1_"': str(pb_r1[2]),
                    '"_pmlb_dt_dx_"': str(pb_dt/pb_dx), '"_pmlb_dt_dy_"': str(pb_dt/pb_dy), '"_pmlb_dt_dz_"': str(pb_dt/pb_dz),
                    '"_pmlb_m_dt_dx_"': str(-pb_dt/pb_dx), '"_pmlb_m_dt_dy_"': str(-pb_dt/pb_dy), '"_pmlb_m_dt_dz_"': str(-pb_dt/pb_dz),
                    '"_pmlb_dt_"': str(pb_dt),
                    '"_pmlb_dt_sigEx_"': str(pb_dt*pb_sig_ex), '"_pmlb_m_dt_sigEx_"': str(-pb_dt*pb_sig_ex),
                    '"_pmlb_dt_sigHx_"': str(pb_dt*pb_sig_hx), '"_pmlb_m_dt_sigHx_"': str(-pb_dt*pb_sig_hx),
                    '"pml_b"': '"' + self.name + '"',
                    '"grid_b"': '"' + self.connections["f"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                pml_b = grids[self.name]
                
                ##----- pml_b --
                return grids

        ##---------------------------------------------- pml_r ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "r":
            if self.blockLevel >= 1:
                json_file = open('pmls/pml_r/pml_r.json')
                file_content = json_file.read()
                json_file.close()
                
                pr_nx, pr_ny, pr_nz = self.nx, self.ny, self.nz
                pr_r0, pr_r1 = self.r0, self.r1
                pr_dt, pr_dx, pr_dy, pr_dz = self.dt, self.dx, self.dy, self.dz 
                assert "lc" in self.connections
                assert "lu" in self.connections
                assert "ld" in self.connections
                assert "lf" in self.connections
                assert "lb" in self.connections
                assert self.connections["lc"].blockLevel == self.blockLevel
                assert self.connections["lu"].blockLevel == self.blockLevel
                assert self.connections["ld"].blockLevel == self.blockLevel
                assert self.connections["lf"].blockLevel == self.blockLevel
                assert self.connections["lb"].blockLevel == self.blockLevel
        
                gr_nx, gr_ny, gr_nz = self.connections["lc"].nx, self.connections["lc"].ny, self.connections["lc"].nz
                gb_nx, gb_nz = self.connections["lb"].nx, self.connections["lb"].nz
                gd_ny, gd_nz = self.connections["ld"].ny, self.connections["ld"].nz
                gf_nz = self.connections["lf"].nz
                gb_nz = self.connections["lb"].nz
                gu_nz = self.connections["lu"].nz
                
                pr_sig_ez, pr_sig_hz = self.sig_e, self.sig_h
                replaceDic = {
                    '"_pmlr_nx_"': str(pr_nx), '"_pmlr_ny_"': str(pr_ny), '"_pmlr_nz_"': str(pr_nz),
                    '"_pmlr_nx_p1_"': str(pr_nx + 1), '"_pmlr_ny_p1_"': str(pr_ny + 1), '"_pmlr_nz_p1_"': str(pr_nz + 1),
                    '"_pmlr_nz_m1_"': str(pr_nz - 1),
                    '"_pmlr_x0_"': str(pr_r0[0]), '"_pmlr_y0_"': str(pr_r0[1]), '"_pmlr_z0_"': str(pr_r0[2]),
                    '"_pmlr_x1_"': str(pr_r1[0]), '"_pmlr_y1_"': str(pr_r1[1]), '"_pmlr_z1_"': str(pr_r1[2]),
                    '"_pmlr_dt_dx_"': str(pr_dt/pr_dx), '"_pmlr_dt_dy_"': str(pr_dt/pr_dy), '"_pmlr_dt_dz_"': str(pr_dt/pr_dz),
                    '"_pmlr_m_dt_dx_"': str(-pr_dt/pr_dx), '"_pmlr_m_dt_dy_"': str(-pr_dt/pr_dy), '"_pmlr_m_dt_dz_"': str(-pr_dt/pr_dz),
                    '"_pmlr_dt_"': str(pr_dt),
                    '"_pmlr_dt_sigEz_"': str(pr_dt*pr_sig_ez), '"_pmlr_m_dt_sigEz_"': str(-pr_dt*pr_sig_ez),
                    '"_pmlr_dt_sigHz_"': str(pr_dt*pr_sig_hz), '"_pmlr_m_dt_sigHz_"': str(-pr_dt*pr_sig_hz),
                    '"_gr_nz_m1_"': str(gr_nz - 1), '"_gu_nz_m1_"': str(gu_nz - 1), '"_gd_nz_m1_"': str(gd_nz - 1), 
                    '"_gf_nz_m1_"': str(gf_nz - 1), '"_gb_nz_m1_"': str(gb_nz - 1),
                    '"_gb_nx_"': str(gb_nx),
                    '"_gb_nx_p_gr_nx_"': str(gb_nx + gr_nx),
                    '"_gd_ny_"': str(gd_ny),
                    '"_gd_ny_p_gr_ny_"': str(gd_ny + gr_ny),
                    '"pml_r"': '"' + self.name + '"',
                    '"grid_r"': '"' + self.connections["lc"].name + '"',
                    '"grid_f"': '"' + self.connections["lf"].name + '"',
                    '"grid_b"': '"' + self.connections["lb"].name + '"',
                    '"grid_u"': '"' + self.connections["lu"].name + '"',
                    '"grid_d"': '"' + self.connections["ld"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                pml_r = grids[self.name]
                
                ##----- pml_r --
                return grids

        ##---------------------------------------------- pml_l ------------------------------------------------
        ##------------------------------------------------------------------------------------------------------
        elif self.blockPosition == "l":
            if self.blockLevel >= 1:
                json_file = open('pmls/pml_l/pml_l.json')
                file_content = json_file.read()
                json_file.close()
                
                pl_nx, pl_ny, pl_nz = self.nx, self.ny, self.nz
                pl_r0, pl_r1 = self.r0, self.r1
                pl_dt, pl_dx, pl_dy, pl_dz = self.dt, self.dx, self.dy, self.dz 
                assert "rc" in self.connections
                assert "ru" in self.connections
                assert "rd" in self.connections
                assert "rf" in self.connections
                assert "rb" in self.connections
                assert self.connections["rc"].blockLevel == self.blockLevel
                assert self.connections["ru"].blockLevel == self.blockLevel
                assert self.connections["rd"].blockLevel == self.blockLevel
                assert self.connections["rf"].blockLevel == self.blockLevel
                assert self.connections["rb"].blockLevel == self.blockLevel
        
                gl_nx, gl_ny, gl_nz = self.connections["rc"].nx, self.connections["rc"].ny, self.connections["rc"].nz
                gr_nx, gr_ny = gl_nx, gl_ny
                gb_nx, gb_nz = self.connections["rb"].nx, self.connections["rb"].nz
                gd_ny, gd_nz = self.connections["rd"].ny, self.connections["rd"].nz
                gf_nz = self.connections["rf"].nz
                gb_nz = self.connections["rb"].nz
                gu_nz = self.connections["ru"].nz
                
                pl_sig_ez, pl_sig_hz = self.sig_e, self.sig_h
                replaceDic = {
                    '"_pmll_nx_"': str(pl_nx), '"_pmll_ny_"': str(pl_ny), '"_pmll_nz_"': str(pl_nz),
                    '"_pmll_nx_p1_"': str(pl_nx + 1), '"_pmll_ny_p1_"': str(pl_ny + 1), '"_pmll_nz_p1_"': str(pl_nz + 1),
                    '"_pmll_nz_m1_"': str(pl_nz - 1),
                    '"_pmll_x0_"': str(pl_r0[0]), '"_pmll_y0_"': str(pl_r0[1]), '"_pmll_z0_"': str(pl_r0[2]),
                    '"_pmll_x1_"': str(pl_r1[0]), '"_pmll_y1_"': str(pl_r1[1]), '"_pmll_z1_"': str(pl_r1[2]),
                    '"_pmll_dt_dx_"': str(pl_dt/pl_dx), '"_pmll_dt_dy_"': str(pl_dt/pl_dy), '"_pmll_dt_dz_"': str(pl_dt/pl_dz),
                    '"_pmll_m_dt_dx_"': str(-pl_dt/pl_dx), '"_pmll_m_dt_dy_"': str(-pl_dt/pl_dy), '"_pmll_m_dt_dz_"': str(-pl_dt/pl_dz),
                    '"_pmll_dt_"': str(pl_dt),
                    '"_pmll_dt_sigEz_"': str(pl_dt*pl_sig_ez), '"_pmll_m_dt_sigEz_"': str(-pl_dt*pl_sig_ez),
                    '"_pmll_dt_sigHz_"': str(pl_dt*pl_sig_hz), '"_pmll_m_dt_sigHz_"': str(-pl_dt*pl_sig_hz),
                    '"_gb_nx_"': str(gb_nx),
                    '"_gb_nx_p_gr_nx_"': str(gb_nx + gr_nx),
                    '"_gd_ny_"': str(gd_ny),
                    '"_gd_ny_p_gr_ny_"': str(gd_ny + gr_ny),
                    '"pml_l"': '"' + self.name + '"',
                    '"grid_l"': '"' + self.connections["rc"].name + '"',
                    '"grid_f"': '"' + self.connections["rf"].name + '"',
                    '"grid_b"': '"' + self.connections["rb"].name + '"',
                    '"grid_u"': '"' + self.connections["ru"].name + '"',
                    '"grid_d"': '"' + self.connections["rd"].name + '"'
                    }
                file_content = MultiWordReplace(file_content, replaceDic)
                grids = json.loads(file_content)
                pml_l = grids[self.name]
                
                ##----- pml_l --
                return grids
        else:
            assert False



    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------- right connections ----------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    def ConnectToRight(self, grid_):
        connections = self.connections
        assert "r" in connections
        if self.blockLevel == 0:
            ##---------------------------------------------- grid_m ------------------------------------------------
            ##------------------------------------------------------------------------------------------------------
            if self.blockPosition == "c":
                json_file = open('layer0/grid_m/connections/grid_r.json')
                file_content = json_file.read()
                json_file.close()
                
                grBlock = connections["r"]
                gr_nx, gr_ny = grBlock.nx, grBlock.ny
                gr_replaceDic = {
                        '"_nx_"': str(self.nx), '"_ny_"': str(self.ny), '"_nz_"': str(self.nz),
                        '"_nx_m1_"': str(self.nx - 1), '"_ny_m1_"': str(self.ny - 1), '"_nz_m1_"': str(self.nz - 1),
                        '"_nx_p1_"': str(self.nx + 1), '"_ny_p1_"': str(self.ny + 1), '"_nz_p1_"': str(self.nz + 1),
                        '"_nx_m2_"': str(self.nx - 2), '"_ny_m2_"': str(self.ny - 2), 
                        '"_gr_nx_m1_"': str(gr_nx - 1), '"_gr_ny_m1_"': str(gr_ny - 1), 
                        '"grid_r"': '"' + connections["r"].name + '"'
                        }
                
                file_content = MultiWordReplace(file_content, gr_replaceDic)
                grid_r_conn = json.loads(file_content)
                
                if "u" in connections or "d" in connections:
                    guBlock = connections["u"]
                    gdBlock = connections["d"]
                    gu_nx, gu_ny, gu_nz = guBlock.nx, guBlock.ny, guBlock.nz
                    gd_ny = gdBlock.ny
                    gl_nz = connections["l"].nz
                    gr_ud_replaceDic = {
                        '"_gd_ny_m1_"': str(gd_ny - 1),
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(self.nz/2)),
                        '"grid_u"': '"' + connections["u"].name + '"',
                        '"grid_d"': '"' + connections["d"].name + '"'
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
                        '"grid_f"': '"' + connections["f"].name + '"',
                        '"grid_b"': '"' + connections["b"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gr_fb_replaceDic)
                    grid_r_conn = json.loads(file_content)
                    
                
                grid_["updateInstructions"].extend(grid_r_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_["updateInstructions"].extend(grid_r_conn["updateInstructions"]["right_only"])
                elif "f" not in connections and "b" not in connections:
                    grid_["updateInstructions"].extend(grid_r_conn["updateInstructions"]["right_up_down"])
                else:
                    grid_["updateInstructions"].extend(grid_r_conn["updateInstructions"]["right_up_down_front_back"])
                
                for gr_updateSequence in grid_r_conn["updateSequences"]:
                    seqName = gr_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gr_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(gr_updateSequence["sequence"]["right_only"])
                            elif "f" not in connections and "b" not in connections:
                                updateSequence["sequence"].extend(gr_updateSequence["sequence"]["right_up_down"])
                            else:
                                updateSequence["sequence"].extend(gr_updateSequence["sequence"]["right_up_down_front_back"])
            else:
                assert False
        elif self.blockLevel >= 1:
            ##---------------------------------------------- grid_r, grid_rr ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            if self.blockPosition == "r":
                json_file = open('layer1/grid_r/connections/grid_rr.json')
                file_content = json_file.read()
                json_file.close()
                
                gb_nx = connections["b"].nx
                gd_ny = connections["d"].ny
                grr_replaceDic = {
                    '"_gr_nx_"': str(self.nx), '"_gr_ny_"': str(self.ny), '"_gr_nz_"': str(self.nz),
                    '"_gr_nx_m1_"': str(self.nx - 1), '"_gr_ny_m1_"': str(self.ny - 1), '"_gr_nz_m1_"': str(self.nz - 1),
                    '"_gr_nx_p1_"': str(self.nx + 1), '"_gr_ny_p1_"': str(self.ny + 1), '"_gr_nz_p1_"': str(self.nz + 1),
                    '"_gd_ny2_"': str(int(gd_ny/2)), '"_gd_ny2_m1_"': str(int(gd_ny/2) - 1), '"_gd_ny2_p1_"': str(int(gd_ny/2) + 1),
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                    '"grid_rr"': '"' + connections["r"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, grr_replaceDic)
                grid_rr_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"])
                
                for grr_updateSequence in grid_rr_conn["updateSequences"]:
                    seqName = grr_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(grr_updateSequence["sequence"])

            ##---------------------------------------------- grid_u, grid_uu ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "u":
                json_file = open('layer1/grid_u/connections/grid_rr.json')
                file_content = json_file.read()
                json_file.close()
                
                grrBlock = connections["r"]
                grr_ny = grrBlock.ny
                gb_nx = connections["b"].nx
                gd_ny = self.ny
                gr_ny = connections["dr"].ny
                grr_replaceDic = {
                    '"_gu_nx_"': str(self.nx), '"_gu_ny_"': str(self.ny), '"_gu_nz_"': str(self.nz),
                    '"_gu_nx_m1_"': str(self.nx - 1), '"_gu_ny_m1_"': str(self.ny - 1), '"_gu_nz_m1_"': str(self.nz - 1),
                    '"_gu_nx_p1_"': str(self.nx + 1), '"_gu_ny_p1_"': str(self.ny + 1), '"_gu_nz_p1_"': str(self.nz + 1),
                    '"_gu_ny_m2_"': str(self.ny - 2),
                    '"_grr_ny_m1_"': str(grr_ny - 1), 
                    '"_gd_ny2_p_gr_ny2_"': str(int(gd_ny/2) + int(gr_ny/2)), 
                    '"_gd_ny2_p_gr_ny2_m1_"': str(int(gd_ny/2) + int(gr_ny/2) - 1), 
                    '"_gd_ny2_p_gr_ny2_p1_"': str(int(gd_ny/2) + int(gr_ny/2) + 1), 
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                    '"grid_rr"': '"' + connections["r"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, grr_replaceDic)
                grid_rr_conn = json.loads(file_content)

                if 'u' in connections:
                    gll_nz = connections['l'].nz
                    gf_nz = connections['f'].nz
                    g_ud_replaceDic = {
                        '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                        '"grid_uu"': '"' + connections['u'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                if "u" not in connections:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                else:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up"])
                
                
                for grr_updateSequence in grid_rr_conn["updateSequences"]:
                    seqName = grr_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                            if "u" not in connections:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                            else:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up"])                                
            

            ##---------------------------------------------- grid_d, grid_dd ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "d":
                json_file = open('layer1/grid_d/connections/grid_rr.json')
                file_content = json_file.read()
                json_file.close()
                
                gb_nx = connections["b"].nx
                grr_replaceDic = {
                    '"_gd_nx_"': str(self.nx), '"_gd_ny_"': str(self.ny), '"_gd_nz_"': str(self.nz),
                    '"_gd_nx_m1_"': str(self.nx - 1), '"_gd_ny_m1_"': str(self.ny - 1), '"_gd_nz_m1_"': str(self.nz - 1),
                    '"_gd_nx_p1_"': str(self.nx + 1), '"_gd_ny_p1_"': str(self.ny + 1), '"_gd_nz_p1_"': str(self.nz + 1),
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                    '"grid_rr"': '"' + connections["r"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, grr_replaceDic)
                grid_rr_conn = json.loads(file_content)
                
                if 'd' in connections:
                    gdd_ny = connections['d'].ny
                    gll_nz = connections['l'].nz
                    gf_nz = connections['f'].nz
                    g_ud_replaceDic = {
                        '"_gdd_ny_m1_"': str(gdd_ny - 1),
                        '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                        '"grid_dd"': '"' + connections['d'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_rr_conn = json.loads(file_content)

                grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                if "d" not in connections:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                else:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_down"])
                
                for grr_updateSequence in grid_rr_conn["updateSequences"]:
                    seqName = grr_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                            if "d" not in connections:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                            else:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_down"])                                


            ##---------------------------------------------- grid_f, grid_ff ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "f":
                json_file = open('layer1/grid_f/connections/grid_rr.json')
                file_content = json_file.read()
                json_file.close()
                
                grrBlock = connections["r"]
                grr_nx, grr_ny = grrBlock.nx, grrBlock.ny
                gb_nx = self.nx
                gr_nx = connections["br"].nx
                grr_replaceDic = {
                    '"_gf_nx_"': str(self.nx), '"_gf_ny_"': str(self.ny), '"_gf_nz_"': str(self.nz), 
                    '"_gf_nx_m1_"': str(self.nx - 1), '"_gf_ny_m1_"': str(self.ny - 1), '"_gf_nz_m1_"': str(self.nz - 1), 
                    '"_gf_nx_p1_"': str(self.nx + 1), '"_gf_ny_p1_"': str(self.ny + 1), '"_gf_nz_p1_"': str(self.nz + 1), 
                    '"_gf_nx_m2_"': str(self.nx - 2), '"_gf_ny_m2_"': str(self.ny - 2),
                    '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)),
                    '"_gb_nx2_p_gr_nx2_p1_"': str(int(gb_nx/2) + int(gr_nx/2) + 1),
                    '"_gb_nx2_p_gr_nx2_m1_"': str(int(gb_nx/2) + int(gr_nx/2) - 1),
                    '"_grr_nx_m1_"': str(grr_nx - 1),
                    '"_grr_ny_m1_"': str(grr_ny - 1),
                    '"grid_rr"': '"' + connections["r"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, grr_replaceDic)
                grid_rr_conn = json.loads(file_content)

                if 'u' in connections or 'd' in connections:
                    gdd_ny = connections['d'].ny
                    gll_nz = connections['l'].nz
                    gf_nz = self.nz
                    g_ud_replaceDic = {
                        '"_gdd_ny_m1_"': str(gdd_ny - 1),
                        '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                        '"grid_uu"': '"' + connections['u'].name + '"',
                        '"grid_dd"': '"' + connections['d'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                if 'f' in connections:
                    gdd_ny = connections['d'].ny
                    g_ff_replaceDic = {
                        '"_gdd_ny_"': str(gdd_ny), 
                        '"grid_ff"': '"' + connections['f'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ff_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                elif "f" not in connections:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down"])
                else:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down_front"])
                
                for grr_updateSequence in grid_rr_conn["updateSequences"]:
                    seqName = grr_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                            elif "f" not in connections:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down"])  
                            else:                              
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down_front"])  
                            

            ##---------------------------------------------- grid_b, grid_bb ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "b":
                json_file = open('layer1/grid_b/connections/grid_rr.json')
                file_content = json_file.read()
                json_file.close()
                
                grrBlock = connections["r"]
                grr_ny = grrBlock.ny
                grr_replaceDic = {
                    '"_gb_nx_"': str(self.nx), '"_gb_ny_"': str(self.ny), '"_gb_nz_"': str(self.nz), 
                    '"_gb_nx_m1_"': str(self.nx - 1), '"_gb_ny_m1_"': str(self.ny - 1), '"_gb_nz_m1_"': str(self.nz - 1), 
                    '"_gb_nx_p1_"': str(self.nx + 1), '"_gb_ny_p1_"': str(self.ny + 1), '"_gb_nz_p1_"': str(self.nz + 1), 
                    '"_gb_ny_m2_"': str(self.ny - 2),
                    '"_grr_ny_m1_"': str(grr_ny - 1),
                    '"grid_rr"': '"' + connections["r"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, grr_replaceDic)
                grid_rr_conn = json.loads(file_content)
                
                if 'u' in connections or 'd' in connections:
                    gdd_ny = connections['d'].ny
                    gll_nz = connections['l'].nz
                    gf_nz = self.nz
                    g_ud_replaceDic = {
                        '"_gdd_ny_m1_"': str(gdd_ny - 1),
                        '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                        '"grid_uu"': '"' + connections['u'].name + '"',
                        '"grid_dd"': '"' + connections['d'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                if 'b' in connections:
                    gbb_nx = connections['b'].nx
                    gdd_ny = connections['d'].ny
                    g_bb_replaceDic = {
                        '"_gbb_nx_m1_"': str(gbb_nx - 1),
                        '"_gdd_ny_"': str(gdd_ny),
                        '"grid_bb"': '"' + connections['b'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_bb_replaceDic)
                    grid_rr_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_only"])
                elif "b" not in connections:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down"])
                else:
                    grid_["updateInstructions"].extend(grid_rr_conn["updateInstructions"]["right_up_down_back"])
                
                for grr_updateSequence in grid_rr_conn["updateSequences"]:
                    seqName = grr_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(grr_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_only"])
                            elif "b" not in connections:
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down"])  
                            else:                              
                                updateSequence["sequence"].extend(grr_updateSequence["sequence"]["right_up_down_back"])  
                            
            else:
                assert False

    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##---------------------------------------------- left connections --------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    def ConnectToLeft(self, grid_):
        connections = self.connections
        assert "l" in connections
        if self.blockLevel == 0:
            ##---------------------------------------------- grid_m ------------------------------------------------
            ##------------------------------------------------------------------------------------------------------
            if self.blockPosition == "c":
                json_file = open('layer0/grid_m/connections/grid_l.json')
                file_content = json_file.read()
                json_file.close()
                
                glBlock = connections["l"]
                gl_nx, gl_ny, gl_nz = glBlock.nx, glBlock.ny, glBlock.nz
                gl_replaceDic = {
                    '"_nx_"': str(self.nx), '"_ny_"': str(self.ny), '"_nz_"': str(self.nz),
                    '"_nx_m1_"': str(self.nx - 1), '"_ny_m1_"': str(self.ny - 1), '"_nz_m1_"': str(self.nz - 1),
                    '"_nx_p1_"': str(self.nx + 1), '"_ny_p1_"': str(self.ny + 1), '"_nz_p1_"': str(self.nz + 1),
                    '"_gl_nx_m1_"': str(gl_nx - 1), '"_gl_ny_m1_"': str(gl_ny - 1), 
                    '"_gl_nz_"': str(gl_nz),
                    '"_nx_m2_"': str(self.nx - 2), '"_ny_m2_"': str(self.ny - 2), 
                    '"grid_l"': '"' + connections["l"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gl_replaceDic)
                grid_l_conn = json.loads(file_content)

                if "u" in connections or "d" in connections:
                    gd_ny = connections["d"].ny
                    gl_ud_replaceDic = {
                        '"_gd_ny_m1_"': str(gd_ny - 1),
                        '"_gd_ny_p_ny2_"': str(gd_ny + int(self.ny/2)),
                        '"grid_u"': '"' + connections["u"].name + '"',
                        '"grid_d"': '"' + connections["d"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gl_ud_replaceDic)
                    grid_l_conn = json.loads(file_content)

                if "f" in connections or "b" in connections:
                    gb_nx = connections["b"].nx
                    gd_ny = connections["d"].ny
                    gl_fb_replaceDic = {
                        '"_gd_ny_"': str(gd_ny),
                        '"_gb_nx_m1_"': str(gb_nx - 1),
                        '"grid_f"': '"' + connections["f"].name + '"',
                        '"grid_b"': '"' + connections["b"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gl_fb_replaceDic)
                    grid_l_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_l_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_["updateInstructions"].extend(grid_l_conn["updateInstructions"]["left_only"])
                elif "f" not in connections and "b" not in connections:
                    grid_["updateInstructions"].extend(grid_l_conn["updateInstructions"]["left_up_down"])
                else:
                    grid_["updateInstructions"].extend(grid_l_conn["updateInstructions"]["left_up_down_front_back"])
                
                for gl_updateSequence in grid_l_conn["updateSequences"]:
                    seqName = gl_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gl_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(gl_updateSequence["sequence"]["left_only"])
                            elif "f" not in connections and "b" not in connections:
                                updateSequence["sequence"].extend(gl_updateSequence["sequence"]["left_up_down"])
                            else:
                                updateSequence["sequence"].extend(gl_updateSequence["sequence"]["left_up_down_front_back"])
            
            else:
                assert False
            
        elif self.blockLevel >= 1:
            ##---------------------------------------------- grid_l, grid_ll ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            if self.blockPosition == "l":
                json_file = open('layer1/grid_l/connections/grid_ll.json')
                file_content = json_file.read()
                json_file.close()
                
                gllBlock = connections["l"]
                gll_nz = gllBlock.nz
                gb_nx = connections["b"].nx
                gd_ny = connections["d"].ny
                gll_replaceDic = {
                    '"_gl_nx_"': str(self.nx), '"_gl_ny_"': str(self.ny), '"_gl_nz_"': str(self.nz),
                    '"_gl_nx_m1_"': str(self.nx - 1), '"_gl_ny_m1_"': str(self.ny - 1), '"_gl_nz_m1_"': str(self.nz - 1),
                    '"_gl_nx_p1_"': str(self.nx + 1), '"_gl_ny_p1_"': str(self.ny + 1), '"_gl_nz_p1_"': str(self.nz + 1),
                    '"_gll_nz_"': str(gll_nz),
                    '"_gd_ny2_"': str(int(gd_ny/2)), '"_gd_ny2_m1_"': str(int(gd_ny/2) - 1), '"_gd_ny2_p1_"': str(int(gd_ny/2) + 1),
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                    '"grid_ll"': '"' + gllBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gll_replaceDic)
                grid_ll_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"])
                
                for gll_updateSequence in grid_ll_conn["updateSequences"]:
                    seqName = gll_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gll_updateSequence["sequence"])
            
            ##---------------------------------------------- grid_u, grid_uu ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "u":
                json_file = open('layer1/grid_u/connections/grid_ll.json')
                file_content = json_file.read()
                json_file.close()
                
                gllBlock = connections["l"]
                gll_ny, gll_nz = gllBlock.ny, gllBlock.nz
                gb_nx = connections["b"].nx
                gd_ny = self.ny
                gr_ny = connections["dr"].ny
                gll_replaceDic = {
                    '"_gu_nx_"': str(self.nx), '"_gu_ny_"': str(self.ny), '"_gu_nz_"': str(self.nz),
                    '"_gu_nx_m1_"': str(self.nx - 1), '"_gu_ny_m1_"': str(self.ny - 1), '"_gu_nz_m1_"': str(self.nz - 1),
                    '"_gu_nx_p1_"': str(self.nx + 1), '"_gu_ny_p1_"': str(self.ny + 1), '"_gu_nz_p1_"': str(self.nz + 1),
                    '"_gu_ny_m2_"': str(self.ny - 2),
                    '"_gll_ny_m1_"': str(gll_ny - 1), 
                    '"_gll_nz_"': str(gll_nz),
                    '"_gd_ny2_p_gr_ny2_"': str(int(gd_ny/2) + int(gr_ny/2)), 
                    '"_gd_ny2_p_gr_ny2_m1_"': str(int(gd_ny/2) + int(gr_ny/2) - 1), 
                    '"_gd_ny2_p_gr_ny2_p1_"': str(int(gd_ny/2) + int(gr_ny/2) + 1), 
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                    '"grid_ll"': '"' + gllBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gll_replaceDic)
                grid_ll_conn = json.loads(file_content)

                if 'u' in connections:
                    g_ud_replaceDic = {
                        '"grid_uu"': '"' + connections['u'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                if "u" not in connections:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                else:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up"])
                
                for gll_updateSequence in grid_ll_conn["updateSequences"]:
                    seqName = gll_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                            if "u" not in connections:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                            else:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up"])                                
            
            ##---------------------------------------------- grid_d, grid_dd ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "d":
                json_file = open('layer1/grid_d/connections/grid_ll.json')
                file_content = json_file.read()
                json_file.close()
                
                gllBlock = connections["l"]
                gll_ny, gll_nz = gllBlock.ny, gllBlock.nz
                gb_nx = connections["b"].nx
                gll_replaceDic = {
                    '"_gd_nx_"': str(self.nx), '"_gd_ny_"': str(self.ny), '"_gd_nz_"': str(self.nz),
                    '"_gd_nx_m1_"': str(self.nx - 1), '"_gd_ny_m1_"': str(self.ny - 1), '"_gd_nz_m1_"': str(self.nz - 1),
                    '"_gd_nx_p1_"': str(self.nx + 1), '"_gd_ny_p1_"': str(self.ny + 1), '"_gd_nz_p1_"': str(self.nz + 1),
                    '"_gll_ny_m1_"': str(gll_ny - 1), 
                    '"_gll_nz_"': str(gll_nz),
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1),
                    '"grid_ll"': '"' + gllBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gll_replaceDic)
                grid_ll_conn = json.loads(file_content)

                if 'd' in connections:
                    gdd_ny = connections['d'].ny
                    g_ud_replaceDic = {
                        '"_gdd_ny_m1_"': str(gdd_ny - 1),
                        '"grid_dd"': '"' + connections['d'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                if "d" not in connections:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                else:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_down"])
                
                for gll_updateSequence in grid_ll_conn["updateSequences"]:
                    seqName = gll_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                            if "d" not in connections:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                            else:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_down"])                                

            ##---------------------------------------------- grid_f, grid_ff ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "f":
                json_file = open('layer1/grid_f/connections/grid_ll.json')
                file_content = json_file.read()
                json_file.close()
                
                gllBlock = connections["l"]
                gll_nx, gll_ny, gll_nz = gllBlock.nx, gllBlock.ny, gllBlock.nz
                gb_nx = self.nx
                gr_nx = connections["br"].nx
                gll_replaceDic = {
                    '"_gf_nx_"': str(self.nx), '"_gf_ny_"': str(self.ny), '"_gf_nz_"': str(self.nz), 
                    '"_gf_nx_m1_"': str(self.nx - 1), '"_gf_ny_m1_"': str(self.ny - 1), '"_gf_nz_m1_"': str(self.nz - 1), 
                    '"_gf_nx_p1_"': str(self.nx + 1), '"_gf_ny_p1_"': str(self.ny + 1), '"_gf_nz_p1_"': str(self.nz + 1), 
                    '"_gf_nx_m2_"': str(self.nx - 2),
                    '"_gf_ny_m2_"': str(self.ny - 2),
                    '"_gll_nx_m1_"': str(gll_nx - 1), 
                    '"_gll_ny_m1_"': str(gll_ny - 1), 
                    '"_gll_nz_"': str(gll_nz),
                    '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)),
                    '"_gb_nx2_p_gr_nx2_p1_"': str(int(gb_nx/2) + int(gr_nx/2) + 1),
                    '"_gb_nx2_p_gr_nx2_m1_"': str(int(gb_nx/2) + int(gr_nx/2) - 1),
                    '"grid_ll"': '"' + gllBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gll_replaceDic)
                grid_ll_conn = json.loads(file_content)

                if 'u' in connections or 'd' in connections:
                    gdd_ny = connections['d'].ny
                    gf_nx = self.nx
                    g_ud_replaceDic = {
                        '"_gf_nx_m1_"': str(gf_nx - 1),
                        '"_gdd_ny_m1_"': str(gdd_ny - 1),
                        '"grid_uu"': '"' + connections['u'].name + '"',
                        '"grid_dd"': '"' + connections['d'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                if 'f' in connections:
                    gdd_ny = connections['d'].ny
                    g_ff_replaceDic = {
                        '"_gdd_ny_"': str(gdd_ny), 
                        '"grid_ff"': '"' + connections['f'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ff_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                elif "f" not in connections:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down"])
                else:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down_front"])
                
                for gll_updateSequence in grid_ll_conn["updateSequences"]:
                    seqName = gll_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                            elif "f" not in connections:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down"])
                            else:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down_front"])

            ##---------------------------------------------- grid_b, grid_bb ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "b":
                json_file = open('layer1/grid_b/connections/grid_ll.json')
                file_content = json_file.read()
                json_file.close()
                
                gllBlock = connections["l"]
                gll_ny, gll_nz = gllBlock.ny, gllBlock.nz
                gll_replaceDic = {
                    '"_gb_nx_"': str(self.nx), '"_gb_ny_"': str(self.ny), '"_gb_nz_"': str(self.nz),
                    '"_gb_nx_m1_"': str(self.nx - 1), '"_gb_ny_m1_"': str(self.ny - 1), '"_gb_nz_m1_"': str(self.nz - 1),
                    '"_gb_nx_p1_"': str(self.nx + 1), '"_gb_ny_p1_"': str(self.ny + 1), '"_gb_nz_p1_"': str(self.nz + 1),
                    '"_gb_ny_m2_"': str(self.ny - 2),
                    '"_gll_ny_m1_"': str(gll_ny - 1), 
                    '"_gll_nz_"': str(gll_nz),
                    '"grid_ll"': '"' + gllBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gll_replaceDic)
                grid_ll_conn = json.loads(file_content)

                if 'u' in connections or 'd' in connections:
                    gdd_ny = connections['d'].ny
                    g_ud_replaceDic = {
                        '"_gdd_ny_m1_"': str(gdd_ny - 1),
                        '"grid_uu"': '"' + connections['u'].name + '"',
                        '"grid_dd"': '"' + connections['d'].name + '"',
                    }
                    file_content = MultiWordReplace(file_content, g_ud_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                if 'b' in connections:
                    gbb_nx = connections['b'].nx
                    gdd_ny = connections['d'].ny
                    g_bb_replaceDic = {
                        '"_gbb_nx_m1_"': str(gbb_nx - 1), 
                        '"_gdd_ny_"': str(gdd_ny),
                        '"grid_bb"': '"' + connections['b'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_bb_replaceDic)
                    grid_ll_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["general"])
                if "u" not in connections and "d" not in connections:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_only"])
                elif "b" not in connections:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down"])
                else:
                    grid_["updateInstructions"].extend(grid_ll_conn["updateInstructions"]["left_up_down_back"])
                
                for gll_updateSequence in grid_ll_conn["updateSequences"]:
                    seqName = gll_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gll_updateSequence["sequence"]["general"])
                            if "u" not in connections and "d" not in connections:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_only"])
                            elif "b" not in connections:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down"])
                            else:
                                updateSequence["sequence"].extend(gll_updateSequence["sequence"]["left_up_down_back"])

            else:
                assert False

    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##---------------------------------------------- up connections ----------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    def ConnectToUp(self, grid_):
        connections = self.connections
        assert "u" in connections
        
        if self.blockLevel == 0:
            ##---------------------------------------------- grid_m ------------------------------------------------
            ##------------------------------------------------------------------------------------------------------
            assert self.blockPosition == "c"
            json_file = open('layer0/grid_m/connections/grid_u.json')
            file_content = json_file.read()
            json_file.close()
            
            guBlock = connections["u"]
            gu_nx, gu_ny = guBlock.nx, guBlock.ny
            gl_nz = connections["l"].nz
            gu_replaceDic = {
                '"_nx_"': str(self.nx), '"_ny_"': str(self.ny), '"_nz_"': str(self.nz), 
                '"_nx_m1_"': str(self.nx - 1), '"_ny_m1_"': str(self.ny - 1), '"_nz_m1_"': str(self.nz - 1), 
                '"_nx_p1_"': str(self.nx + 1), '"_ny_p1_"': str(self.ny + 1), '"_nz_p1_"': str(self.nz + 1), 
                '"_nx_m2_"': str(self.nx - 2),
                '"_nz_m2_"': str(self.nz - 2), 
                '"_gl_nz_"': str(gl_nz), '"_gl_nz_m1_"': str(gl_nz - 1), '"_gl_nz_p1_"': str(gl_nz + 1), 
                '"_gu_nx_m1_"': str(gu_nx - 1),
                '"grid_u"': '"' + guBlock.name + '"'
                }
            
            file_content = MultiWordReplace(file_content, gu_replaceDic)
            grid_u_conn = json.loads(file_content)
            
            if "f" in connections or "b" in connections:
                gb_nx = connections["b"].nx
                gd_ny = connections["d"].ny
                gu_fb_replaceDic = {
                    '"_gd_ny_"': str(gd_ny),
                    '"_gd_ny_p_ny2_"': str(gd_ny + int(self.ny/2)),
                    '"_gb_nx_m1_"': str(gb_nx - 1),
                    '"grid_f"': '"' + connections["f"].name + '"',
                    '"grid_b"': '"' + connections["b"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gu_fb_replaceDic)
                grid_u_conn = json.loads(file_content)
            

            grid_["updateInstructions"].extend(grid_u_conn["updateInstructions"]["general"])
            if "f" not in connections and "b" not in connections:
                grid_["updateInstructions"].extend(grid_u_conn["updateInstructions"]["up_only"])
            else:
                grid_["updateInstructions"].extend(grid_u_conn["updateInstructions"]["up_front_back"])
            
                                        
            for gu_updateSequence in grid_u_conn["updateSequences"]:
                seqName = gu_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(gu_updateSequence["sequence"]["general"])
                        if "f" not in connections and "b" not in connections:
                            updateSequence["sequence"].extend(gu_updateSequence["sequence"]["up_only"])
                        else:
                            updateSequence["sequence"].extend(gu_updateSequence["sequence"]["up_front_back"])

        elif self.blockLevel >= 1:
            if self.blockPosition == "r":    
                ##---------------------------------------------- grid_r ------------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                if self.blockLevel == 1:
                    json_file = open('layer1/grid_r/connections/grid_u.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gl_nz = self.nz
                    nz = connections["l"].nz
                    gu_replaceDic = {
                        '"_gr_nx_"': str(self.nx), '"_gr_ny_"': str(self.ny), '"_gr_nz_"': str(self.nz),
                        '"_gr_nx_m1_"': str(self.nx - 1), '"_gr_ny_m1_"': str(self.ny - 1), '"_gr_nz_m1_"': str(self.nz - 1),
                        '"_gr_nx_p1_"': str(self.nx + 1), '"_gr_ny_p1_"': str(self.ny + 1), '"_gr_nz_p1_"': str(self.nz + 1),
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                        '"grid_u"': '"' + connections["u"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gu_replaceDic)
                    grid_u_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_u_conn["updateInstructions"])
                    
                    for gu_updateSequence in grid_u_conn["updateSequences"]:
                        seqName = gu_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gu_updateSequence["sequence"])

                ##---------------------------------------------- grid_rr -----------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                else:
                    json_file = open('layer2/grid_rr/connections/grid_uu.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gll_nz = self.nz
                    gf_nz = connections['lf'].nz
                    guu_replaceDic = {
                        '"_grr_nx_"': str(self.nx), '"_grr_ny_"': str(self.ny), '"_grr_nz_"': str(self.nz),
                        '"_grr_nx_m1_"': str(self.nx - 1), '"_grr_ny_m1_"': str(self.ny - 1), '"_grr_nz_m1_"': str(self.nz - 1),
                        '"_grr_nx_p1_"': str(self.nx + 1), '"_grr_ny_p1_"': str(self.ny + 1), '"_grr_nz_p1_"': str(self.nz + 1),
                        '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                        '"grid_uu"': '"' + connections["u"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, guu_replaceDic)
                    grid_uu_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"])
                    
                    for guu_updateSequence in grid_uu_conn["updateSequences"]:
                        seqName = guu_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(guu_updateSequence["sequence"])

            elif self.blockPosition == "l":    
                ##---------------------------------------------- grid_l ------------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                if self.blockLevel == 1:
                    json_file = open('layer1/grid_l/connections/grid_u.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gu_replaceDic = {
                        '"_gl_nx_"': str(self.nx), '"_gl_ny_"': str(self.ny), '"_gl_nz_"': str(self.nz),
                        '"_gl_nx_m1_"': str(self.nx - 1), '"_gl_ny_m1_"': str(self.ny - 1), '"_gl_nz_m1_"': str(self.nz - 1),
                        '"_gl_nx_p1_"': str(self.nx + 1), '"_gl_ny_p1_"': str(self.ny + 1), '"_gl_nz_p1_"': str(self.nz + 1),
                        '"grid_u"': '"' + connections["u"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gu_replaceDic)
                    grid_u_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_u_conn["updateInstructions"])
                    
                    for gu_updateSequence in grid_u_conn["updateSequences"]:
                        seqName = gu_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gu_updateSequence["sequence"])

                ##---------------------------------------------- grid_ll ------------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                else:
                    json_file = open('layer2/grid_ll/connections/grid_uu.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    guu_replaceDic = {
                        '"_gll_nx_"': str(self.nx), '"_gll_ny_"': str(self.ny), '"_gll_nz_"': str(self.nz),
                        '"_gll_nx_m1_"': str(self.nx - 1), '"_gll_ny_m1_"': str(self.ny - 1), '"_gll_nz_m1_"': str(self.nz - 1),
                        '"_gll_nx_p1_"': str(self.nx + 1), '"_gll_ny_p1_"': str(self.ny + 1), '"_gll_nz_p1_"': str(self.nz + 1),
                        '"grid_uu"': '"' + connections["u"].name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, guu_replaceDic)
                    grid_uu_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"])
                    
                    for guu_updateSequence in grid_uu_conn["updateSequences"]:
                        seqName = guu_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(guu_updateSequence["sequence"])

            ##---------------------------------------------- grid_u, grid_uu ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "u":    
                json_file = open('layer1/grid_u/connections/grid_uu.json')
                file_content = json_file.read()
                json_file.close()
                
                gll_nz = connections['l'].nz
                gb_nx = connections['b'].nx
                guu_replaceDic = {
                    '"_gu_nx_"': str(self.nx), '"_gu_ny_"': str(self.ny), '"_gu_nz_"': str(self.nz),
                    '"_gu_nx_m1_"': str(self.nx - 1), '"_gu_ny_m1_"': str(self.ny - 1), '"_gu_nz_m1_"': str(self.nz - 1),
                    '"_gu_nx_p1_"': str(self.nx + 1), '"_gu_ny_p1_"': str(self.ny + 1), '"_gu_nz_p1_"': str(self.nz + 1),
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_m1_"': str(gll_nz - 1), '"_gll_nz_p1_"': str(gll_nz + 1), 
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1), 
                    '"grid_uu"': '"' + connections["u"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, guu_replaceDic)
                grid_uu_conn = json.loads(file_content)

                grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"])
                
                for guu_updateSequence in grid_uu_conn["updateSequences"]:
                    seqName = guu_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(guu_updateSequence["sequence"])

            ##---------------------------------------------- grid_f, grid_ff ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "f":    
                json_file = open('layer1/grid_f/connections/grid_uu.json')
                file_content = json_file.read()
                json_file.close()
                
                guu_nx = connections["u"].nx
                gll_nz = connections["l"].nz
                gb_nx = self.nx
                gr_nx = connections['br'].nx
                guu_replaceDic = {
                    '"_gf_nx_"': str(self.nx), '"_gf_ny_"': str(self.ny), '"_gf_nz_"': str(self.nz), 
                    '"_gf_nx_m1_"': str(self.nx - 1), '"_gf_ny_m1_"': str(self.ny - 1), '"_gf_nz_m1_"': str(self.nz - 1), 
                    '"_gf_nx_p1_"': str(self.nx + 1), '"_gf_ny_p1_"': str(self.ny + 1), '"_gf_nz_p1_"': str(self.nz + 1), 
                    '"_gf_nx_m2_"': str(self.nx - 2),
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_m1_"': str(gll_nz - 1), '"_gll_nz_p1_"': str(gll_nz + 1), 
                    '"_guu_nx_m1_"': str(guu_nx - 1), 
                    '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)), 
                    '"_gb_nx2_p_gr_nx2_m1_"': str(int(gb_nx/2) + int(gr_nx/2) - 1), 
                    '"_gb_nx2_p_gr_nx2_p1_"': str(int(gb_nx/2) + int(gr_nx/2) + 1), 
                    '"grid_uu"': '"' + connections["u"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, guu_replaceDic)
                grid_uu_conn = json.loads(file_content)
                
                if 'f' in connections:
                    gdd_ny = connections['d'].ny
                    gf_ny = self.ny
                    g_ff_replaceDic = {
                        '"_gdd_ny_p_gf_ny2_"': str(gdd_ny + int(gf_ny/2)),
                        '"grid_ff"': '"' + connections['f'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ff_replaceDic)
                    grid_uu_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"]["general"])
                if "f" not in connections:
                    grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"]["up_only"])
                else:
                    grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"]["up_front"])
                
                for guu_updateSequence in grid_uu_conn["updateSequences"]:
                    seqName = guu_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(guu_updateSequence["sequence"]["general"])
                            if "f" not in connections:
                                updateSequence["sequence"].extend(guu_updateSequence["sequence"]["up_only"])
                            else:
                                updateSequence["sequence"].extend(guu_updateSequence["sequence"]["up_front"])

            ##---------------------------------------------- grid_b, grid_bb ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "b":    
                json_file = open('layer1/grid_b/connections/grid_uu.json')
                file_content = json_file.read()
                json_file.close()
                
                gll_nz = connections['l'].nz
                guu_replaceDic = {
                    '"_gb_nx_"': str(self.nx), '"_gb_ny_"': str(self.ny), '"_gb_nz_"': str(self.nz), 
                    '"_gb_nx_m1_"': str(self.nx - 1), '"_gb_ny_m1_"': str(self.ny - 1), '"_gb_nz_m1_"': str(self.nz - 1), 
                    '"_gb_nx_p1_"': str(self.nx + 1), '"_gb_ny_p1_"': str(self.ny + 1), '"_gb_nz_p1_"': str(self.nz + 1), 
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_m1_"': str(gll_nz - 1), '"_gll_nz_p1_"': str(gll_nz + 1), 
                    '"grid_uu"': '"' + connections["u"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, guu_replaceDic)
                grid_uu_conn = json.loads(file_content)
                
                if 'b' in connections:
                    gbb_nx = connections['b'].nx
                    gdd_ny = connections['d'].ny
                    gf_ny = self.ny
                    g_bb_replaceDic = {
                        '"_gbb_nx_m1_"': str(gbb_nx - 1), 
                        '"_gdd_ny_p_gf_ny2_"': str(gdd_ny + int(gf_ny/2)),
                        '"grid_bb"': '"' + connections['b'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_bb_replaceDic)
                    grid_uu_conn = json.loads(file_content)
                
                
                grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"]["general"])
                if "b" not in connections:
                    grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"]["up_only"])
                else:
                    grid_["updateInstructions"].extend(grid_uu_conn["updateInstructions"]["up_back"])
                
                
                for guu_updateSequence in grid_uu_conn["updateSequences"]:
                    seqName = guu_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(guu_updateSequence["sequence"]["general"])
                            if "b" not in connections:
                                updateSequence["sequence"].extend(guu_updateSequence["sequence"]["up_only"])
                            else:
                                updateSequence["sequence"].extend(guu_updateSequence["sequence"]["up_back"])
        else:
            assert False 

    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##-------------------------------------------- down connections ----------------------------------------
    ##------------------------------------------------------------------------------------------------------
    ##------------------------------------------------------------------------------------------------------
    def ConnectToDown(self, grid_):
        connections = self.connections
        assert "d" in connections
        
        if self.blockLevel == 0:
            ##---------------------------------------------- grid_m ------------------------------------------------
            ##------------------------------------------------------------------------------------------------------
            assert self.blockPosition == "c"
            json_file = open('layer0/grid_m/connections/grid_d.json')
            file_content = json_file.read()
            json_file.close()
            
            gdBlock = connections["d"]
            gd_nx, gd_ny = gdBlock.nx, gdBlock.ny
            gl_nz = connections["l"].nz
            gd_replaceDic = {
                '"_nx_"': str(self.nx), '"_ny_"': str(self.ny), '"_nz_"': str(self.nz), 
                '"_nx_m1_"': str(self.nx - 1), '"_ny_m1_"': str(self.ny - 1), '"_nz_m1_"': str(self.nz - 1), 
                '"_nx_p1_"': str(self.nx + 1), '"_ny_p1_"': str(self.ny + 1), '"_nz_p1_"': str(self.nz + 1), 
                '"_nx_m2_"': str(self.nx - 2),
                '"_gl_nz_"': str(gl_nz), '"_gl_nz_m1_"': str(gl_nz - 1), '"_gl_nz_p1_"': str(gl_nz + 1),
                '"_gd_ny_"': str(gd_ny), '"_gd_nx_m1_"': str(gd_nx - 1),
                '"grid_d"': '"' + gdBlock.name + '"'
                }
            
            file_content = MultiWordReplace(file_content, gd_replaceDic)
            grid_d_conn = json.loads(file_content)
            
            if "f" in connections or "b" in connections:
                gb_nx = connections["b"].nx
                gd_ny = connections["d"].ny
                gd_fb_replaceDic = {
                    '"_gd_ny_"': str(gd_ny),
                    '"_gb_nx_m1_"': str(gb_nx - 1),
                    '"grid_f"': '"' + connections["f"].name + '"',
                    '"grid_b"': '"' + connections["b"].name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gd_fb_replaceDic)
                grid_d_conn = json.loads(file_content)
            

            grid_["updateInstructions"].extend(grid_d_conn["updateInstructions"]["general"])
            if "f" not in connections and "b" not in connections:
                grid_["updateInstructions"].extend(grid_d_conn["updateInstructions"]["down_only"])
            else:
                grid_["updateInstructions"].extend(grid_d_conn["updateInstructions"]["down_front_back"])
            
                                        
            for gd_updateSequence in grid_d_conn["updateSequences"]:
                seqName = gd_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(gd_updateSequence["sequence"]["general"])
                        if "f" not in connections and "b" not in connections:
                            updateSequence["sequence"].extend(gd_updateSequence["sequence"]["down_only"])
                        else:
                            updateSequence["sequence"].extend(gd_updateSequence["sequence"]["down_front_back"])

        elif self.blockLevel >= 1:
            if self.blockPosition == "r":    
                ##---------------------------------------------- grid_r ------------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                if self.blockLevel == 1:
                    json_file = open('layer1/grid_r/connections/grid_d.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gdBlock = connections["d"]
                    gd_ny = gdBlock.ny
                    gl_nz = self.nz
                    nz = connections["l"].nz
                    gd_replaceDic = {
                        '"_gr_nx_"': str(self.nx), '"_gr_ny_"': str(self.ny), '"_gr_nz_"': str(self.nz), 
                        '"_gr_nx_m1_"': str(self.nx - 1), '"_gr_ny_m1_"': str(self.ny - 1), '"_gr_nz_m1_"': str(self.nz - 1), 
                        '"_gr_nx_p1_"': str(self.nx + 1), '"_gr_ny_p1_"': str(self.ny + 1), '"_gr_nz_p1_"': str(self.nz + 1), 
                        '"_gd_ny_"': str(gd_ny),
                        '"_gl_nz_p_nz2_"': str(gl_nz + int(nz/2)),
                        '"grid_d"': '"' + gdBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gd_replaceDic)
                    grid_d_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_d_conn["updateInstructions"])
                    
                    for gd_updateSequence in grid_d_conn["updateSequences"]:
                        seqName = gd_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gd_updateSequence["sequence"])
                    
                ##---------------------------------------------- grid_rr -----------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                else:
                    json_file = open('layer2/grid_rr/connections/grid_dd.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gddBlock = connections["d"]
                    gdd_ny = gddBlock.ny
                    gll_nz = self.nz
                    gf_nz = connections['lf'].nz
                    gdd_replaceDic = {
                        '"_grr_nx_"': str(self.nx), '"_grr_ny_"': str(self.ny), '"_grr_nz_"': str(self.nz), 
                        '"_grr_nx_m1_"': str(self.nx - 1), '"_grr_ny_m1_"': str(self.ny - 1), '"_grr_nz_m1_"': str(self.nz - 1), 
                        '"_grr_nx_p1_"': str(self.nx + 1), '"_grr_ny_p1_"': str(self.ny + 1), '"_grr_nz_p1_"': str(self.nz + 1), 
                        '"_gll_nz_p_gf_nz2_"': str(gll_nz + int(gf_nz/2)),
                        '"_gdd_ny_"': str(gdd_ny),
                        '"grid_dd"': '"' + gddBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gdd_replaceDic)
                    grid_dd_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"])
                    
                    for gdd_updateSequence in grid_dd_conn["updateSequences"]:
                        seqName = gdd_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gdd_updateSequence["sequence"])

            elif self.blockPosition == "l":    
                ##---------------------------------------------- grid_l ------------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                if self.blockLevel == 1:
                    json_file = open('layer1/grid_l/connections/grid_d.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gdBlock = connections["d"]
                    gd_ny = gdBlock.ny
                    gl_nz = self.nz
                    gd_replaceDic = {
                        '"_gl_nx_"': str(self.nx), '"_gl_ny_"': str(self.ny), '"_gl_nz_"': str(self.nz), 
                        '"_gl_nx_m1_"': str(self.nx - 1), '"_gl_ny_m1_"': str(self.ny - 1), '"_gl_nz_m1_"': str(self.nz - 1), 
                        '"_gl_nx_p1_"': str(self.nx + 1), '"_gl_ny_p1_"': str(self.ny + 1), '"_gl_nz_p1_"': str(self.nz + 1), 
                        '"_gd_ny_"': str(gd_ny),
                        '"grid_d"': '"' + gdBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gd_replaceDic)
                    grid_d_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_d_conn["updateInstructions"])
                    
                    for gd_updateSequence in grid_d_conn["updateSequences"]:
                        seqName = gd_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gd_updateSequence["sequence"])

                ##---------------------------------------------- grid_ll -----------------------------------------------
                ##------------------------------------------------------------------------------------------------------
                else:
                    json_file = open('layer2/grid_ll/connections/grid_dd.json')
                    file_content = json_file.read()
                    json_file.close()
                    
                    gddBlock = connections["d"]
                    gdd_ny = gddBlock.ny
                    gdd_replaceDic = {
                        '"_gll_nx_"': str(self.nx), '"_gll_ny_"': str(self.ny), '"_gll_nz_"': str(self.nz), 
                        '"_gll_nx_m1_"': str(self.nx - 1), '"_gll_ny_m1_"': str(self.ny - 1), '"_gll_nz_m1_"': str(self.nz - 1), 
                        '"_gll_nx_p1_"': str(self.nx + 1), '"_gll_ny_p1_"': str(self.ny + 1), '"_gll_nz_p1_"': str(self.nz + 1), 
                        '"_gdd_ny_"': str(gdd_ny),
                        '"grid_dd"': '"' + gddBlock.name + '"'
                        }
                    
                    file_content = MultiWordReplace(file_content, gdd_replaceDic)
                    grid_dd_conn = json.loads(file_content)
                    
                    grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"])
                    
                    for gdd_updateSequence in grid_dd_conn["updateSequences"]:
                        seqName = gdd_updateSequence["name"]
                        for updateSequence in grid_["updateSequences"]:
                            if updateSequence["name"] == seqName:
                                updateSequence["sequence"].extend(gdd_updateSequence["sequence"])

            ##---------------------------------------------- grid_d, grid_dd ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "d":
                json_file = open('layer1/grid_d/connections/grid_dd.json')
                file_content = json_file.read()
                json_file.close()
                
                gddBlock = connections["d"]
                gdd_ny = gddBlock.ny
                gll_nz = connections['l'].nz
                gb_nx =  connections['b'].nx
                gdd_replaceDic = {
                    '"_gd_nx_"': str(self.nx), '"_gd_ny_"': str(self.ny), '"_gd_nz_"': str(self.nz), 
                    '"_gd_nx_m1_"': str(self.nx - 1), '"_gd_ny_m1_"': str(self.ny - 1), '"_gd_nz_m1_"': str(self.nz - 1), 
                    '"_gd_nx_p1_"': str(self.nx + 1), '"_gd_ny_p1_"': str(self.ny + 1), '"_gd_nz_p1_"': str(self.nz + 1), 
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_m1_"': str(gll_nz - 1), '"_gll_nz_p1_"': str(gll_nz + 1),
                    '"_gb_nx2_"': str(int(gb_nx/2)), '"_gb_nx2_p1_"': str(int(gb_nx/2) + 1), '"_gb_nx2_m1_"': str(int(gb_nx/2) - 1),
                    '"_gdd_ny_"': str(gdd_ny),
                    '"grid_dd"': '"' + gddBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gdd_replaceDic)
                grid_dd_conn = json.loads(file_content)
                                    
                grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"])
                
                for gdd_updateSequence in grid_dd_conn["updateSequences"]:
                    seqName = gdd_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gdd_updateSequence["sequence"])
            
            
            ##---------------------------------------------- grid_f, grid_ff ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "f":
                json_file = open('layer1/grid_f/connections/grid_dd.json')
                file_content = json_file.read()
                json_file.close()
                
                gddBlock = connections["d"]
                gdd_nx, gdd_ny = gddBlock.nx, gddBlock.ny
                gll_nz = connections["l"].nz
                gb_nx = self.nx
                gr_nx = connections['br'].nx
                gdd_replaceDic = {
                    '"_gf_nx_"': str(self.nx), '"_gf_ny_"': str(self.ny), '"_gf_nz_"': str(self.nz), 
                    '"_gf_nx_m1_"': str(self.nx - 1), '"_gf_ny_m1_"': str(self.ny - 1), '"_gf_nz_m1_"': str(self.nz - 1), 
                    '"_gf_nx_p1_"': str(self.nx + 1), '"_gf_ny_p1_"': str(self.ny + 1), '"_gf_nz_p1_"': str(self.nz + 1), 
                    '"_gf_nx_m2_"': str(self.nx - 2),
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_m1_"': str(gll_nz - 1), '"_gll_nz_p1_"': str(gll_nz + 1), 
                    '"_gdd_ny_"': str(gdd_ny), 
                    '"_gdd_nx_m1_"': str(gdd_nx - 1), 
                    '"_gb_nx2_p_gr_nx2_"': str(int(gb_nx/2) + int(gr_nx/2)), 
                    '"_gb_nx2_p_gr_nx2_m1_"': str(int(gb_nx/2) + int(gr_nx/2) - 1), 
                    '"_gb_nx2_p_gr_nx2_p1_"': str(int(gb_nx/2) + int(gr_nx/2) + 1), 
                    '"grid_dd"': '"' + gddBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gdd_replaceDic)
                grid_dd_conn = json.loads(file_content)

                if 'f' in connections:
                    g_ff_replaceDic = {
                        '"grid_ff"': '"' + connections['f'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_ff_replaceDic)
                    grid_dd_conn = json.loads(file_content)
                
                grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"]["general"])
                if "f" not in connections:
                    grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"]["down_only"])
                else:
                    grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"]["down_front"])
                
                for gdd_updateSequence in grid_dd_conn["updateSequences"]:
                    seqName = gdd_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gdd_updateSequence["sequence"]["general"])
                            if "f" not in connections:
                                updateSequence["sequence"].extend(gdd_updateSequence["sequence"]["down_only"])
                            else:
                                updateSequence["sequence"].extend(gdd_updateSequence["sequence"]["down_front"])

            ##---------------------------------------------- grid_b, grid_bb ---------------------------------------
            ##------------------------------------------------------------------------------------------------------
            elif self.blockPosition == "b":
                json_file = open('layer1/grid_b/connections/grid_dd.json')
                file_content = json_file.read()
                json_file.close()
                
                gddBlock = connections["d"]
                gdd_ny = gddBlock.ny
                gll_nz = connections['l'].nz
                gdd_replaceDic = {
                    '"_gb_nx_"': str(self.nx), '"_gb_ny_"': str(self.ny), '"_gb_nz_"': str(self.nz), 
                    '"_gb_nx_m1_"': str(self.nx - 1), '"_gb_ny_m1_"': str(self.ny - 1), '"_gb_nz_m1_"': str(self.nz - 1), 
                    '"_gb_nx_p1_"': str(self.nx + 1), '"_gb_ny_p1_"': str(self.ny + 1), '"_gb_nz_p1_"': str(self.nz + 1), 
                    '"_gll_nz_"': str(gll_nz), '"_gll_nz_m1_"': str(gll_nz - 1), '"_gll_nz_p1_"': str(gll_nz + 1), 
                    '"_gdd_ny_"': str(gdd_ny),
                    '"grid_dd"': '"' + gddBlock.name + '"'
                    }
                
                file_content = MultiWordReplace(file_content, gdd_replaceDic)
                grid_dd_conn = json.loads(file_content)
                
                if 'b' in connections:
                    gbb_nx = connections['b'].nx
                    g_bb_replaceDic = {
                        '"_gbb_nx_m1_"': str(gbb_nx - 1), 
                        '"grid_bb"': '"' + connections['b'].name + '"'
                    }
                    file_content = MultiWordReplace(file_content, g_bb_replaceDic)
                    grid_dd_conn = json.loads(file_content)
                
                
                grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"]["general"])
                if "b" not in connections:
                    grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"]["down_only"])
                else:
                    grid_["updateInstructions"].extend(grid_dd_conn["updateInstructions"]["down_back"])
                
                for gdd_updateSequence in grid_dd_conn["updateSequences"]:
                    seqName = gdd_updateSequence["name"]
                    for updateSequence in grid_["updateSequences"]:
                        if updateSequence["name"] == seqName:
                            updateSequence["sequence"].extend(gdd_updateSequence["sequence"]["general"])
                            if "b" not in connections:
                                updateSequence["sequence"].extend(gdd_updateSequence["sequence"]["down_only"])
                            else:
                                updateSequence["sequence"].extend(gdd_updateSequence["sequence"]["down_back"])


            else:
                assert False


    def ConnectToFrontPml(self, grid_):
        connections = self.connections
        assert "pml-f" in connections
        assert self.blockLevel >= 1
        assert self.blockPosition == "f"
        assert "f" not in connections
        
        json_file = open('layer1/grid_f/connections/pml_f.json')
        file_content = json_file.read()
        json_file.close()
        
        pmlf_replaceDic = {
            '"_gf_nx_"': str(self.nx), '"_gf_ny_"': str(self.ny), '"_gf_nz_"': str(self.nz), 
            '"_gf_nx_p1_"': str(self.nx + 1), '"_gf_ny_p1_"': str(self.ny + 1), '"_gf_nz_p1_"': str(self.nz + 1), 
            '"pml_f"': '"' + connections["pml-f"].name + '"'
            }
        
        file_content = MultiWordReplace(file_content, pmlf_replaceDic)
        pml_f_conn = json.loads(file_content)
                            
        grid_["updateInstructions"].extend(pml_f_conn["updateInstructions"])
        
        for pf_updateSequence in pml_f_conn["updateSequences"]:
            seqName = pf_updateSequence["name"]
            for updateSequence in grid_["updateSequences"]:
                if updateSequence["name"] == seqName:
                    updateSequence["sequence"].extend(pf_updateSequence["sequence"])
            

    def ConnectToBackPml(self, grid_):
        connections = self.connections
        assert "pml-b" in connections
        assert self.blockLevel >= 1
        assert self.blockPosition == "b"
        assert "b" not in connections
        
        json_file = open('layer1/grid_b/connections/pml_b.json')
        file_content = json_file.read()
        json_file.close()
        
        pmlb_replaceDic = {
            '"_gb_nx_"': str(self.nx), '"_gb_ny_"': str(self.ny), '"_gb_nz_"': str(self.nz), 
            '"_gb_nx_p1_"': str(self.nx + 1), '"_gb_ny_p1_"': str(self.ny + 1), '"_gb_nz_p1_"': str(self.nz + 1), 
            '"pml_b"': '"' + connections["pml-b"].name + '"',
            '"_pmlb_nx_"': str(connections["pml-b"].nx)
            }
        
        file_content = MultiWordReplace(file_content, pmlb_replaceDic)
        pml_b_conn = json.loads(file_content)
                            
        grid_["updateInstructions"].extend(pml_b_conn["updateInstructions"])
        
        for pb_updateSequence in pml_b_conn["updateSequences"]:
            seqName = pb_updateSequence["name"]
            for updateSequence in grid_["updateSequences"]:
                if updateSequence["name"] == seqName:
                    updateSequence["sequence"].extend(pb_updateSequence["sequence"])
            
    
    def ConnectToRightPml(self, grid_):
        connections = self.connections
        assert "pml-r" in connections
        assert self.blockLevel >= 1
        assert self.blockPosition in ["r", "u", "d", "f", "b"]
        assert "r" not in connections
        
        if self.blockPosition == "r":
            json_file = open('layer1/grid_r/connections/pml_r.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.connections["b"].nx
            gd_ny = self.connections["d"].ny

            pmlr_replaceDic = {
                '"_gr_nx_"': str(self.nx), '"_gr_ny_"': str(self.ny), '"_gr_nz_"': str(self.nz), 
                '"_gr_nx_p1_"': str(self.nx + 1), '"_gr_ny_p1_"': str(self.ny + 1), '"_gr_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_"': str(gb_nx), '"_gd_ny_"': str(gd_ny),
                '"pml_r"': '"' + connections["pml-r"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmlr_replaceDic)
            pml_r_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_r_conn["updateInstructions"])
            
            for pr_updateSequence in pml_r_conn["updateSequences"]:
                seqName = pr_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pr_updateSequence["sequence"])

        elif self.blockPosition == "u":
            json_file = open('layer1/grid_u/connections/pml_r.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.connections["b"].nx
            gd_ny = self.ny
            gr_ny = self.connections["dr"].ny
            pmlr_replaceDic = {
                '"_gu_nx_"': str(self.nx), '"_gu_ny_"': str(self.ny), '"_gu_nz_"': str(self.nz), 
                '"_gu_nx_p1_"': str(self.nx + 1), '"_gu_ny_p1_"': str(self.ny + 1), '"_gu_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_"': str(gb_nx), '"_gd_ny_p_gr_ny_"': str(gd_ny + gr_ny),
                '"pml_r"': '"' + connections["pml-r"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmlr_replaceDic)
            pml_r_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_r_conn["updateInstructions"])
            
            for pr_updateSequence in pml_r_conn["updateSequences"]:
                seqName = pr_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pr_updateSequence["sequence"])

        elif self.blockPosition == "d":
            json_file = open('layer1/grid_d/connections/pml_r.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.connections["b"].nx
            pmlr_replaceDic = {
                '"_gd_nx_"': str(self.nx), '"_gd_ny_"': str(self.ny), '"_gd_nz_"': str(self.nz), 
                '"_gd_nx_p1_"': str(self.nx + 1), '"_gd_ny_p1_"': str(self.ny + 1), '"_gd_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_"': str(gb_nx), 
                '"pml_r"': '"' + connections["pml-r"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmlr_replaceDic)
            pml_r_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_r_conn["updateInstructions"])
            
            for pr_updateSequence in pml_r_conn["updateSequences"]:
                seqName = pr_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pr_updateSequence["sequence"])

        elif self.blockPosition == "f":
            json_file = open('layer1/grid_f/connections/pml_r.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.nx
            gr_nx = self.connections["br"].nx
            pmlr_replaceDic = {
                '"_gf_nx_"': str(self.nx), '"_gf_ny_"': str(self.ny), '"_gf_nz_"': str(self.nz), 
                '"_gf_nx_p1_"': str(self.nx + 1), '"_gf_ny_p1_"': str(self.ny + 1), '"_gf_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_p_gr_nx_"': str(gb_nx + gr_nx), 
                '"pml_r"': '"' + connections["pml-r"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmlr_replaceDic)
            pml_r_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_r_conn["updateInstructions"])
            
            for pr_updateSequence in pml_r_conn["updateSequences"]:
                seqName = pr_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pr_updateSequence["sequence"])

        elif self.blockPosition == "b":
            json_file = open('layer1/grid_b/connections/pml_r.json')
            file_content = json_file.read()
            json_file.close()
            
            pmlr_replaceDic = {
                '"_gb_nx_"': str(self.nx), '"_gb_ny_"': str(self.ny), '"_gb_nz_"': str(self.nz), 
                '"_gb_nx_p1_"': str(self.nx + 1), '"_gb_ny_p1_"': str(self.ny + 1), '"_gb_nz_p1_"': str(self.nz + 1), 
                '"pml_r"': '"' + connections["pml-r"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmlr_replaceDic)
            pml_r_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_r_conn["updateInstructions"])
            
            for pr_updateSequence in pml_r_conn["updateSequences"]:
                seqName = pr_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pr_updateSequence["sequence"])
        else:
            assert False


    def ConnectToLeftPml(self, grid_):
        connections = self.connections
        assert "pml-l" in connections
        assert self.blockLevel >= 1
        assert self.blockPosition in ["l", "u", "d", "f", "b"]
        assert "l" not in connections
        
        if self.blockPosition == "l":
            json_file = open('layer1/grid_l/connections/pml_l.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.connections["b"].nx
            gd_ny = self.connections["d"].ny

            pmll_replaceDic = {
                '"_gl_nx_"': str(self.nx), '"_gl_ny_"': str(self.ny), '"_gl_nz_"': str(self.nz), 
                '"_gl_nx_p1_"': str(self.nx + 1), '"_gl_ny_p1_"': str(self.ny + 1), '"_gl_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_"': str(gb_nx), '"_gd_ny_"': str(gd_ny),
                '"_gll_nz_"': str(connections["pml-l"].nz),
                '"pml_l"': '"' + connections["pml-l"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmll_replaceDic)
            pml_l_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_l_conn["updateInstructions"])
            
            for pl_updateSequence in pml_l_conn["updateSequences"]:
                seqName = pl_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pl_updateSequence["sequence"])

        elif self.blockPosition == "u":
            json_file = open('layer1/grid_u/connections/pml_l.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.connections["b"].nx
            gd_ny = self.ny
            gr_ny = self.connections["dr"].ny
            pmll_replaceDic = {
                '"_gu_nx_"': str(self.nx), '"_gu_ny_"': str(self.ny), '"_gu_nz_"': str(self.nz), 
                '"_gu_nx_p1_"': str(self.nx + 1), '"_gu_ny_p1_"': str(self.ny + 1), '"_gu_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_"': str(gb_nx), '"_gd_ny_p_gr_ny_"': str(gd_ny + gr_ny),
                '"_gll_nz_"': str(connections["pml-l"].nz),
                '"pml_l"': '"' + connections["pml-l"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmll_replaceDic)
            pml_l_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_l_conn["updateInstructions"])
            
            for pl_updateSequence in pml_l_conn["updateSequences"]:
                seqName = pl_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pl_updateSequence["sequence"])

        elif self.blockPosition == "d":
            json_file = open('layer1/grid_d/connections/pml_l.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.connections["b"].nx
            pmll_replaceDic = {
                '"_gd_nx_"': str(self.nx), '"_gd_ny_"': str(self.ny), '"_gd_nz_"': str(self.nz), 
                '"_gd_nx_p1_"': str(self.nx + 1), '"_gd_ny_p1_"': str(self.ny + 1), '"_gd_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_"': str(gb_nx), 
                '"_gll_nz_"': str(connections["pml-l"].nz),
                '"pml_l"': '"' + connections["pml-l"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmll_replaceDic)
            pml_l_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_l_conn["updateInstructions"])
            
            for pl_updateSequence in pml_l_conn["updateSequences"]:
                seqName = pl_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pl_updateSequence["sequence"])

        elif self.blockPosition == "f":
            json_file = open('layer1/grid_f/connections/pml_l.json')
            file_content = json_file.read()
            json_file.close()
            
            gb_nx = self.nx
            gr_nx = self.connections["br"].nx
            pmll_replaceDic = {
                '"_gf_nx_"': str(self.nx), '"_gf_ny_"': str(self.ny), '"_gf_nz_"': str(self.nz), 
                '"_gf_nx_p1_"': str(self.nx + 1), '"_gf_ny_p1_"': str(self.ny + 1), '"_gf_nz_p1_"': str(self.nz + 1), 
                '"_gb_nx_p_gr_nx_"': str(gb_nx + gr_nx), 
                '"_gll_nz_"': str(connections["pml-l"].nz),
                '"pml_l"': '"' + connections["pml-l"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmll_replaceDic)
            pml_l_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_l_conn["updateInstructions"])
            
            for pl_updateSequence in pml_l_conn["updateSequences"]:
                seqName = pl_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pl_updateSequence["sequence"])

        elif self.blockPosition == "b":
            json_file = open('layer1/grid_b/connections/pml_l.json')
            file_content = json_file.read()
            json_file.close()
            
            pmll_replaceDic = {
                '"_gb_nx_"': str(self.nx), '"_gb_ny_"': str(self.ny), '"_gb_nz_"': str(self.nz), 
                '"_gb_nx_p1_"': str(self.nx + 1), '"_gb_ny_p1_"': str(self.ny + 1), '"_gb_nz_p1_"': str(self.nz + 1), 
                '"_gll_nz_"': str(connections["pml-l"].nz),
                '"pml_l"': '"' + connections["pml-l"].name + '"'
                }
            
            file_content = MultiWordReplace(file_content, pmll_replaceDic)
            pml_l_conn = json.loads(file_content)
                                
            grid_["updateInstructions"].extend(pml_l_conn["updateInstructions"])
            
            for pl_updateSequence in pml_l_conn["updateSequences"]:
                seqName = pl_updateSequence["name"]
                for updateSequence in grid_["updateSequences"]:
                    if updateSequence["name"] == seqName:
                        updateSequence["sequence"].extend(pl_updateSequence["sequence"])
        else:
            assert False

            


        
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



        
