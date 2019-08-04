
__all__ = ["MultilevelGridsParameters"]

import numpy as np
import json
import pickle

from ParameterFileGenerator import GridBlock, GridCollectionner
from Geometries import Hyperboloid, Cone, Cylinder

class MultilevelGridsParameters:
    
    def __init__(self):
        self.layersNumOfCells = []
        self.geometries = []
        self.materials = []
        self.sources = []
        self.views = []
        self.boxes = []
        self.grids = []
        self.pmls = {}
        self.gridCollectionner = GridCollectionner()
        self.gridsIntersectingGeometries = {}
        self.viewFolder = "3D/auto/"
        self.maxBufferSize = 1024*1024*100
        self.minBufferSize = 1024*1024*1
        
    def SetCenterGridDimensions(self, r0, r1, dr, S, nFactor=8, nIsPowerOf2=False):
        n_cells = ((r1 - r0)/dr).astype("int")
        #n_cells = 2**(np.ceil(np.log2(n_cells)).astype("int"))   ## n = 2**m
        
        if nIsPowerOf2:
            n_cells = 2**(np.ceil(np.log2(n_cells)).astype("int"))   ## n = 2**m
        elif np.any(n_cells % nFactor != 0):
            #make  n_cell a multiple of nFactor
            n_cells += -(n_cells % nFactor) + nFactor
            assert np.all(n_cells % nFactor == 0)
        
        dr = (r1 - r0)/n_cells
        dt = S / (np.linalg.norm(1.0/dr)) 
        self.boxes.append({"r0":r0, "r1":r1, "n":n_cells, "dr":dr, "dt":dt})
        self.layersNumOfCells.append(n_cells)
                        
    def AddLayer(self, numOfCells):
        n_lev = len(self.boxes)
        assert n_lev >= 1
        box_prev = self.boxes[-1]

        #print("n_prev: ", box_prev["n"], "   n_prev%8: ", box_prev["n"]%8)
        assert np.all(box_prev["n"] % 2 == 0)
        assert np.all(numOfCells % 2 == 0)

        assert np.all(box_prev["n"] % 8 == 0)
        n_prev_8 = box_prev["n"] / 8
        #print("n_prev_8: ", n_prev_8)
        numOfCells += (4 - numOfCells%4)*(n_prev_8 % 2 == 0) + (2 - numOfCells%4)*(n_prev_8 % 2 != 0)
        #print("numOfCells: ", numOfCells)

        n = box_prev["n"]/2 + 2*numOfCells
        assert np.all(n % 8 == 0)
        
        self.layersNumOfCells.append(numOfCells)
        dr = box_prev["dr"]*2
        dt = box_prev["dt"]*2
        r0 = box_prev["r0"] - numOfCells*dr
        r1 = box_prev["r1"] + numOfCells*dr
        
        self.boxes.append({"r0":r0, "r1":r1, "n":n, "dr":dr, "dt":dt})
        
    def AddPML(self, face, numOfCells):
        self.pmls[face] = {"n": numOfCells}
                
    def AddGeometry(self, geomParams):
        self.geometries.append(geomParams)
        
    def AddMaterial(self, matParams):
        self.materials.append(matParams)

    def AddSource(self, srcParams):
        self.sources.append(srcParams)
    
    def AddView(self, viewParams):
        self.views.append(viewParams)

    def SetupGeometries(self):
        for geomParams in self.geometries:
            if geomParams["type"] == "hyperboloid":
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        r0 , r1 = grid.r0, grid.r1
                        
                        geom = Hyperboloid(geomParams["coneAngle"], geomParams["tipRadius"], 
                                           geomParams["height"], geomParams["apexPosition"])
                        boundingBox = geom.GetBoundingBox(r0[1], r1[1])
                        
                        if grid.name not in self.gridsIntersectingGeometries:
                            self.gridsIntersectingGeometries[grid.name] = {}

                        if boundingBox is None:
                            continue
                        bb_r0, bb_r1 = boundingBox
                        
                        ## check if the grid overlaps the boundingbox
                        overlap_r0 = np.maximum(r0, bb_r0)
                        overlap_r1 = np.minimum(r1, bb_r1)
                        
                        if np.all(overlap_r0 <= overlap_r1):
                            grid.AddGeometry(
                                {"type":"hyperboloid", 
                                "geometryName": geomParams["geometryName"], "coneAngle": geomParams["coneAngle"],
                                "tipRadius": geomParams["tipRadius"], "height": geomParams["height"],
                                "apexPosition": geomParams["apexPosition"],
                                "boundingBox":[bb_r0, bb_r1]
                                })
                            self.gridsIntersectingGeometries[grid.name][geomParams["geometryName"]] = boundingBox
            elif geomParams["type"] == "cone":
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        r0 , r1 = grid.r0, grid.r1
                        
                        geom = Cone(geomParams["coneAngle"], geomParams["tipRadius"], 
                                    geomParams["height"], geomParams["apexPosition"])
                        boundingBox = geom.GetBoundingBox(r0[1], r1[1])
                        
                        if grid.name not in self.gridsIntersectingGeometries:
                            self.gridsIntersectingGeometries[grid.name] = {}

                        if boundingBox is None:
                            continue
                        bb_r0, bb_r1 = boundingBox
                        #print(grid.name, bb_r0, bb_r1)
                        
                        ## check if the grid overlaps the boundingbox
                        overlap_r0 = np.maximum(r0, bb_r0)
                        overlap_r1 = np.minimum(r1, bb_r1)
                        
                        if np.all(overlap_r0 <= overlap_r1):
                            grid.AddGeometry(
                                {"type":"cone", 
                                "geometryName": geomParams["geometryName"], "coneAngle": geomParams["coneAngle"],
                                "tipRadius": geomParams["tipRadius"], "height": geomParams["height"],
                                "apexPosition": geomParams["apexPosition"],
                                "boundingBox":[bb_r0, bb_r1]
                                })
                            self.gridsIntersectingGeometries[grid.name][geomParams["geometryName"]] = boundingBox

            elif geomParams["type"] == "cylinder":
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        r0 , r1 = grid.r0, grid.r1
                        
                        geom = Cylinder(geomParams["radius"], 
                                        geomParams["height"], geomParams["topCenter"])
                        boundingBox = geom.GetBoundingBox(r0[1], r1[1])
                        
                        if grid.name not in self.gridsIntersectingGeometries:
                            self.gridsIntersectingGeometries[grid.name] = {}

                        if boundingBox is None:
                            continue
                        bb_r0, bb_r1 = boundingBox
                        
                        ## check if the grid overlaps the boundingbox
                        overlap_r0 = np.maximum(r0, bb_r0)
                        overlap_r1 = np.minimum(r1, bb_r1)
                                                
                        if np.all(overlap_r0 <= overlap_r1):
                            alignEven = "no"
                            r_tc = geomParams["topCenter"]
                            if  r0[1] < r_tc[1] < r1[1]  and r_tc[1] - geomParams["height"] < r0[1]:
                                alignEven = "yes"

                            grid.AddGeometry(
                                {"type":"cylinder", 
                                "geometryName": geomParams["geometryName"],
                                "radius": geomParams["radius"], "height": geomParams["height"],
                                "topCenter": geomParams["topCenter"],
                                "alignEven": alignEven,
                                "boundingBox":[bb_r0, bb_r1]
                                })
                            self.gridsIntersectingGeometries[grid.name][geomParams["geometryName"]] = boundingBox

            else:
                assert False            
        print("gridsIntersectingGeometries ", self.gridsIntersectingGeometries)

    def SetupMaterials(self):
        for matParams in self.materials:
            if matParams["type"] == "DrudeMetal_PureScattered":
                geomName = matParams["geometryName"]
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        
                        if geomName in self.gridsIntersectingGeometries[grid.name]:
                            boundingBox = self.gridsIntersectingGeometries[grid.name][geomName]
                            gridMat = {"type": "DrudeMetal_PureScattered",
                            "boundingBox": boundingBox,
                            "geometryName": geomName, "plasmaFrequency": matParams["plasmaFrequency"],
                            "scatteringRate": matParams["scatteringRate"]
                            }
                            if "wireMeshAlong" in matParams:
                                gridMat["wireMeshAlong"] = matParams["wireMeshAlong"]
                            grid.AddMaterial(gridMat)
   
            elif matParams["type"] == "pec_PureScattered":
                geomName = matParams["geometryName"]
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        
                        if geomName in self.gridsIntersectingGeometries[grid.name]:
                            boundingBox = self.gridsIntersectingGeometries[grid.name][geomName]
                            grid.AddMaterial(
                            {"type": "pec_PureScattered",
                            "boundingBox": boundingBox,
                            "geometryName": geomName
                            })
            elif matParams["type"] == "pec":
                geomName = matParams["geometryName"]
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        
                        if geomName in self.gridsIntersectingGeometries[grid.name]:
                            boundingBox = self.gridsIntersectingGeometries[grid.name][geomName]
                            grid.AddMaterial(
                            {"type": "pec",
                            "boundingBox": boundingBox,
                            "geometryName": geomName
                            })
            else:
                assert False            


    def SetupSources(self):
        for srcParams in self.sources:
            if srcParams["type"] == "GaussianPointSource":
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        j_r = srcParams["position"]
                        r0 , r1 = grid.r0, grid.r1
                        dr = np.array([grid.dx, grid.dy, grid.dz])
                        
                        if np.all(np.logical_and(j_r >= r0, j_r < r1)):
                            j_inds = ((j_r - r0)/dr).astype("int")
                            grid.AddSource(
                                {"type":"GaussianPointSource",
                                "ind_x": j_inds[0], "ind_y": j_inds[1], "ind_z": j_inds[2],
                                "polarization": srcParams["polarization"], 
                                "amplitude": srcParams["amplitude"], 
                                "t_center": srcParams["t_center"], 
                                "t_decay": srcParams["t_decay"],
                                "modulationFrequency":srcParams["modulationFrequency"], 
                                "modulationPhase": srcParams["modulationPhase"], 
                                "timeOffsetFraction": srcParams["timeOffsetFraction"]
                                })
            
            elif srcParams["type"] == "GaussianLineSource_y":       ## directed along y
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        j_r = srcParams["position"]     ## center point
                        height = srcParams["height"]
                        
                        r0 , r1 = grid.r0, grid.r1
                        dr = np.array([grid.dx, grid.dy, grid.dz])
                        
                        j_r0 = np.array([j_r[0], j_r[1] - height/2.0, j_r[2]])
                        j_r1 = np.array([j_r[0], j_r[1] + height/2.0, j_r[2]])

                        overlap_r0 = np.maximum(r0, j_r0)
                        overlap_r1 = np.minimum(r1, j_r1)
                        
                        if np.all(overlap_r0 <= overlap_r1):
                            j_inds = ((j_r - r0)/dr).astype("int")
                            grid.AddSource(
                                {"type":"GaussianLineSource_y",
                                "j_r": j_r,
                                "height": height,
                                "polarization": srcParams["polarization"], 
                                "amplitude": srcParams["amplitude"], 
                                "t_center": srcParams["t_center"], 
                                "t_decay": srcParams["t_decay"],
                                "modulationFrequency":srcParams["modulationFrequency"], 
                                "modulationPhase": srcParams["modulationPhase"], 
                                "timeOffsetFraction": srcParams["timeOffsetFraction"]
                                })

            elif srcParams["type"] == "GaussianSheetSource_z":      ## normal to z
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        j_r = srcParams["position"]     ## center point
                        y_width = srcParams["y_width"]    ## dimensions of the sheet
                        x_width = srcParams["x_width"]    
                        
                        r0 , r1 = grid.r0, grid.r1
                        dr = np.array([grid.dx, grid.dy, grid.dz])
                        
                        j_r0 = np.array([j_r[0] - x_width/2.0, j_r[1] - y_width/2.0, j_r[2]])
                        j_r1 = np.array([j_r[0] + x_width/2.0, j_r[1] + y_width/2.0, j_r[2]])

                        overlap_r0 = np.maximum(r0, j_r0)
                        overlap_r1 = np.minimum(r1, j_r1)
                        
                        if np.all(overlap_r0 <= overlap_r1):
                            j_inds = ((j_r - r0)/dr).astype("int")
                            grid.AddSource(
                                {"type":"GaussianSheetSource_z",
                                "j_r": j_r,
                                "y_width": y_width,
                                "x_width": x_width,
                                "polarization": srcParams["polarization"], 
                                "amplitude": srcParams["amplitude"], 
                                "t_center": srcParams["t_center"], 
                                "t_decay": srcParams["t_decay"],
                                "modulationFrequency":srcParams["modulationFrequency"], 
                                "modulationPhase": srcParams["modulationPhase"], 
                                "timeOffsetFraction": srcParams["timeOffsetFraction"]
                                })
            
            elif srcParams["type"] == "PureScatteredRectPlaneWave":
                geomName = srcParams["geometryName"]
                for grid_dic in self.grids:
                    for grid in grid_dic.values():
                        
                        if geomName in self.gridsIntersectingGeometries[grid.name]:
                            boundingBox = self.gridsIntersectingGeometries[grid.name][geomName]

                            grid.AddSource(
                                {"type":"PureScatteredRectPlaneWave",
                                "boundingBox":boundingBox,
                                "polarization": srcParams["polarization"],
                                "amplitude": srcParams["amplitude"],
                                "propagationDirection": srcParams["propagationDirection"],
                                "propagationVelocity": srcParams["propagationVelocity"],
                                "t_center": srcParams["t_center"], 
                                "rectWidth": srcParams["rectWidth"], "rectEdgeWidth": srcParams["rectEdgeWidth"],
                                "modulationFrequency": srcParams["modulationFrequency"], 
                                "modulationPhase": srcParams["modulationPhase"],
                                "geometryName": geomName
                                })
            else:
                assert False            
                
    def SetupViews(self):
        n_levels = len(self.grids)
        for viewParams in self.views:
            if viewParams["type"] == "partial":
                for ind_lev in range(n_levels):
                    grid_dic = self.grids[ind_lev]
                    saveRate = 2**(n_levels - ind_lev - 1)
                    bufferSize = int(self.maxBufferSize / (2*ind_lev + 1))
                    if bufferSize < self.minBufferSize:
                        bufferSize = self.minBufferSize
                    for grid in grid_dic.values():
                        plane = viewParams["plane"]
                        at = viewParams["at"]
                        r0 , r1 = grid.r0, grid.r1
                        plane_dir = {"x":0, "y":1, "z":2}
                        
                        if r0[plane_dir[plane]] <= at < r1[plane_dir[plane]]:
                            fileName = self.viewFolder + grid.name + "_" + \
                                    viewParams["arrayName"] + "_" + viewParams["direction"] + \
                                    "_@" + plane + "=" + str(at)
                            grid.AddView(
                                {"type":"partial", 
                                "plane":plane, 
                                "at":at, 
                                "direction": viewParams["direction"], 
                                "arrayName": viewParams["arrayName"],
                                "fileName": fileName, 
                                "saveRate": saveRate,
                                "bufferSize": bufferSize
                                })
            elif viewParams["type"] == "entire":
                for ind_lev in range(n_levels):
                    grid_dic = self.grids[ind_lev]
                    saveRate = 2**(n_levels - ind_lev - 1)
                    for grid in grid_dic.values():                        
                        fileName = self.viewFolder + grid.name + "_" + \
                                viewParams["arrayName"] + "_" + viewParams["direction"] 
                        grid.AddView(
                            {"type":"entire", 
                            "direction": viewParams["direction"], 
                            "arrayName": viewParams["arrayName"],
                            "fileName": fileName, 
                            "saveRate": saveRate
                            })
            else: 
                assert False
        

    def SetupCenterGrid(self):
        box = self.boxes[0]
        r0, r1 = box["r0"], box["r1"]
        dr = box["dr"]
        dt = box["dt"]
        n = box["n"]
        grid_m = GridBlock(name = "grid_m0", blockLevel = 0, blockPosition = "c")
        grid_m.SetCorners(r0, r1)
        grid_m.SetNumOfCells(n[0], n[1], n[2])
        grid_m.SetCellDimentions(dr[0], dr[1], dr[2])
        grid_m.SetTimeStep(dt)
        
        self.grids = [{"c": grid_m}]
        self.gridCollectionner.AddGrid(grid_m)
        assert len(self.gridCollectionner.grids) == 1
        
    def SetupFirstLayerGrid(self):
        if len(self.boxes) <= 1:
            return

        box_0 = self.boxes[0]
        box_1 = self.boxes[1]

        r0_lev0, r1_lev0 = box_0["r0"], box_0["r1"]
        r0_lev1, r1_lev1 = box_1["r0"], box_1["r1"]
        dr_lev1 = box_1["dr"]
        dt_lev1 = box_1["dt"]
        
        n_lev1 = box_0["n"]/2
        dn_lev1 = self.layersNumOfCells[1]
        
        grid_r = GridBlock(name = "grid_r1", blockLevel = 1, blockPosition = "r")
        grid_r.SetCorners([r0_lev0[0], r0_lev0[1], r1_lev0[2]], [r1_lev0[0], r1_lev0[1], r1_lev1[2]])
        grid_r.SetNumOfCells(n_lev1[0], n_lev1[1], dn_lev1[2])
        grid_r.SetCellDimentions(dr_lev1[0], dr_lev1[1], dr_lev1[2])
        grid_r.SetTimeStep(dt_lev1)

        grid_l = GridBlock(name = "grid_l1", blockLevel = 1, blockPosition = "l")
        grid_l.SetCorners([r0_lev0[0], r0_lev0[1], r0_lev1[2]], [r1_lev0[0], r1_lev0[1], r0_lev0[2]])
        grid_l.SetNumOfCells(n_lev1[0], n_lev1[1], dn_lev1[2])
        grid_l.SetCellDimentions(dr_lev1[0], dr_lev1[1], dr_lev1[2])
        grid_l.SetTimeStep(dt_lev1)

        grid_u = GridBlock(name = "grid_u1", blockLevel = 1, blockPosition = "u")
        grid_u.SetCorners([r0_lev0[0], r1_lev0[1], r0_lev1[2]], [r1_lev0[0], r1_lev1[1], r1_lev1[2]])
        grid_u.SetNumOfCells(n_lev1[0], dn_lev1[1], n_lev1[2] + 2*dn_lev1[2])
        grid_u.SetCellDimentions(dr_lev1[0], dr_lev1[1], dr_lev1[2])
        grid_u.SetTimeStep(dt_lev1)

        grid_d = GridBlock(name = "grid_d1", blockLevel = 1, blockPosition = "d")
        grid_d.SetCorners([r0_lev0[0], r0_lev1[1], r0_lev1[2]], [r1_lev0[0], r0_lev0[1], r1_lev1[2]])
        grid_d.SetNumOfCells(n_lev1[0], dn_lev1[1], n_lev1[2] + 2*dn_lev1[2])
        grid_d.SetCellDimentions(dr_lev1[0], dr_lev1[1], dr_lev1[2])
        grid_d.SetTimeStep(dt_lev1)

        grid_f = GridBlock(name = "grid_f1", blockLevel = 1, blockPosition = "f")
        grid_f.SetCorners([r1_lev0[0], r0_lev1[1], r0_lev1[2]], [r1_lev1[0], r1_lev1[1], r1_lev1[2]])
        grid_f.SetNumOfCells(dn_lev1[0], n_lev1[1] + 2*dn_lev1[1], n_lev1[2] + 2*dn_lev1[2])
        grid_f.SetCellDimentions(dr_lev1[0], dr_lev1[1], dr_lev1[2])
        grid_f.SetTimeStep(dt_lev1)

        grid_b = GridBlock(name = "grid_b1", blockLevel = 1, blockPosition = "b")
        grid_b.SetCorners([r0_lev1[0], r0_lev1[1], r0_lev1[2]], [r0_lev0[0], r1_lev1[1], r1_lev1[2]])
        grid_b.SetNumOfCells(dn_lev1[0], n_lev1[1] + 2*dn_lev1[1], n_lev1[2] + 2*dn_lev1[2])
        grid_b.SetCellDimentions(dr_lev1[0], dr_lev1[1], dr_lev1[2])
        grid_b.SetTimeStep(dt_lev1)

        ## connections   
        grid_m = self.grids[0]['c']
         
        grid_m.AddConnection("r", grid_r)
        grid_r.AddConnection("l", grid_m)

        grid_m.AddConnection("l", grid_l)
        grid_l.AddConnection("r", grid_m)

        grid_m.AddConnection("u", grid_u)
        grid_r.AddConnection("u", grid_u)
        grid_l.AddConnection("u", grid_u)
        grid_u.AddConnection("dc", grid_m)
        grid_u.AddConnection("dr", grid_r)
        grid_u.AddConnection("dl", grid_l)

        grid_m.AddConnection("d", grid_d)
        grid_r.AddConnection("d", grid_d)
        grid_l.AddConnection("d", grid_d)
        grid_d.AddConnection("uc", grid_m)
        grid_d.AddConnection("ur", grid_r)
        grid_d.AddConnection("ul", grid_l)

        grid_m.AddConnection("f", grid_f)
        grid_r.AddConnection("f", grid_f)
        grid_l.AddConnection("f", grid_f)
        grid_u.AddConnection("f", grid_f)
        grid_d.AddConnection("f", grid_f)
        grid_f.AddConnection("bc", grid_m)
        grid_f.AddConnection("br", grid_r)
        grid_f.AddConnection("bl", grid_l)
        grid_f.AddConnection("bu", grid_u)
        grid_f.AddConnection("bd", grid_d)

        grid_m.AddConnection("b", grid_b)
        grid_r.AddConnection("b", grid_b)
        grid_l.AddConnection("b", grid_b)
        grid_u.AddConnection("b", grid_b)
        grid_d.AddConnection("b", grid_b)
        grid_b.AddConnection("fc", grid_m)
        grid_b.AddConnection("fr", grid_r)
        grid_b.AddConnection("fl", grid_l)
        grid_b.AddConnection("fu", grid_u)
        grid_b.AddConnection("fd", grid_d)

        self.grids.append({"r": grid_r, "l": grid_l, "u": grid_u, "d": grid_d, "f": grid_f, "b": grid_b})
        assert len(self.grids) == 2

        self.gridCollectionner.AddGrid(grid_r)
        self.gridCollectionner.AddGrid(grid_l)
        self.gridCollectionner.AddGrid(grid_u)
        self.gridCollectionner.AddGrid(grid_d)
        self.gridCollectionner.AddGrid(grid_f)
        self.gridCollectionner.AddGrid(grid_b)
        assert len(self.gridCollectionner.grids) == 7
        

    def SetupSecondLayerGrid(self, layerIndex):
        if len(self.boxes) <= 2:
            return
        assert layerIndex >= 2

        box_1 = self.boxes[layerIndex - 1]
        box_2 = self.boxes[layerIndex]

        r0_lev1, r1_lev1 = box_1["r0"], box_1["r1"]
        r0_lev2, r1_lev2 = box_2["r0"], box_2["r1"]
        dr_lev2 = box_2["dr"]
        dt_lev2 = box_2["dt"]
        
        n_lev2 = box_1["n"]/2
        dn_lev2 = self.layersNumOfCells[layerIndex]

        grid_rr = GridBlock(name = "grid_r{}".format(layerIndex), blockLevel = layerIndex, blockPosition = "r")
        grid_rr.SetCorners([r0_lev1[0], r0_lev1[1], r1_lev1[2]], [r1_lev1[0], r1_lev1[1], r1_lev2[2]])
        grid_rr.SetNumOfCells(n_lev2[0], n_lev2[1], dn_lev2[2])
        grid_rr.SetCellDimentions(dr_lev2[0], dr_lev2[1], dr_lev2[2])
        grid_rr.SetTimeStep(dt_lev2)

        grid_ll = GridBlock(name = "grid_l{}".format(layerIndex), blockLevel = layerIndex, blockPosition = "l")
        grid_ll.SetCorners([r0_lev1[0], r0_lev1[1], r0_lev2[2]], [r1_lev1[0], r1_lev1[1], r0_lev1[2]])
        grid_ll.SetNumOfCells(n_lev2[0], n_lev2[1], dn_lev2[2])
        grid_ll.SetCellDimentions(dr_lev2[0], dr_lev2[1], dr_lev2[2])
        grid_ll.SetTimeStep(dt_lev2)

        grid_uu = GridBlock(name = "grid_u{}".format(layerIndex), blockLevel = layerIndex, blockPosition = "u")
        grid_uu.SetCorners([r0_lev1[0], r1_lev1[1], r0_lev2[2]], [r1_lev1[0], r1_lev2[1], r1_lev2[2]])
        grid_uu.SetNumOfCells(n_lev2[0], dn_lev2[1], n_lev2[2] + 2*dn_lev2[2])
        grid_uu.SetCellDimentions(dr_lev2[0], dr_lev2[1], dr_lev2[2])
        grid_uu.SetTimeStep(dt_lev2)

        grid_dd = GridBlock(name = "grid_d{}".format(layerIndex), blockLevel = layerIndex, blockPosition = "d")
        grid_dd.SetCorners([r0_lev1[0], r0_lev2[1], r0_lev2[2]], [r1_lev1[0], r0_lev1[1], r1_lev2[2]])
        grid_dd.SetNumOfCells(n_lev2[0], dn_lev2[1], n_lev2[2] + 2*dn_lev2[2])
        grid_dd.SetCellDimentions(dr_lev2[0], dr_lev2[1], dr_lev2[2])
        grid_dd.SetTimeStep(dt_lev2)

        grid_ff = GridBlock(name = "grid_f{}".format(layerIndex), blockLevel = layerIndex, blockPosition = "f")
        grid_ff.SetCorners([r1_lev1[0], r0_lev2[1], r0_lev2[2]], [r1_lev2[0], r1_lev2[1], r1_lev2[2]])
        grid_ff.SetNumOfCells(dn_lev2[0], n_lev2[1] + 2*dn_lev2[1], n_lev2[2] + 2*dn_lev2[2])
        grid_ff.SetCellDimentions(dr_lev2[0], dr_lev2[1], dr_lev2[2])
        grid_ff.SetTimeStep(dt_lev2)

        grid_bb = GridBlock(name = "grid_b{}".format(layerIndex), blockLevel = layerIndex, blockPosition = "b")
        grid_bb.SetCorners([r0_lev2[0], r0_lev2[1], r0_lev2[2]], [r0_lev1[0], r1_lev2[1], r1_lev2[2]])
        grid_bb.SetNumOfCells(dn_lev2[0], n_lev2[1] + 2*dn_lev2[1], n_lev2[2] + 2*dn_lev2[2])
        grid_bb.SetCellDimentions(dr_lev2[0], dr_lev2[1], dr_lev2[2])
        grid_bb.SetTimeStep(dt_lev2)    
            
        ## connections
        grids_1 = self.grids[layerIndex - 1]
        grid_r = grids_1['r']
        grid_l = grids_1['l']
        grid_u = grids_1['u']
        grid_d = grids_1['d']
        grid_f = grids_1['f']
        grid_b = grids_1['b']
        
        grid_r.AddConnection("r", grid_rr)
        grid_u.AddConnection("r", grid_rr)
        grid_d.AddConnection("r", grid_rr)
        grid_f.AddConnection("r", grid_rr)
        grid_b.AddConnection("r", grid_rr)
        grid_rr.AddConnection("lc", grid_r)
        grid_rr.AddConnection("lu", grid_u)
        grid_rr.AddConnection("ld", grid_d)
        grid_rr.AddConnection("lf", grid_f)
        grid_rr.AddConnection("lb", grid_b)

        grid_l.AddConnection("l", grid_ll)
        grid_u.AddConnection("l", grid_ll)
        grid_d.AddConnection("l", grid_ll)
        grid_f.AddConnection("l", grid_ll)
        grid_b.AddConnection("l", grid_ll)
        grid_ll.AddConnection("rc", grid_l)
        grid_ll.AddConnection("ru", grid_u)
        grid_ll.AddConnection("rd", grid_d)
        grid_ll.AddConnection("rf", grid_f)
        grid_ll.AddConnection("rb", grid_b)

        grid_rr.AddConnection("u", grid_uu)
        grid_ll.AddConnection("u", grid_uu)
        grid_u.AddConnection("u", grid_uu)
        grid_f.AddConnection("u", grid_uu)
        grid_b.AddConnection("u", grid_uu)
        grid_uu.AddConnection("dc", grid_u)
        grid_uu.AddConnection("dr", grid_rr)
        grid_uu.AddConnection("dl", grid_ll)
        grid_uu.AddConnection("df", grid_f)
        grid_uu.AddConnection("db", grid_b)

        grid_rr.AddConnection("d", grid_dd)
        grid_ll.AddConnection("d", grid_dd)
        grid_d.AddConnection("d", grid_dd)
        grid_f.AddConnection("d", grid_dd)
        grid_b.AddConnection("d", grid_dd)
        grid_dd.AddConnection("uc", grid_d)
        grid_dd.AddConnection("ur", grid_rr)
        grid_dd.AddConnection("ul", grid_ll)
        grid_dd.AddConnection("uf", grid_f)
        grid_dd.AddConnection("ub", grid_b)

        grid_f.AddConnection("f", grid_ff)
        grid_rr.AddConnection("f", grid_ff)
        grid_ll.AddConnection("f", grid_ff)
        grid_uu.AddConnection("f", grid_ff)
        grid_dd.AddConnection("f", grid_ff)
        grid_ff.AddConnection("bc", grid_f)
        grid_ff.AddConnection("br", grid_rr)
        grid_ff.AddConnection("bl", grid_ll)
        grid_ff.AddConnection("bu", grid_uu)
        grid_ff.AddConnection("bd", grid_dd)

        grid_b.AddConnection("b", grid_bb)
        grid_rr.AddConnection("b", grid_bb)
        grid_ll.AddConnection("b", grid_bb)
        grid_uu.AddConnection("b", grid_bb)
        grid_dd.AddConnection("b", grid_bb)
        grid_bb.AddConnection("fc", grid_b)
        grid_bb.AddConnection("fr", grid_rr)
        grid_bb.AddConnection("fl", grid_ll)
        grid_bb.AddConnection("fu", grid_uu)
        grid_bb.AddConnection("fd", grid_dd)
            
        self.grids.append({"r": grid_rr, "l": grid_ll, "u": grid_uu, "d": grid_dd, "f": grid_ff, "b": grid_bb})
        assert len(self.grids) == layerIndex + 1
        
        self.gridCollectionner.AddGrid(grid_rr)
        self.gridCollectionner.AddGrid(grid_ll)
        self.gridCollectionner.AddGrid(grid_uu)
        self.gridCollectionner.AddGrid(grid_dd)
        self.gridCollectionner.AddGrid(grid_ff)
        self.gridCollectionner.AddGrid(grid_bb)
        assert len(self.gridCollectionner.grids) == layerIndex*6 + 1
        
        
    def SetupPMLs(self):
        if len(self.pmls) == 0:
            return None

        blockLevel = len(self.boxes) - 1
        assert len(self.boxes) >= 2
        box = self.boxes[-1]
        r0, r1 = box["r0"], box["r1"]
        dr = box["dr"]
        dt = box["dt"]
        n = box["n"]
                
        grids = self.grids[-1]
        grid_r = grids['r']
        grid_l = grids['l']
        grid_u = grids['u']
        grid_d = grids['d']
        grid_f = grids['f']
        grid_b = grids['b']
        
            
        if "f" in self.pmls:
            pml_params = self.pmls["f"]
            n_pml_cells = pml_params["n"]
            Dx_pml = n_pml_cells*dr[0]
            
            pml_f = GridBlock(name = "pml_f_", blockLevel = blockLevel, blockPosition = "f", gridType = "pml")
            pml_f.SetCorners([r1[0], r0[1], r0[2]], [r1[0] + Dx_pml, r1[1], r1[2]])
            pml_f.SetNumOfCells(n_pml_cells, n[1], n[2])
            pml_f.SetCellDimentions(dr[0], dr[1], dr[2])
            pml_f.SetTimeStep(dt)

            grid_f.AddConnection("pml-f", pml_f)
            pml_f.AddConnection("b", grid_f)

            self.gridCollectionner.AddGrid(pml_f)

        if "b" in self.pmls:
            pml_params = self.pmls["b"]
            n_pml_cells = pml_params["n"]
            Dx_pml = n_pml_cells*dr[0]
            
            pml_b = GridBlock(name = "pml_b_", blockLevel = blockLevel, blockPosition = "b", gridType = "pml")
            pml_b.SetCorners([r0[0] - + Dx_pml, r0[1], r0[2]], [r0[0], r1[1], r1[2]])
            pml_b.SetNumOfCells(n_pml_cells, n[1], n[2])
            pml_b.SetCellDimentions(dr[0], dr[1], dr[2])
            pml_b.SetTimeStep(dt)

            grid_b.AddConnection("pml-b", pml_b)
            pml_b.AddConnection("f", grid_b)

            self.gridCollectionner.AddGrid(pml_b)

            
        if "r" in self.pmls:
            pml_params = self.pmls["r"]
            n_pml_cells = pml_params["n"]
            Dz_pml = n_pml_cells*dr[2]
            
            pml_r = GridBlock(name = "pml_r_", blockLevel = blockLevel, blockPosition = "r", gridType = "pml")
            pml_r.SetCorners([r0[0], r0[1], r1[2]], [r1[0], r1[1], r1[2] + Dz_pml])
            pml_r.SetNumOfCells(n[0], n[1], n_pml_cells)
            pml_r.SetCellDimentions(dr[0], dr[1], dr[2])
            pml_r.SetTimeStep(dt)
            
            grid_r.AddConnection("pml-r", pml_r)
            grid_u.AddConnection("pml-r", pml_r)
            grid_d.AddConnection("pml-r", pml_r)
            grid_f.AddConnection("pml-r", pml_r)
            grid_b.AddConnection("pml-r", pml_r)
            pml_r.AddConnection("lc", grid_r)
            pml_r.AddConnection("lu", grid_u)
            pml_r.AddConnection("ld", grid_d)
            pml_r.AddConnection("lf", grid_f)
            pml_r.AddConnection("lb", grid_b)
            
            self.gridCollectionner.AddGrid(pml_r)


        if "l" in self.pmls:
            pml_params = self.pmls["l"]
            n_pml_cells = pml_params["n"]
            Dz_pml = n_pml_cells*dr[2]
            
            pml_l = GridBlock(name = "pml_l_", blockLevel = blockLevel, blockPosition = "l", gridType = "pml")
            pml_l.SetCorners([r0[0], r0[1], r0[2] - Dz_pml], [r1[0], r1[1], r0[2]])
            pml_l.SetNumOfCells(n[0], n[1], n_pml_cells)
            pml_l.SetCellDimentions(dr[0], dr[1], dr[2])
            pml_l.SetTimeStep(dt)

            grid_l.AddConnection("pml-l", pml_l)
            grid_u.AddConnection("pml-l", pml_l)
            grid_d.AddConnection("pml-l", pml_l)
            grid_f.AddConnection("pml-l", pml_l)
            grid_b.AddConnection("pml-l", pml_l)
            pml_l.AddConnection("rc", grid_l)
            pml_l.AddConnection("ru", grid_u)
            pml_l.AddConnection("rd", grid_d)
            pml_l.AddConnection("rf", grid_f)
            pml_l.AddConnection("rb", grid_b)

            self.gridCollectionner.AddGrid(pml_l)

        
    def SetupGrids(self):
        n_layers = len(self.layersNumOfCells)
        assert len(self.boxes) == n_layers
        
        for i in range(n_layers):
            if i == 0:
                self.SetupCenterGrid()
            elif i == 1:
                self.SetupFirstLayerGrid()
            else:
                self.SetupSecondLayerGrid(i)
                
        self.SetupPMLs()
        

    def GetGridParamsDic(self):
        params_dic = {"grids":{}, "views":self.views, "geometries":self.geometries,
                      "materials":self.materials, "sources":self.sources,
                      "boxes":self.boxes}
        for lev in range(len(self.grids)):
            grid_dic = self.grids[lev]
            for grid in grid_dic.values():
                r0 , r1 = grid.r0, grid.r1
                dr = np.array([grid.dx, grid.dy, grid.dz])
                n = np.array([grid.nx, grid.ny, grid.nz])
                dt = grid.dt
                params_dic["grids"][grid.name] = {"r0":r0, "r1":r1, "dr":dr, "dt":dt, "n":n, "level":lev}
        #print("\n\n", params_dic)
        return params_dic
        
        
    def SetupCollectionAndRun(self, t_max, filename, paramfileName):
        dt_coarse = self.boxes[-1]["dt"]
        nt = int(t_max/dt_coarse)
        self.gridCollectionner.SetNumOfCoarseTimeSteps(nt)
        print("Num of time steps: ", nt)
        
        self.SetupGrids()
        self.SetupGeometries()
        self.SetupMaterials()
        self.SetupSources()
        self.SetupViews()
        
        gridCollection = self.gridCollectionner.GenerateGridCollection()
        outfile = open(filename, "w")
        json.dump(gridCollection, outfile, indent=4)
        outfile.close()
                
        params_dic = self.GetGridParamsDic()
        paramfile = open(paramfileName, "wb")
        pickle.dump(params_dic, paramfile)
        paramfile.close()
        
            
