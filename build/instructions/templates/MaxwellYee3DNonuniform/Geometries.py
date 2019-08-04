
__all__ = ["Hyperboloid", "Cone", "Cylinder"]


import numpy as np


class Hyperboloid:
    def __init__(self, coneAngle, apexRadius, height, apexPosition):
        self.coneAngle = coneAngle
        self.apexRadius = apexRadius
        self.height = height
        self.apexPosition = apexPosition
        a, b = self.GetCanonicalScaleFators()
        self.coneTipPosition = np.array([apexPosition[0], apexPosition[1]+a, apexPosition[2]])
        
    def GetCanonicalScaleFators(self):
        b_a = np.arctan(self.coneAngle/2.0)
        b = self.apexRadius / b_a
        a = b / b_a
        return a, b
        
    def GetBoundingBox(self, y0, y1):
        a, b = self.GetCanonicalScaleFators()
        x_coneTip, y_coneTip, z_coneTip = self.coneTipPosition
        #print(a, b, self.coneTipPosition)
        y0_rel = y0 - y_coneTip
        y1_rel = y1 - y_coneTip
        y_min_rel = min(y0_rel, y1_rel)
        if y_min_rel < -a - self.height:
            y_min_rel = -a - self.height
        if y_min_rel <= -a:
            #x_bb_rel = abs(y_min_rel)*np.tan(self.coneAngle/2.0)  ## bb: boundinx box
            x_bb_rel = b/a*np.sqrt(y_min_rel**2 - a**2)  ## bb: boundinx box
            z_bb_rel = x_bb_rel
            y_bb_top = max(y0, y1)
            if y_bb_top > self.apexPosition[1]:
                y_bb_top = self.apexPosition[1]
            r0_bb = np.array([x_coneTip - x_bb_rel, min(y0, y1), z_coneTip - z_bb_rel])
            r1_bb = np.array([x_coneTip + x_bb_rel, y_bb_top, z_coneTip + z_bb_rel])
            return [r0_bb, r1_bb]
        else:
            return None

class Cone:
    def __init__(self, coneAngle, apexRadius, height, apexPosition):
        self.coneAngle = coneAngle
        self.apexRadius = apexRadius
        self.height = height
        self.apexToBaseDistance = height
        self.apexPosition = apexPosition
        tip_to_roundedtop = self.GetTipToApexDistance()[0]
        self.coneTipPosition = np.array([apexPosition[0], apexPosition[1]+tip_to_roundedtop, apexPosition[2]])
        
    def GetTipToApexDistance(self):
        tip_to_flattop = self.apexRadius * np.cos(self.coneAngle/2.0) / np.tan(self.coneAngle/2.0)
        tip_to_circleCenter = tip_to_flattop + self.apexRadius*np.sin(self.coneAngle/2.0)    # flattop : rounded cap removed
        tip_to_roundedtop = tip_to_circleCenter - self.apexRadius;

        sharpConeHeight = tip_to_roundedtop + self.apexToBaseDistance
        return tip_to_roundedtop, tip_to_flattop, tip_to_circleCenter
        
    def GetBoundingBox(self, y0, y1):
        tip_to_roundedtop, tip_to_flattop, tip_to_circleCenter = self.GetTipToApexDistance()
        x_coneTip, y_coneTip, z_coneTip = self.coneTipPosition

        y0_rel = y0 - y_coneTip
        y1_rel = y1 - y_coneTip
        y_min_rel = min(y0_rel, y1_rel)
        if y_min_rel < -tip_to_roundedtop - self.height:
            y_min_rel = -tip_to_roundedtop - self.height
        #print("y_min_rel:", y_min_rel, "y_min:", y_min_rel+tip_to_roundedtop, " self.height:", self.height)
        if y_min_rel <= -tip_to_roundedtop:
            x_bb_rel = None
            if y_min_rel > -tip_to_flattop:
                x_bb_rel_sq = self.apexRadius**2 - (tip_to_circleCenter + y_min_rel)**2  ## bb: boundinx box
                assert x_bb_rel_sq >= 0.0
                x_bb_rel = np.sqrt(x_bb_rel_sq)
            else:
                x_bb_rel = np.abs(y_min_rel)*np.tan(self.coneAngle/2)  ## bb: boundinx box
            z_bb_rel = x_bb_rel
            
            y_bb_top = max(y0, y1)
            if y_bb_top > self.apexPosition[1]:
                y_bb_top = self.apexPosition[1]

            y_bb_bot = min(y0, y1)
            if y_bb_bot < self.apexPosition[1] - self.height:
                y_bb_bot = self.apexPosition[1] - self.height
                
            if y_bb_bot > y_bb_top:
                return None

            r0_bb = np.array([x_coneTip - x_bb_rel, y_bb_bot, z_coneTip - z_bb_rel])
            r1_bb = np.array([x_coneTip + x_bb_rel, y_bb_top, z_coneTip + z_bb_rel])
            
            #print(y0, y1, x_bb_rel)
            return [r0_bb, r1_bb]
        else:
            return None
            

class Cylinder:
    def __init__(self, radius, height, topCenter):
        self.radius = radius
        self.height = height
        self.topCenter = topCenter
        
        
    def GetBoundingBox(self, y0, y1):
        y_top = self.topCenter[1]
        y_bot = y_top - self.height
        if y_top > max(y0, y1):
            y_top = max(y0, y1)
        if y_bot < min(y0, y1):
            y_bot = min(y0, y1)
        
        x_tc, z_tc = self.topCenter[0], self.topCenter[2]
        
        if y_bot < y_top:
            r0_bb = np.array([x_tc - self.radius, y_bot, z_tc - self.radius ])
            r1_bb = np.array([x_tc + self.radius, y_top, z_tc + self.radius ])
        
            return [r0_bb, r1_bb]
            
        else:
            return None
        

        
        
        
        
