import bisect
from datetime import datetime
import numpy as np
import math
import matplotlib

Inf = float('inf')

#b = range()
#g = range()
synthesis = [700, 50]
conv = 0.0000065
frag = 0.006

rh = 31.03 
lambd_d = 0.21*60                                                                     
lambd = 1.16*60                                                                       
shape_d = rh*lambd_d                                                                 
shape = rh*lambd

dc = 1
threshold = Inf

def divtime(mother):
        if mother:
            return np.random.gamma(shape, 1/rh)
        else:
            return (np.random.gamma(shape, 1/rh) + np.random.gamma(shape_d, 1/rh))
            
class yeastcell:
    def __init__(self, no, time, div, mother, status, aggregates, Sup35, Hsp104):
        self.number = no
        self.time = time
        self.division = div
        self.mother = mother
        self.psi = status
        self.aggregates = aggregates
        self.Sup35 = Sup35
        self.Hsp104 = Hsp104
            
    def divide(self):
        global dc
        dc += 1
        self.time += 1
        self.mother = True
        self.Sup35 = math.floor(0.6*self.Sup35)
        self.Hsp104 = math.floor(0.6*self.Hsp104)
        probs = np.random.random(len(self.aggregates))
        remaining = []
        going = []
        for i in range(len(self.aggregates)):
            if self.aggregates[i] < threshold and probs[i] >= 0.6:
                going.append(self.aggregates[i])
            else:
                remaining.append(self.aggregates[i])        
        self.aggregates = remaining        
        return yeastcell(dc, self.time, self.division + divtime(False), False, len(going) != 0, 
               going, math.ceil(self.Sup35*0.4), math.ceil(self.Hsp104 * 0.4))
    
    def upProt(self, time):
        self.Sup35 = self.Sup35 + synthesis[0]*(time - self.time)
        self.Hsp104 = self.Hsp104 + synthesis[1]*(time - self.time)
        
    def copy(self):
        return yeastcell(self.number, self.time, self.division, self.mother, self.psi,
                         self.aggregates, self.Sup35, self.Hsp104)
               
stand_cell = yeastcell(1, 0, 0, False, True, [20, 20, 20, 40, 40, 40, 40, 50], 45326, 3855)

def simCult(termination = Inf):

    culture = [stand_cell]
    psiplus = 0
    ratios = []

    for cell in culture:
        
        while True:
            if cell.psi == False: cell.upProt(min(1000, cell.division))
            elif cell.time >= cell.division: 
                culture.append(cell.divide())
                cell.division += divtime(True)
            elif cell.time >= 1000: break
            else: 
            
                synt = sum(synthesis)
                aggregation = cell.Sup35*len(cell.aggregates)*conv
                Z = sum(cell.aggregates) - len(cell.aggregates)
                fragmentation = frag*cell.Hsp104*Z/(cell.Hsp104/2 + Z)
                total = synt + aggregation + fragmentation
                
                interval = np.random.exponential(1/total)
                reaction = float(np.random.random(1)*total)
                
                if reaction <= synthesis[0]:
                    cell.Sup35 += 1
                elif reaction <= sum(synthesis):
                    cell.Hsp104 += 1
                elif reaction <= (synt + fragmentation):
                    dagg = math.ceil((reaction - synt)*cell.Hsp104/(2*(frag*cell.Hsp104 - (reaction - synt))))
                    rsum = np.cumsum([alen - 1 for alen in cell.aggregates])
                    agg = bisect.bisect_left(rsum, dagg)
                    rplace = rsum[agg] - dagg + 1
                    remlen = cell.aggregates[agg] - rplace
                    cell.aggregates.remove(cell.aggregates[agg])
                    for a in [rplace, remlen]:
                        if a >= 6:
                            cell.aggregates.append(a)
                        else:
                            cell.Sup35 += a
                    if len(cell.aggregates) == 0:
                        cell.psi = False
                        continue        
                else:
                    agg = np.random.randint(len(cell.aggregates))
                    cell.aggregates[agg] += 1
                    cell.Sup35 -= 1
                            
                cell.time += interval  
                
        ratios.append((cell.Sup35/(cell.Sup35 + sum(cell.aggregates))))
        if cell.psi:
            psiplus += 1
        if cell.number >= termination: break    

    PR = psiplus/len(culture)
    SR = np.mean(ratios)
    
    return PR, SR           
