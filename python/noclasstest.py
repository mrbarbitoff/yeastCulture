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
            
def upProt(cell, time):
    cell[6] += math.floor(synthesis[0]*(time - cell[1]))
    cell[7] += math.floor(synthesis[1]*(time - cell[1]))
    return cell

def divide(cell):
    global dc
    dc += 1
    probs = np.random.random(len(cell[5]))
    remaining = []
    going = []
    for i in range(len(cell[5])):
        if cell[5][i] < threshold and probs[i] >= 0.6:
            going.append(cell[5][i])
        else:
            remaining.append(cell[5][i])  
    mother = [cell[0], cell[1] + 1, cell[2] + divtime(True), True, len(remaining) != 0, remaining, 
              math.floor(cell[6]*0.6), math.floor(cell[7]*0.6)]
    daughter = [dc, cell[1] + 1, cell[2] + divtime(False), False, len(going) != 0, going, 
              math.ceil(cell[6]*0.4), math.ceil(cell[7]*0.4)]          
    return mother, daughter          
              
stand_cell = [1, 0, divtime(False), False, True, [20, 20, 20, 40, 40, 40, 40, 50], 45326, 3855]

def simCult(termination = Inf):

    culture = [stand_cell]
    psiplus = 0
    ratios = []

    for cell in culture:
        
        while True:
            if cell[4] == False: cell = upProt(cell, min(1000, cell[2]))
            elif cell[1] >= cell[2]: 
                cell, daughter = divide(cell)
                culture.append(daughter)
            elif cell[1] >= 1000: break
            else: 
            
                synt = sum(synthesis)
                aggregation = cell[6]*len(cell[5])*conv
                Z = sum(cell[5]) - len(cell[5])
                fragmentation = frag*cell[7]*Z/(cell[7]/2 + Z)
                total = synt + aggregation + fragmentation
                
                interval = np.random.exponential(1/total)
                reaction = float(np.random.random(1)*total)
                
                if reaction <= synthesis[0]:
                    cell[6] += 1
                elif reaction <= sum(synthesis):
                    cell[7] += 1
                elif reaction <= (synt + fragmentation):
                    dagg = math.ceil((reaction - synt)*cell[7]/(2*(frag*cell[7] - (reaction - synt))))
                    rsum = np.cumsum([alen - 1 for alen in cell[5]])
                    agg = bisect.bisect_left(rsum, dagg)
                    rplace = rsum[agg] - dagg + 1
                    remlen = cell[5][agg] - rplace
                    cell[5].remove(cell[5][agg])
                    for a in [rplace, remlen]:
                        if a >= 6:
                            cell[5].append(a)
                        else:
                            cell[6] += a
                    if len(cell[5]) == 0:
                        cell[4] = False
                        continue        
                else:
                    agg = np.random.randint(len(cell[5]))
                    cell[5][agg] += 1
                    cell[6] -= 1
                            
                cell[1] += interval  
                
        ratios.append((cell[6]/(cell[6] + sum(cell[5]))))
        if cell[4]:
            psiplus += 1
        if cell[0] >= termination: break    

    PR = psiplus/len(culture)
    SR = np.mean(ratios)
    
    return PR, SR           
