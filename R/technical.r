# Reaction rates
synthesis = c(700, 50) 
synt = sum(synthesis)
b = seq(from=0.0000065, to=0.0003, length.out=16)                                                                            
g = seq(from=0.001, to=0.013, length.out=25)
conv = 0.0000065
frag = 0.006

# Division parameters 
rho = 31.03                                                                                                                                    
shape_d = rho*0.21*60                                                                 
shape = rho*1.16*60
div_count = 1
threshold = Inf

# Simulating time until next mother division
mTime = function() {
  return(rgamma(1, shape=shape, scale=1/rho))
}

# Simulating time until becoming a mother
dTime = function() {
  return(rgamma(1, shape=shape_d, scale=1/rho) + mTime())
}

# Typical cell for simulation
standard_cell = c(1, 0, dTime(), 1, 0, 0, 46680, 3326)

# Updating protein levels
upProt = function(cell, time){
    cell[7] = cell[7] + floor((time - cell[2])*synthesis[1])
    cell[8] = cell[8] + floor((time - cell[2])*synthesis[2])
    cell[2] = time
    return(cell)
}

# Writing short outputs to files
writeOut = function() {    
    write.table(prmatr, file='psi.txt', col.names=F, row.names=F, sep='\t')
    write.table(srmatr, file='sup.txt', col.names=F, row.names=F, sep='\t')
    write.table(almatr, file='agl.txt', col.names=F, row.names=F, sep='\t')
}

# Plotting a heat map for simulation output
library(colorRamps)
library(lattice)
jet = matlab.like2(100)

plotHeatMap = function(psi, sup){    
    #Rotating input matrices for proper imaging
    psi = t(apply(psi, 2, rev))
    sup = t(apply(sup, 2, rev))    
    #Constructing image for the first input matrix 
    psi = levelplot(psi, col.regions=jet, at=seq(0, 1, length.out=100), aspect='fill',
        colorkey=list(width=3, labels=list(cex=1.0)), 
        xlab=list(label='Conversion (ascending)', cex=1.2),
        ylab=list(label='Fragmentation (ascending)', cex=1.2), 
        main=list(label='Prion stability', cex=1.35, font=2), 
        xlim=c(0.5, nrow(psi) + 0.5), ylim=c(0.5, ncol(psi) + 0.5), scales=list(draw=F))    
    #Constructing image for the second input matrix    
    sup = levelplot(sup, col.regions=jet, at=seq(0, 1, length.out=100), aspect = 'fill',
        colorkey=list(width=3, labels=list(cex=1.0)), 
        xlab=list(label='Conversion (ascending)', cex=1.2),
        ylab=list(label='Fragmentation (ascending)', cex=1.2), 
        main=list(label='Soluble Sup35p fraction', cex=1.35, font=2), 
        xlim=c(0.5, nrow(sup) + 0.5), ylim=c(0.5, ncol(sup) + 0.5), scales=list(draw=F))    
    #Plotting images side-by-side to a PNG file
    png(file='Simulation.png', width=1200, height=500)
    print(psi, more=T, pos=c(0, 0, 0.48, 1))
    print(sup, pos=c(0.52, 0, 1, 1))
    dev.off()    
}

