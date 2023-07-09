# Simulate single culture
simCult = function(termination=Inf, cell=standard_cell, aggrs=c(20, 20, 20, 40, 40, 40, 40, 50)){ 
	# Initializing culture matrix and aggregate variable
    assign(paste0('agg', '1'), aggrs)
    culture = t(as.matrix(cell))
    k = 1
    proceed = T
	# Simulating
    while (T) {
        mycell = culture[k,]
        aggrs = get(paste0('agg', k))        
        while(mycell[2] <= 1000) {                                                       
            if(length(aggrs) == 0) {
                mycell = upProt(mycell, min(mycell[3], 1000))
            }
            if(mycell[2]>= mycell[3]) {                
				# Dividing the cell here
                div_count <<- div_count + 1
                probs = runif(length(aggrs))                 
                assign(paste0('agg', div_count), aggrs[probs < 0.4&aggrs <= threshold])
                aggrs = aggrs[probs > 0.4|aggrs > threshold]  
                dcell = c(div_count, mycell[2] + 1, mycell[3] + dTime(), as.numeric(length(agg_d) != 0), 0, 0, 
                          ceiling(mycell[7:8]*0.4))                          
                mycell = c(mycell[1], mycell[2] + 1, mycell[3] + mTime(), as.numeric(length(aggrs) != 0), 0, 0, 
                           floor(mycell[7:8]*0.6))     
                culture = rbind(culture, dcell, deparse.level=0)
            } else if(mycell[2] >= 1000) {
                break
            } else {              
				# Running Gillespie's algorithm
                aggregation = mycell[7]*length(aggrs)*conv
                Z = sum(aggrs) - length(aggrs)
                fragmentation = (Z*frag*mycell[8])/((mycell[8]/2) + Z)
                total = aggregation + fragmentation + synt                
                interval = rexp(1, total)
                reaction = runif(1)*total                
                if(reaction <= synthesis[1]) {mycell[7] = mycell[7] + 1}                    
                else if(reaction <= synt) {mycell[8] = mycell[8] + 1}                      
                else if(reaction <= synt + fragmentation) {
                    dagg = ceiling(((reaction - synt)*mycell[8]/2)/(mycell[8]*frag - (reaction - synt)))
                    adef = cumsum(aggrs - 1)
                    cagg = findInterval(dagg, adef, rightmost.closed=T) + 1                    
                    place = adef[cagg] - dagg + 1
                    remlen = aggrs[cagg] - place                        
                    # Sorting protein from the fragmented aggregate to soluble/insoluble pools    
                    tmp = c(place, remlen)
                    newsol = sum(tmp[tmp < 6])
                    newaggrs = tmp[tmp >= 6]
                    aggrs = c(aggrs[-cagg], newaggrs)
                    mycell[7] = mycell[7] + newsol
                } else {
                    cur_agg = sample(1:length(aggrs), 1)
                    aggrs[cur_agg] = aggrs[cur_agg] + 1
                    mycell[7] = mycell[7] - 1
                }
                mycell[2] = mycell[2] + interval
            }
        }
		# Calculating final parameters
        if(length(aggrs) >= 1) {mycell[4] = 1}
        mycell[5] = mycell[7]/(mycell[7] + sum(aggrs))
        mycell[6] = mean(aggrs)        
        culture[k,] = mycell            
        if(k < nrow(culture)&k < termination) {k = k + 1} else {break}
    }
    finale = sapply(culture, 2, mean)
    output = finale[4:6]    
    return(output)    
}

library(doParallel)

# Simulate multiple cultures
simPar = function(ncores=13, columns=length(b), aggregates=c(20, 20, 20, 40, 40, 40, 40, 50)){    
    # Output matrices
    prmatr = matrix(0, length(g), length(b))
    srmatr = prmatr
    almatr = prmatr    
    # Actual simulation here
    for(i in 1:columns) {
        conv = b[i]
        registerDoParallel(cores=ncores)
        stime = Sys.time()        
        column = foreach(j = 1:length(g), .export = ls(envir=.GlobalEnv)) %dopar% {
            frag = g[j]
            ratios = simCult(aggrs=aggregates)
        }        
        PRC = sapply(column, function(elt) elt[1])
        SRC = sapply(column, function(elt) elt[2])
        ALC = sapply(column, function(elt) elt[3])        
        prmatr[,i] = rev(PRC)
        srmatr[,i] = rev(SRC)
        almatr[,i] = rev(ALC)        
        duration = round(as.numeric(difftime(Sys.time(), stime, units="hours")), 1)
        message('I`ve finished simulating column numero', i, 'at', as.character(Sys.time()), 
            'MSK! \nIt took me', duration, 'hours to simulate it. \nRegards, \nSlowbro-code. (yup, I`ve got a lvl-up!)',
            '\n===================\n')
        writeOut()
    }    
    plotHeatMap(prmatr, srmatr)    
}
