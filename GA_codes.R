#-----------------------------------------------------------------
#                           GA - genalg
#-----------------------------------------------------------------

fitness_rbga <- function(idx){
  genes = which(idx == 1)
  dm_fit = new("DistMap", raw.data=raw.data,data=norm.data,
               insitu.matrix=as.matrix(insitu.matrix[,genes]),
               geometry=geometry)
  dm_fit <- binarizeSingleCellData(dm_fit, seq(0.15, 0.5, 0.01))
  dm_fit <- mapCells(dm_fit)
  dm_fit_mcc = dm_fit@mcc.scores
  eval = mean((dm_mcc - dm_fit_mcc)**2, na.rm = TRUE)
  #return(eval)
  if(!is.nan(eval)){return(1*eval)}
  else{return(1)}
}
library(genalg)
GAmodel_rbga1 <- rbga.bin(size = 84, popSize = 500 , suggestions = NULL,
                    iters = 100, elitism = NA, evalFunc = fitness_rbga,
                    verbose = TRUE,showSettings=TRUE,zeroToOneRatio = 4)

#-----------------------------------------------------------------
#                             GA - GA package
#-----------------------------------------------------------------

fitness_ga <- function(idx){
  genes = which(idx == 1)
  dm_fit = new("DistMap", raw.data=raw.data,data=norm.data,
               insitu.matrix=as.matrix(insitu.matrix[,genes]),
               geometry=geometry)
  dm_fit <- binarizeSingleCellData(dm_fit, seq(0.15, 0.5, 0.01))
  dm_fit <- mapCells(dm_fit)
  dm_fit_mcc = dm_fit@mcc.scores
  eval = mean((dm_mcc - dm_fit_mcc)**2, na.rm = TRUE)
  #return(eval)
  if(!is.nan(eval)){return(-1*eval)}
  else{return(-1)}
}

library(GA)

GA_model1 <- ga(type = "binary",
                fitness = fitness_ga ,
                nBits = 84,
                popSize = 100,
                pcrossover = 0.9,
                pmutation = 0.1,
                maxiter = 10,
                run = 10,
                keepBest = TRUE,
                parallel = TRUE
)

GA_model2 <- ga(type = "permutation",
                lower = 1,
                upper = 22,
                seed = 12345,
                fitness = fitness_func ,
                nBits = 84,
                popSize = 100,
                pcrossover = 0.8,
                pmutation = 0.1,
                maxiter = 100,
                run = 10,
                keepBest = TRUE,
                parallel = TRUE
)

# controlx = list(fnscale=1)
# GAmodel <- ga(type = "binary",fitness = fitness_func,pmutation = 0.05,
#               popSize = 50,maxiter=10,run = 3,parallel=TRUE, elitism = 0.1,
#               lower = 1,upper = 84,nBits = 84,
#               optimArgs = list(method = "L-BFGS-B",
#                                poptim = 0.05,pressel = 0.5,control = controlx),
#               keepBest = T)



