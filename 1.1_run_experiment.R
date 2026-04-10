# #--read the arguments --- #
# args = commandArgs(trailingOnly = TRUE)
# p = as.numeric(args[1])
# n = as.numeric(args[2])
# graph = args[3]
# iter_rep = as.numeric(args[4])
p = 300
n_list = c(100)
prob_list = c(0.05)
graph_list = c("cluster")
iter_rep = 30

#--download the necessary libraries  --
library(doParallel)

#--load additional files --
# dir = "/home/tn325/HCP/sbggm_final/"
dir = "~/Documents/Cornell/Projects/GSEMMS/Code/Experiments/"

out_dir = paste0(dir, "out/")

source(paste0(dir, "1.0_experiment_helper.R"))


#--set parameters to simulate data --#
# prob = 0.2
type = "Gaussian"
size = NULL
vis = FALSE
seed = 1

#--set parameters to solve data --#
iter = 5000
burnin= 2000
save = TRUE
verbose = FALSE
g.start = "empty"
cores = 1

jump = 1 
var1 = 0.02 #ss parameter
var2 = 2 #ss parameter
lambda = 2 #ss parameter

g.prior = 0.2

#-- run experiments in parallel --#
ncores <- 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

outlist = foreach(i = 1:iter_rep) %:% 
  foreach(j = 1:length(n_list)) %dopar% {
    
    # for (i in 1:iter_rep) {
    #   print(i)
    #   for (j in 1:length(n_list)) {
    library( BDgraph )
    library (pROC)
    library (ssgraph)
    library(glasso)
    library(huge)
    library(parallel)
    library(broom)
    library(GEMMS)
    
    n <- n_list[j]
    for (prob in prob_list) {
      for (graph in graph_list) {
        result = run_experiments( p = p, n = n, prob = prob, graph = graph, type = type, vis = vis, 
                                  jump = jump, iter = iter, burnin = burnin, save = save, cores = cores,
                                  verbose = verbose,g.start =g.start,var1=var1,var2=var2,lambda=lambda,
                                  g.prior = g.prior, seed=i+2023) 
        
        #--print data to a Rdata file
        filename = paste0("result_p",p,"_n",n,"_",graph,"_prob", prob,"_rep",i,".Rdata")
        save( result, file = paste0(out_dir, filename ))
      }

      
    }
  }



# Stop the cluster
stopCluster(cl)



