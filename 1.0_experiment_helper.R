library(SEMMS) # Bar(2015)
library(parallel)
library(BayesS5)
library(glasso)

all_cores <<- detectCores()



#### 3. Run experiments ####

run_experiments = function( p , n , prob , graph = "random", type = "Gaussian", 
                            vis = FALSE, jump = 1, iter = 10000, burnin = 7000, 
                            save = TRUE, cores = 1, verbose = TRUE,g.start="empty",var1=var1,
                            var2=var2,lambda=lambda, seed=1, g.prior=0.2)
{
  
  ##-----simulate data-------------##
  set.seed(10*seed)
  if (graph=="random") {
    data.sim = bdgraph.sim( p = p, n = n,
                            graph = graph, prob=prob,
                            type="Gaussian", vis=FALSE)
    true_g = as.matrix( data.sim $ G )
    true_K =  as.matrix( data.sim $ K )
    true_sigma =  as.matrix( data.sim $ sigma )
  } else if (graph=="cluster") {
    data.sim = bdgraph.sim( p = p, n = n,
                            graph = graph, prob=prob,
                            type="Gaussian", vis=FALSE, 
                            class=3)
    true_g = as.matrix( data.sim $ G )
    true_K =  as.matrix( data.sim $ K )
    true_sigma =  as.matrix( data.sim $ sigma )
  } else if (graph=="scale-free") {
    data.sim = bdgraph.sim( p = p, n = n,
                            graph = graph, vis=FALSE)
    true_g = as.matrix( data.sim $ G )
    true_K =  as.matrix( data.sim $ K )
    true_sigma =  as.matrix( data.sim $ sigma )
  } else if (graph=="band") {
    K=huge.generator(graph=graph,g=3)$omega
    data.sim = bdgraph.sim( p = p, n = n,
                            graph = "fixed",K=K, vis=FALSE)
    # data.sim = bdgraph.sim( p = p, n = n,
    #                         graph = graph, prob=prob,
    #                         type="Gaussian", vis=FALSE)
    # true_g = as.matrix( data.sim $ G )
    # true_K =  as.matrix( data.sim $ K )
    # true_sigma =  as.matrix( data.sim $ sigma )
    # data.sim = huge.generator(p = p, n = n,
    #                           graph = graph, g = 5)
    # true_g = as.matrix( data.sim $ theta )
    # true_K =  as.matrix( data.sim $ omega )
    # true_sigma =  as.matrix( data.sim $ sigma )
  }
  
  
  # #----solve data using bdmcmc method with app ratio of norm constants---#
  # set.seed(seed+3)
  # t1_bd      = proc.time()
  # sample_bd  = bdgraph.mpl( data = data.sim, algorithm = "bdmcmc", iter=100000, burnin = burnin, jump = jump, 
  #                       save = save,cores=cores,
  #                       g.start=g.start,g.prior=g.prior) 
  # time_bd   = as.numeric( ( proc.time() - t1_bd )[ 3 ] )
  
  # all_weights_bd = round(sample_bd$all_weights,10) #obtain weights for AUC vs iteration graphh
  # all_graphs_bd = sample_bd$all_graphs #obtain graphs for AUC vs iteration graphh
  # sample_graphs_bd = sample_bd$sample_graphs #obtain graphs for AUC vs iteration graphh
  
  #----solve data using ss method --------#
  set.seed(seed+5)
  t1_ss      = proc.time()
  if (n <= p) {
    # ensure that precision matrix is positive definite
    shift <- min(eigen(cov(data.sim$data))$values)
    incov <- cov(data.sim$data) - diag(shift,nrow=p,ncol=p) + diag(10e-6,nrow=p,ncol=p)
    
    sample_ss  = ssgraph( data = data.sim, iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,
                          g.start=g.start,save=save,cores=cores,g.prior=g.prior,
                          sig.start=incov)
  } else {
    sample_ss  = ssgraph( data = data.sim, iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,
                          g.start=g.start,save=save,cores=cores,g.prior=g.prior,
                          sig.start=cov(data.sim$data))
  }
  time_ss    = as.numeric( ( proc.time() - t1_ss )[ 3 ] )
  
  all_weights_ss = round(sample_ss$all_weights,10) #obtain weights for AUC vs iteration graphh
  all_graphs_ss = sample_ss$all_graphs #obtain graphs for AUC vs iteration graphh
  sample_graphs_ss = sample_ss$sample_graphs #obtain graphs for AUC vs iteration graphh
  
  
  # #----solve data using semms method---#
  # t1_semms      = proc.time()
  # sample_semms <- run_nodewise_reg(data=data.sim$data, method="SEMMS",
  #                                  rule="or", cores=cores, nn=30)
  # time_semms   = as.numeric( ( proc.time() - t1_semms )[ 3 ] )
  
  # #----solve data using cv glasso method---#
  # t1_glasso      = proc.time()
  # out.glasso = huge(data.sim$data, method = "glasso")
  # sample_glasso <- huge.select(out.glasso, criterion="ric")$opt.icov
  # time_glasso   = as.numeric( ( proc.time() - t1_glasso )[ 3 ] )
  
  
  
  return(list( 
    # # bd mcmc
    # time_bd = time_bd, all_weights_bd = all_weights_bd,
    # all_graphs_bd = all_graphs_bd, sample_graphs_bd = sample_graphs_bd,
    # # K_bd = as.matrix(sample_bd $ K_hat),
    # ss mcmc
    time_ss = time_ss, all_weights_ss = all_weights_ss,
    all_graphs_ss = all_graphs_ss, sample_graphs_ss = sample_graphs_ss,
    K_ss = as.matrix(sample_ss $ K_hat),
    # # semms
    # K_semms = as.matrix(sample_semms $ K_hat),
    # plinks_semms = as.matrix(sample_semms $ p_links),
    # t_semms = as.matrix(sample_semms $ t_stat),
    # time_semms = time_semms,
    # # cv glasso
    # K_glasso = as.matrix(sample_glasso),
    # time_glasso = time_glasso,
    # true graph
    true_g = true_g,
    true_K =  true_K,
    true_sigma =  true_sigma,
    true_data = as.matrix( data.sim $ data ) )
  )
}
