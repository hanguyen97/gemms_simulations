
calc_AUC_ROC = function(predictor,response){
  
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #order the vectors so that the predictor is increasing
  predictor.order = order(predictor,decreasing=FALSE)
  predictor.sorted = predictor[predictor.order]
  response.sorted = response[predictor.order]
  
  #determine amount of zeroes and ones
  ones = sum(response)
  zeroes = length(response)-ones
  
  #if there are duplicates
  if (sum(duplicated(predictor.sorted))>0){
    #create a vector with one index for every group of duplicates
    dup_index = cumsum(duplicated(predictor.sorted)==0)
    
    #create a vector sum_vec that sums the true positives in each group of duplicates   
    df <- data.frame(duplicates=dup_index,response.sorted=response.sorted)
    sum_vec = aggregate(response.sorted ~ duplicates, data=df, sum)[,2]
    
    #create a vector that averages the maximum amount of false positives of the current group with the previous group
    fp = cumsum(response.sorted==0)
    df <- data.frame(duplicates=dup_index,fp=fp)
    max_vec = aggregate(fp ~ duplicates, data=df, max)[,2]
    top = c(0,max_vec)
    bottom = c(max_vec,0)
    average_vec = head((top+bottom)/2,-1)
    
    #AUC is the dot product of the two vectors divided by the normalizing constant
    AUC = (sum_vec%*%average_vec)/(ones*zeroes)
    
  }
  
  #if there are no duplicates
  if (sum(duplicated(predictor.sorted))==0){
    fp = cumsum(response.sorted==0)
    AUC = sum(fp * response.sorted)
    AUC = AUC/(zeroes*ones)  
    
  }
  
  return(AUC)
}

calc_AUC_ROC_calibrated = function(predictor,response,cut_AUC_calibrated = 200)
{
  output_tp_fp = BDgraph:::compute_tp_fp(actual = response, pred = predictor, cut = cut_AUC_calibrated, smooth = FALSE,calibrate=TRUE)
  fp = output_tp_fp$fp
  tp = output_tp_fp$tp
  diffs_x = fp[-length(fp)] - fp[-1]
  means_vert = (tp[-1] + tp[-length(tp)])/2
  auc_calibrated = sum(means_vert * diffs_x)
  return(auc_calibrated)
}

calc_AUC_PR = function(predictor,response){
  
  #when the predictor only contains zeroes or ones, then AUC_pr will return zero
  if ((sum(predictor)==0)|(sum(predictor)==length(predictor))){
    AUC_pr = 0
  } else {
    
    fg = predictor[response==1]
    bg = predictor[response==0]
    pr = pr.curve(scores.class0 = fg, scores.class1 = bg,curve=FALSE)
    
    AUC_pr = pr$auc.davis.goadrich
  }
  return(AUC_pr)
}

calc_CE_MSE = function(predictor,response){
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #determine length of vectors
  qp = length(response)
  
  #determine amount of zeroes and ones
  ones = sum(response)
  zeroes = length(response)-ones
  
  #calculate CE, MSE 
  CE = sum(abs(predictor-response))/qp
  MSE = sum((predictor-response)^2)/qp
  
  #calculate the weighed CE and weighed MSE
  ones_index = which(response==1)
  zeroes_index = which(response==0)
  predictor_ones =  predictor[ones_index]
  predictor_zeroes = predictor[zeroes_index]
  CE_ones = sum(1-predictor_ones)/ones
  CE_zeroes = sum(predictor_zeroes)/zeroes
  MSE_ones = sum((1-predictor_ones)^2)/ones
  MSE_zeroes = sum(predictor_zeroes^2)/zeroes
  CE_weighed = 0.5*CE_ones + 0.5*CE_zeroes
  MSE_weighed = 0.5*MSE_ones + 0.5*MSE_zeroes
  
  #return the metrics
  return(list(CE=CE,CE_weighed=CE_weighed,MSE=MSE,MSE_weighed=MSE_weighed))
}

calc_TP_FP = function(predictor,response,cutoff){
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #determine length of vectors
  qp = length(response)
  
  #determine amount of zeroes and ones
  ones = sum(response)
  zeroes = length(response)-ones
  
  #determine estimated graph from predictor and cutoff point
  positive_index = which(predictor>=cutoff)
  estimated_graph = c(rep(0,qp))
  estimated_graph[positive_index] = 1
  
  #obtain TP,FP,FN,TN
  tp = sum((response == 1) * (estimated_graph == 1))
  tn = sum((response == 0) * (estimated_graph == 0))
  fp = sum((response == 0) * (estimated_graph == 1))
  fn = sum((response == 1) * (estimated_graph == 0))
  
  #obtain the rate of TP,FP,FN,TN
  tpr = tp/ones
  tnr = tn/zeroes
  fpr = fp/zeroes
  fnr = fn/ones
  
  #return the metrics
  return(list(tp=tp,tn=tn,fp=fp,fn=fn,tpr=tpr,tnr=tnr,fpr=fpr,fnr=fnr)) 
  
}

calculate_metrics = function (all_graphs=NULL,sample_graphs=NULL,all_weights=NULL,response, 
                              thin = 100, plinks_diff = 1000, cut_AUC_calibrated,verbose = TRUE,cutoff=0.5,
                              K_hat=NULL, true_K=NULL, type="mcmc", p_links=NULL, threshold_stat=NULL) 
{
  if (is.null(response))
    stop("please include the argument response")
  # Stein's loss
  # SteinL <- Stein(S=true_K, H=K_hat)
  SteinL <- NULL
  # F norm
  Fnorm <- norm(true_K - K_hat, type = c("F"))
  
  if (type=="mcmc") {
    #check if required information is given
    if (is.null(all_graphs))
      stop("please include the argument all_graphs")
    if (is.null(sample_graphs))
      stop("please include the argument sample_graphs")
    if (is.null(all_weights))
      stop("please include the argument all_weights")
    
    #check if the thin value is in the right format
    if (is.null(thin)) 
      stop("'thin' must be a number")
    if (!is.numeric(thin)) 
      stop("'thin' must be a number")
    if (is.matrix(thin)) 
      stop("'thin' must be a number")
    
    #obtain amount of entries
    qp = length(response)
    
    #create output -- iteration vectors
    iter = length(all_graphs)
    iter_vec = c(1:iter)
    iter_vec_thin = c(thin * (1:floor(iter/thin)))
    
    #create output -- plinks vectors
    result_plinks = matrix(0,qp,iter)
    totaltime_per_edge = c(rep(0,qp))
    
    #create thinned output vectors
    AUC_ROC_vec = c(rep(0,length(iter_vec_thin)))
    AUC_PR_vec = c(rep(0,length(iter_vec_thin)))
    AUC_ROC_calibrated_vec = c(rep(0,length(iter_vec_thin)))
    CE_vec = c(rep(0,length(iter_vec_thin)))
    CE_weighed_vec = c(rep(0,length(iter_vec_thin)))
    MSE_vec = c(rep(0,length(iter_vec_thin)))
    MSE_weighed_vec = c(rep(0,length(iter_vec_thin)))
    TP_vec = c(rep(0,length(iter_vec_thin)))
    FP_vec = c(rep(0,length(iter_vec_thin)))
    TN_vec = c(rep(0,length(iter_vec_thin)))
    FN_vec = c(rep(0,length(iter_vec_thin)))
    plinks_max_diff = c(rep(0,length(iter_vec_thin)))
    
    #we compute the edge inclusion matrix at every single iteration
    for (g in 1:iter) {
      # for (g in c(iter)) {
      # # print progress
      # if (verbose == TRUE) {
      #   mes = paste(c("Calculating plinks ... in progress : ", floor(100 * 
      #                                                                   g/iter), "%"), collapse = "")
      #    cat(mes, "\r")
      #    utils::flush.console()
      #  }
      
      #compute and save posterior inclusion probability for every edge
      which_edge = which(unlist(strsplit(as.character(sample_graphs[all_graphs[g]]), "")) == 1)
      totaltime_per_edge[which_edge] = totaltime_per_edge[which_edge] + all_weights[g]
      result_plinks[,g] = totaltime_per_edge/sum(all_weights[c(1:g)])
      
    }
    
    #we calculate the AUC, CE and SVAUC values every x=thin iterations
    i = 0
    for (g in iter_vec_thin){
      
      # #print progress
      # if (verbose == TRUE) {
      #   mes = paste(c("Calculating metrics ... in progress : ", floor(100 * 
      #                                                                   g/iter), "%"), collapse = "")
      #   cat(mes, "\r")
      #   utils::flush.console()
      # }
      
      #count
      i = i + 1
      
      #obtain predictor
      predictor = result_plinks[,g]
      
      #calculate AUC metrics
      AUC_ROC_vec[i] = calc_AUC_ROC(predictor=predictor,response=response)   
      AUC_PR_vec[i] = calc_AUC_PR(predictor=predictor,response=response) 
      AUC_ROC_calibrated_vec[i] = calc_AUC_ROC_calibrated(predictor=predictor,response=response,cut_AUC_calibrated=cut_AUC_calibrated)
      
      #calculate CE and MSE metric
      CE_MSE_obj = calc_CE_MSE(predictor,response)
      CE_vec[i] = CE_MSE_obj$CE
      CE_weighed_vec[i] = CE_MSE_obj$CE_weighed
      MSE_vec[i] = CE_MSE_obj$MSE
      MSE_weighed_vec[i] = CE_MSE_obj$MSE_weighed
      
      #calculate TP,FP,TN,FN
      TP_FP_obj = calc_TP_FP(predictor,response,cutoff)
      TP_vec[i] =  TP_FP_obj$tp
      FP_vec[i] =  TP_FP_obj$fp
      TN_vec[i] =  TP_FP_obj$tn
      FN_vec[i] =  TP_FP_obj$fn
      
      #calculate plinks convergence
      if (g > plinks_diff){
        diff = abs(result_plinks[,g]-result_plinks[,g-plinks_diff])
        plinks_max_diff[i] = max(diff)
      }
      
      
    }
    
    #calculate final predictor
    predictor = result_plinks[,iter] 
    
    #calculate final AUC metrics
    AUC_ROC = calc_AUC_ROC(predictor=predictor,response=response)   
    AUC_PR = calc_AUC_PR(predictor=predictor,response=response) 
    AUC_ROC_calibrated = calc_AUC_ROC_calibrated(predictor=predictor,response=response,cut_AUC_calibrated=cut_AUC_calibrated)
    
    #calculate final CE and MSE metric
    CE_MSE_obj = calc_CE_MSE(predictor,response)
    CE = CE_MSE_obj$CE
    CE_weighed = CE_MSE_obj$CE_weighed
    MSE = CE_MSE_obj$MSE
    MSE_weighed = CE_MSE_obj$MSE_weighed
    
    #calculate final TP,FP,TN,FN
    TP_FP_obj = calc_TP_FP(predictor,response,cutoff)
    TP =  TP_FP_obj$tp
    FP =  TP_FP_obj$fp
    TN =  TP_FP_obj$tn
    FN =  TP_FP_obj$fn  
    
    #return all necessary values
    out = list(
      #final values
      plinks=predictor,
      AUC_ROC=AUC_ROC,
      AUC_PR=AUC_PR,
      # AUC_ROC_calibrated = AUC_ROC_calibrated,
      # CE = CE,
      # CE_weighed = CE_weighed,
      # MSE = MSE,
      # MSE_weighed = MSE_weighed,
      TP =  TP,
      FP =  FP,
      TN =  TN,
      FN =  FN,
      SteinL = SteinL,
      Fnorm = Fnorm
      
      # #values per iteration
      # plinks_vec=result_plinks[iter_vec_thin],
      # AUC_ROC_vec=AUC_ROC_vec,
      # AUC_PR_vec=AUC_PR_vec,
      # AUC_ROC_calibrated_vec = AUC_ROC_calibrated_vec,
      # CE_vec = CE_vec,
      # CE_weighed_vec = CE_weighed_vec,
      # MSE_vec = MSE_vec,
      # MSE_weighed_vec = MSE_weighed_vec,
      # TP_vec =  TP_vec,
      # FP_vec =  FP_vec,
      # TN_vec =  TN_vec,
      # FN_vec =  FN_vec,
      # plinks_max_diff = plinks_max_diff,
      # 
      # #iteration vector
      # iter_vec_thin=iter_vec_thin
    )
    
  } else if (type=="point_est") {
    
    predictor = abs(K_hat[upper.tri(K_hat)])
    
    #calculate final AUC metrics
    AUC_ROC = calc_AUC_ROC(predictor=predictor,response=response)   
    AUC_PR = calc_AUC_PR(predictor=predictor,response=response) 
    
    binary_predictor = K_hat[upper.tri(K_hat)]
    binary_predictor[binary_predictor != 0] = 1 
    
    #calculate final TP,FP,TN,FN
    TP_FP_obj = calc_TP_FP(binary_predictor,response,cutoff)
    TP =  TP_FP_obj$tp
    FP =  TP_FP_obj$fp
    TN =  TP_FP_obj$tn
    FN =  TP_FP_obj$fn  
    
    #return all necessary values
    out = list(
      #final values
      plinks=predictor,
      AUC_ROC=AUC_ROC,
      AUC_PR=AUC_PR,
      # AUC_ROC_calibrated = NULL,
      # CE = NULL,
      # CE_weighed = NULL,
      # MSE = NULL,
      # MSE_weighed = NULL,
      TP =  TP,
      FP =  FP,
      TN =  TN,
      FN =  FN,
      SteinL = SteinL,
      Fnorm = Fnorm
      
      # #values per iteration
      # plinks_vec=NULL,
      # AUC_ROC_vec=NULL,
      # AUC_PR_vec=NULL,
      # AUC_ROC_calibrated_vec = NULL,
      # CE_vec = NULL,
      # CE_weighed_vec = NULL,
      # MSE_vec = NULL,
      # MSE_weighed_vec = NULL,
      # TP_vec =  NULL,
      # FP_vec =  NULL,
      # TN_vec =  NULL,
      # FN_vec =  NULL,
      # plinks_max_diff = NULL,
      # 
      # #iteration vector
      # iter_vec_thin=NULL
    )
    
  } else if (type=="threshold_stat") {
    
    predictor = abs(threshold_stat[upper.tri(threshold_stat)])
    
    #calculate final AUC metrics based on z stat
    AUC_ROC = calc_AUC_ROC(predictor=predictor,response=response)   
    AUC_PR = calc_AUC_PR(predictor=predictor,response=response) 
    
    binary_predictor = K_hat[upper.tri(K_hat)]
    binary_predictor[binary_predictor != 0] = 1 
    
    #calculate final TP,FP,TN,FN
    TP_FP_obj = calc_TP_FP(binary_predictor,response,cutoff)
    TP =  TP_FP_obj$tp
    FP =  TP_FP_obj$fp
    TN =  TP_FP_obj$tn
    FN =  TP_FP_obj$fn  
    
    #return all necessary values
    out = list(
      #final values
      plinks=predictor,
      AUC_ROC=AUC_ROC,
      AUC_PR=AUC_PR,
      # AUC_ROC_calibrated = NULL,
      # CE = NULL,
      # CE_weighed = NULL,
      # MSE = NULL,
      # MSE_weighed = NULL,
      TP =  TP,
      FP =  FP,
      TN =  TN,
      FN =  FN,
      SteinL = SteinL,
      Fnorm = Fnorm
      
      # #values per iteration
      # plinks_vec=NULL,
      # AUC_ROC_vec=NULL,
      # AUC_PR_vec=NULL,
      # AUC_ROC_calibrated_vec = NULL,
      # CE_vec = NULL,
      # CE_weighed_vec = NULL,
      # MSE_vec = NULL,
      # MSE_weighed_vec = NULL,
      # TP_vec =  NULL,
      # FP_vec =  NULL,
      # TN_vec =  NULL,
      # FN_vec =  NULL,
      # plinks_max_diff = NULL,
      # 
      # #iteration vector
      # iter_vec_thin=NULL
    )
    
  } else if (type=="post_density") {
    #check if required information is given
    if (is.null(p_links))
      stop("please include the argument p_links")
    
    #obtain amount of entries
    qp = length(response)
    
    #calculate final predictor
    predictor = p_links
    
    #calculate final AUC metrics
    AUC_ROC = calc_AUC_ROC(predictor=predictor,response=response)   
    AUC_PR = calc_AUC_PR(predictor=predictor,response=response) 
    # AUC_ROC_calibrated = calc_AUC_ROC_calibrated(predictor=predictor,response=response,cut_AUC_calibrated=cut_AUC_calibrated)
    
    # #calculate final CE and MSE metric
    # CE_MSE_obj = calc_CE_MSE(predictor,response)
    # CE = CE_MSE_obj$CE
    # CE_weighed = CE_MSE_obj$CE_weighed
    # MSE = CE_MSE_obj$MSE
    # MSE_weighed = CE_MSE_obj$MSE_weighed
    
    #calculate final TP,FP,TN,FN
    binary_predictor <- as.numeric(predictor > 0.01) # for semms and gsemms
    TP_FP_obj = calc_TP_FP(binary_predictor,response,cutoff)
    TP =  TP_FP_obj$tp
    FP =  TP_FP_obj$fp
    TN =  TP_FP_obj$tn
    FN =  TP_FP_obj$fn  
    
    #return all necessary values
    out = list(
      #final values
      plinks=predictor,
      AUC_ROC=AUC_ROC,
      AUC_PR=AUC_PR,
      # AUC_ROC_calibrated = AUC_ROC_calibrated,
      # CE = CE,
      # CE_weighed = CE_weighed,
      # MSE = MSE,
      # MSE_weighed = MSE_weighed,
      TP =  TP,
      FP =  FP,
      TN =  TN,
      FN =  FN,
      SteinL = SteinL,
      Fnorm = Fnorm
      
      # #values per iteration
      # plinks_vec=result_plinks[iter_vec_thin],
      # AUC_ROC_vec=AUC_ROC_vec,
      # AUC_PR_vec=AUC_PR_vec,
      # AUC_ROC_calibrated_vec = AUC_ROC_calibrated_vec,
      # CE_vec = CE_vec,
      # CE_weighed_vec = CE_weighed_vec,
      # MSE_vec = MSE_vec,
      # MSE_weighed_vec = MSE_weighed_vec,
      # TP_vec =  TP_vec,
      # FP_vec =  FP_vec,
      # TN_vec =  TN_vec,
      # FN_vec =  FN_vec,
      # plinks_max_diff = plinks_max_diff,
      # 
      # #iteration vector
      # iter_vec_thin=iter_vec_thin
    )
  }
  
  return(out)
}

read_data = function(n=100,p=10,graph="random",rep=1,thin=100,
                     plinks_diff=1000,cut_AUC_calibrated=200,cutoff=0.5,
                     input_dir=NULL)
{
  
  #load data
  filename = paste0("result_p",p,"_n",n,"_",graph,"_prob", prob,"_rep",rep,".Rdata")
  load(file = paste0(input_dir, filename ))
  
  #obtain the true underlying data and save in the lsit obj_true
  actual = result$true_g 
  response = actual[upper.tri(actual)] #transform the adjacency matrix to a vector
  true_K = result$true_K
  true_sigma = result$true_sigma 
  true_data = result$true_data
  obj_true = list(actual=actual,response=response,true_K=true_K,true_sigma=true_sigma,true_data=true_data)
  
  # # ##########################
  # # ### ---bd method---###
  # # ##########################
  
  # #calculate all metrics and store them in the list obj_bd
  # # ignore K hat of bd method
  # obj_bd = calculate_metrics(all_graphs=result$all_graphs_bd,sample_graphs=result$sample_graphs_bd,
  #                            all_weights=result$all_weights_bd,response=response, thin = thin,
  #                            plinks_diff = plinks_diff, cut_AUC_calibrated=cut_AUC_calibrated,verbose = TRUE,
  #                            cutoff=cutoff, K_hat=true_K, true_K = true_K, type="mcmc")
  
  # #store information of all iterations to the list obj_bd
  # obj_bd$all_graphs = result$all_graphs_bd
  # obj_bd$sample_graphs = result$sample_graphs_bd
  # obj_bd$all_weights = result$all_weights_bd
  # obj_bd$K_hat = result$K_bd
  
  # #store time and interation data to the list obj_bd
  # obj_bd$time = result$time_bd
  # obj_bd$iter = length(obj_bd$all_graphs)
  
  # ##########################
  # ### ---ss method---###
  # ##########################
  # 
  # #calculate all metrics and store them in the list obj_ss
  # obj_ss = calculate_metrics(all_graphs=result$all_graphs_ss,sample_graphs=result$sample_graphs_ss,
  #                            all_weights=result$all_weights_ss,response=response, thin = thin,
  #                            plinks_diff = plinks_diff, cut_AUC_calibrated=cut_AUC_calibrated,verbose = TRUE,
  #                            cutoff=cutoff, K_hat=result$K_ss, true_K =  true_K, type="mcmc")
  # 
  # #store information of all iterations to the list obj_ss
  # obj_ss$all_graphs = result$all_graphs_ss
  # obj_ss$sample_graphs = result$sample_graphs_ss
  # obj_ss$all_weights = result$all_weights_ss
  # obj_ss$K_hat = result$K_ss
  # 
  # #store time and interation data to the list obj_ss
  # obj_ss$time = result$time_ss
  # obj_ss$iter = length(obj_ss$all_graphs)
  
  
  ##########################
  ### ---SEMMS method---###
  ##########################
  semms_plinks = result$plinks_semms[upper.tri(result$plinks_semms)]
  #calculate all metrics and store them in the list obj_bd
  obj_semms = calculate_metrics(response=response, cutoff=cutoff,
                                K_hat=result$K_semms, true_K=true_K,
                                type="post_density",
                                p_links=semms_plinks)
  obj_semms$K_hat = result$K_semms
  obj_semms$time = result$time_semms
  
  # ##########################
  # ### ---glasso method---###
  # ##########################
  
  # #calculate all metrics and store them in the list obj_glasso
  # obj_glasso = calculate_metrics(response=response, cutoff=cutoff,
  #                               K_hat=result$K_glasso, true_K =  true_K,
  #                               type="point_est")
  # obj_glasso$K_hat = result$K_glasso
  
  # #store time and interation data to the list obj_bd
  # obj_glasso$time = result$time_glasso
  
  ##########################
  ### ---Output method---###
  ##########################
  
  #create list to store output
  output_list = list( obj_true = obj_true,
                      # obj_bd=obj_bd,
                      # obj_ss=obj_ss,
                      obj_semms=obj_semms
                      # obj_glasso=obj_glasso
  )
  
  return(output_list)                      
}