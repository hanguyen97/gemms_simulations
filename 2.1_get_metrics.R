library(BDgraph)
library(PRROC)
library(data.table)
library(StatPerMeCo)

dir = "/home/tn325/HCP/sbggm_final/"
# dir = "~/Documents/Cornell/Projects/HCP/Code/Simulations/sbggm_final/"

source(paste0(dir, "2.0_metric_helper.R"))
out_dir = paste0(dir, "outTEST/")


#--set parameters
graph_list = c("cluster") 
n_list = c(50, 100)
p_list = c(300)
rep_list = 1:30
thin = 100
plinks_diff = 1000
cut_AUC_calibrated = 200
cutoff=0.5


##---------------------------------------------------------------------------|
for (prob in c(0.05)) {
  for (p in p_list)
  {
    for( n in n_list )
    {
      for( graph in graph_list )
      {
        performance <- data.table(n=numeric(), 
                                  p=numeric(), 
                                  n_link=numeric(),
                                  rep=numeric(), 
                                  method=character(), 
                                  TP=numeric(), 
                                  FP=numeric(), 
                                  TN=numeric(), 
                                  FN=numeric(), 
                                  AUC_ROC=numeric(),
                                  AUC_PR=numeric(),
                                  Fnorm=numeric(), 
                                  time=numeric())
        
        selected_metrics <- c("TP", "FP", "TN", "FN", "AUC_ROC", "AUC_PR", "Fnorm", "time")
        # selected_metrics <- c("TP", "FP", "TN", "FN", "AUC_ROC", "AUC_PR", "time")
        
        for (rep in rep_list){
          
          # Save metrics 
          output_list = read_data(n=n,p=p,graph=graph,rep=rep,thin=thin,
                                  plinks_diff=plinks_diff,cut_AUC_calibrated=cut_AUC_calibrated,
                                  cutoff=cutoff, input_dir=out_dir)
          filename = paste0("ssmetrics_p",p,"_n",n,"_",graph,"_prob", prob,"_rep",rep,".Rdata")
          save(output_list, file = paste0(out_dir, filename )) 
          
          # Concatenate into a table
          result_list <- output_list[-1]
          filtered_list <- lapply(result_list, function(x) x[selected_metrics])
          
          # Convert the filtered list of lists to a data frame
          perf <- do.call(rbind, lapply(filtered_list, as.data.table))
          perf[, c("method", "n", "p", "rep", "n_link") := list(names(filtered_list), n, p, rep, sum(output_list$obj_true$response))]
          performance <- rbindlist(list(performance, perf), use.names=TRUE)
        }
        
        # Save all metrics 
        performance[, TPR := TP/(TP+FN)]
        performance[, FPR := FP/(FP+TN)]
        performance[, MCC := (TN*TP - FN*FP)/(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))^(1/2)) ]
        performance[, F1 := 2*TP/(2*TP +FP + FN)]
        write.csv(performance, paste0(out_dir, "performance_p",p,"_n",n,"_",graph,"_prob", prob,"_rep",rep,"_", Sys.Date(), ".csv"))
        
        # Save summary of metrics 
        perf_summary <- performance[, .(mean_density = mean(n_link)/choose(p, 2),
                                        n_sims = max(rep),
                                        # mean_SteinL = mean(SteinL),
                                        # sd_SteinL = sd(SteinL),
                                        # mean_Fnorm = mean(Fnorm),
                                        # sd_Fnorm = sd(Fnorm),
                                        mean_TPR = mean(TPR),
                                        sd_TPR = sd(TPR),
                                        mean_FPR = mean(FPR),
                                        sd_FPR = sd(FPR),
                                        mean_AUC_ROC = mean(AUC_ROC),
                                        sd_AUC_ROC = sd(AUC_ROC),
                                        mean_AUC_PR = mean(AUC_PR),
                                        sd_AUC_PR = sd(AUC_PR),
                                        mean_MCC = mean(MCC),
                                        sd_MCC = sd(MCC),
                                        # mean_F1 = mean(F1),
                                        # sd_F1 = sd(F1),
                                        mean_time = mean(time),
                                        sd_time = sd(time)), by=c("method", "n", "p")]
        
        write.csv(perf_summary, paste0(out_dir, "sssummary_p",p,"_n",n,"_",graph,"_prob", prob,"_rep",rep,"_", Sys.Date(), ".csv"))
        
      }
    }    
  }     
}

