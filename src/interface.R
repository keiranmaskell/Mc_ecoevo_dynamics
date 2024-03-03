#interface.R

#define working direcory
#this is the directory where the barcode count summary csv files are stored - not yet combined
if(readline("Use Default WD? (y/n)")=='n'){
  wd <- readline("Specify working directory filepath")
  setwd(wd)
}else{
  message("Using default working directory")
  system("cd */Mc_ecoevo_dynamics-main")
}

sprintf("Working directory set to %s",getwd())
system("ls")



source("src/init.R")

data_loader()


data_df <- make_df()

#apply the centered-logratio transform
data_df <- data_df %>%
    group_by(timerep_ID) %>%
    #mutate(CLR_transf = clr(Freqs))
    mutate(CLR_transf = yield_clr(Freqs))

#make a new df containing the means and sds across replicates of CLR-transformed and raw frequencies
summary_df <- summary_df(data_df)



test_unity_sum()