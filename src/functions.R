#functions.R

#Centered log ratio transformation
yield_clr <- function(x){
    #zeroes must be replaced with a small value
    if(0 %in% x){
        epsilon <- 10^(-10)
        x_val <- x + epsilon
    }else{
        x_val <- x
    }
    values <- log(x_val/(exp(mean(log(x_val)))))
    return(values)
}


data_loader <- function(){

    files <- list.files()
    files_R1 <- files[grep(files, pattern="R1_paired")]
    files_R2 <- files[grep(files, pattern="R2_paired")]

    if(length(files_R1) != length(files_R2)){
        stop("number of forward and reverse reads must be equal; R1_paired = R2_paired")
    }
    

    if(!dir.exists("combined_csvs")){
        system(paste0('mkdir ',sprintf("%s/combined_csvs",getwd())))
    }
    output_dir <- sprintf("%s/combined_csvs",getwd())



# Combine files using lapply
    combined_files <- lapply(files_R1, function(file_R1) {
        index_current <- substr(file_R1, 1, (nchar(file_R1) - 13))

        # Read current R1 file and matching R2 file
        current_file_R1 <- read.csv(file_R1)
        matching_file_R2 <- read.csv(files_R2[files_R2 == sprintf("%sR2_paired.csv", index_current)])

        # Combine the R1 and R2 files
        combined_file <- matching_file_R2
        combined_file$Directionality <- "Both"
        combined_file$Read <- paste0(index_current, "R1_paired.csv & ", files_R2[files_R2 == sprintf("%sR2_paired.csv", index_current)])
        combined_file$Matches <- current_file_R1$Matches + matching_file_R2$Matches

        # Write the combined file to the output directory
        write.csv(combined_file, sprintf("%s/%scombined_R1R2.csv", getwd(), index_current))
        #write.csv(combined_file, sprintf("%s/%scombined_R1R2.csv", "/Users/keiranmaskell/Desktop/newcode_testdir", index_current))

        return(combined_file)
    })


    files <- list.files(output_dir)

    #test to ensure same number of rows in each output csv
    files <- list.files(output_dir, pattern = ".csv", full.names = TRUE)
    data_test <- lapply(files, read.csv)

    row_counts <- sapply(data_test, function(x) nrow(x))
    indices_of_exceptions <- which(row_counts != as.numeric(test_set))

    # Print results
    if (length(indices_of_exceptions) > 0) {
      cat("Warning: Dimensions differ within set. Exceptions found in files:", indices_of_exceptions, "\n")
    } else {
      cat("Pass: Dimensionality is equivalent.\n")
    }

    list.files(output_dir)


}

make_df <- function(){

    input_dir <- sprintf("%s/combined_csvs",getwd())

    files <- list.files(input_dir)[grep(list.files(input_dir), pattern = ".csv")]

    input_csvs <- lapply(sprintf("%s/%s",input_dir,files),read.csv)

    #the features, the timepoints, and the reps can be obtained at this stage
    #barcode_key <- data.frame(data_test[[1]]$Barcode)
    #timepoints <- unique(substr(files,1,2))
    #reps <- length(files[grep(files, pattern = timepoints[1])])

    #make one big df out of all the relevant csvs
    data_df <- do.call(rbind, input_csvs)

    #clean up df
    data_df <- data_df %>%
      dplyr::select(Effector_label, Barcode, Read, Matches) %>%
      #mutate(Read = substr(Read, 2, 2)) %>%
      mutate(Day = substr(Read, 1, 2)) %>%
      group_by(substr(Read, 2, 7)) %>%
      mutate(Freqs = Matches / sum(Matches)) %>%
      select(-Read)

    names(data_df) <- c("Effector", "BC","Abs_counts","Day","timerep_ID","Freqs")

    #replace timerep_ID with simple rep IDs 
    timepoints <- unique(substr(files,1,2))
    reps <- length(files[grep(files, pattern = timepoints[1])])
    #works bc data_df is ordered by timerep_ID already 
    data_df$pch <- rep(rep(seq(1,reps,1),each=length(unique(data_df$Effector))),length(unique(timepoints)))


    return(data_df)

}


summary_df <- function(datadf){

    mean_sd_combined <- data_df %>%
    group_by(Effector, Day)%>%
    #mutate(mean_CLR_by_rep = mean(CLR_transf))
    summarize(
            mean_CLR_byrep = mean(CLR_transf),
            sd_CLR_byrep = sd(CLR_transf),
            mean_freqs_byrep = mean(Freqs),
            sd_freqs_byrep = sd(Freqs)
            )
    return(mean_sd_combined)
}


# test_unity_sum <- function(data_df){
#     pf <- "PASS"
#     for(i in data_df$Day){
#         for(j in unique(data_df$timerep_ID[data_df$Day==i])){
#             if(as.numeric(sum(data_df$Freqs[data_df$Day==i & data_df$timerep_ID==j]))!=1 & all.equal(sum(data_df$Freqs[data_df$Day=="D0" & data_df$timerep_ID=="0-36-4"]), 1, tolerance = 1e-15)==FALSE){
#                 pf <- "FAIL"
#                 print(as.numeric(sum(data_df$Freqs[data_df$Day==i & data_df$timerep_ID==j])))
#                 print(sprintf("warning: frequencies do not sum to 1 at %s for rep %s",i,j))
#             }else{
#                 next
#                 print(sprintf("pass: frequencies sum to 1 at %s for rep %s",i, j))
#             }
#         }
#     }
#     print(sprintf("result of test for unity sum of frequencies across all samples is %s",pf))
# }