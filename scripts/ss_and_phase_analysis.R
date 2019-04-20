# ---
# Title: iCTC Project
# Description: Network creation and other funtions to be used in Network_script
# Contact: Kishore Hari <kishorehari@iisc.ac.in>
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in>, Kishore Hari <kishorehari@iisc.ac.in>, Mohit Kumar Jolly <mkjolly@iisc.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

library(tidyverse)
library(stringr)
library(magrittr)
library(ggthemes)

### This is the one used throughout in all sorts of analyses
# calculates the z-scores using the log values of all parameters combined together
z_score_calculator_log <- function(folder){#browser()
    curr <- getwd() # to reset the wd
    setwd(folder)
    
    solution_filez <- list.files(pattern = "solution.*dat")
    solution_filez <- solution_filez[!str_detect(solution_filez, "z")]
    
    infos <- file.info(solution_filez)$size
    empties <- which(is.na(infos) | infos == 0)
    if(length(empties) != 0)
        solution_filez <- solution_filez[-empties]
    
    prs <- read.delim(list.files(pattern = "prs", recursive = T)[1])
    comps <- prs$Parameter[which(str_detect(prs$Parameter, "Prod"))] %>% str_remove("Prod_of_")
    # Get the solutions, reorder them and normalize them
    df_list <- sapply(solution_filez, function(x){#browser()
        df <- read.delim(x, header = F)[, -c(1,2)]
        n_sol <- str_extract(x, "_\\d+\\.") %>% str_remove("_") %>% str_remove("\\.") %>% as.numeric 
        n_comp <- ncol(df)/n_sol
        namez <- paste0(rep(comps, each = n_sol), "_", 1:n_sol)
        col_order <- rep(1:n_comp, each = n_sol) + (0:(n_sol-1))*n_comp
        df <- df[, col_order]
        #df <- sapply(df, function(k){2^k}) %>% data.frame
        colnames(df) <- namez
        df
    })
    #browser()
    # get the means and variance of all componenets
    
    df_mv <- sapply(comps, function(x){
        d <- lapply(df_list, function(y){
            y[, which(str_detect(colnames(y),x))] %>% unlist
        }) %>% unlist %>% as.numeric
        c(mean(d, na.rm = T), sd(d, na.rm = T))
    })
    #browser()
    df_full <- lapply(comps, function(x){
        d <- lapply(df_list, function(y){
            y[, which(str_detect(colnames(y),x))] %>% unlist
        }) %>% unlist %>% as.numeric
    }) %>% reduce(cbind.data.frame)
    colnames(df_full) <- comps
    #browser()
    N_states <- str_extract(names(df_list), "_\\d+\\.") %>% str_remove("_") %>% str_remove("\\.") %>% as.integer()
    N_lengths <- sapply(df_list, nrow)*(N_states)
    df_full$n_states <- rep(N_states, N_lengths) %>% as.character()
    write_delim(df_full, paste0(folder, "_rfull.txt"), delim = "\t")
    
    # Calculate z-scores 
    df_zs <- sapply(df_list, function(x){#browser()
        x_names <- colnames(x)
        x <- sapply(x_names, function(y){ #browser()
            centrs <- df_mv[, which(str_detect(y, colnames(df_mv)))] %>% unlist
            (x[, y] - centrs[1])/centrs[2]
        }) %>% data.frame
        if(ncol(x) == 1) x %<>% t %<>% data.frame
        colnames(x) <- x_names
        x
    })
    
    z_files <- solution_filez %>% str_remove("\\.dat") %>% paste0("_lz.dat")
    sapply(1:length(df_zs), function(x){#browser()
        write_delim(df_zs[[x]], path = z_files[x], col_names = T, delim = "\t")
    })
    
    all_sols <- sapply(comps, function(x){
        d <- lapply(df_zs, function(y){
            y[, which(str_detect(colnames(y),x))] %>% unlist
        }) %>% unlist %>% as.numeric
    }) %>% data.frame
    
    colnames(all_sols) <- comps
    all_sols$n_states <- rep(N_states, N_lengths) %>% as.character()
    write_delim(all_sols, paste0(folder, "_lfull.txt"), delim = "\t")
    
    setwd(curr)
    #all_sols
}

# calculates the z-scores by normalizing all the steady states using production, degradation and lambda.
z_score_normalized_calculator <- function(folder){#browser()
    curr <- getwd() # to set it back later
    setwd(folder) 
    #browser()
    solution_filez <- list.files(pattern = "solution.*dat")
    solution_filez <- solution_filez[!str_detect(solution_filez, "z")]
    # To take care of empty files
    infos <- file.info(solution_filez)$size
    empties <- which(is.na(infos) | infos == 0)
    if(length(empties) != 0)
        solution_filez <- solution_filez[-empties]
    
    prs <- read.delim(list.files(pattern = "prs"))
    prs$Parameter <- as.character(prs$Parameter)
    n_comp <- sum(str_detect(prs$Parameter, "Prod")) 
    prods <- which(str_detect(prs$Parameter, "Prod"))
    degs <- which(str_detect(prs$Parameter, "Deg"))
    comp_names <- prs$Parameter[prods] %>% str_remove("Prod_of_")
    comps <- comp_names
    lamdas <- lapply(comp_names, function(x){#browser()
        c(which(str_detect(prs$Parameter, paste0("Inh.*",x, "$"))), which(str_detect(prs$Parameter, paste0("Act.*",x, "$"))))
    })
    
    
    # Get the g/k ratios for all components
    parameters <- read.delim(list.files(pattern = "parameter"), header = F)[, -(1:2)]
    parameters <- apply(parameters, 1, function(x){#browser()
        final <- numeric(length(lamdas))
        x <- as.numeric(x)
        for (i in 1:length(lamdas)){
            k <- lamdas[[i]]
            final[i] <- ifelse(length(k), prod(x[k]), 1)
            final[i] <- final[i]*x[prods[i]]/x[degs[i]]
        }
        final
    }) %>% t %>% data.frame
    
    # get the solutions, convert the log values to normal, divide with g/k and reorder them 
    df_list <- sapply(solution_filez, function(x){#browser()
        n_sol <- str_extract(x, "_\\d+\\.") %>% str_remove("_") %>% str_remove("\\.") %>% as.numeric 
        
        df <- read.delim(x, header = F)
        par_rowz <- df[,1] %>% as.numeric
        parz <- parameters[par_rowz, ]
        
        df <- df[, -(1:2)]
        df <- sapply(df, function(k){2^k}) %>% data.frame
        if(ncol(df) == 1) df %<>% t %<>% data.frame
        
        for (i in 1:n_sol)
        {
            cols <- n_comp*(i-1) + 1:n_comp
            df[,cols] <- df[,cols]/parz
            df[,cols] <- log2(df[,cols])
        }
        
        namez <- paste0(rep(comps, each = n_sol), "_", 1:n_sol)
        col_order <- rep(1:n_comp, each = n_sol) + (0:(n_sol-1))*n_comp
        df <- df[, col_order]
        
        colnames(df) <- namez
        df
    })
    
    # collect all the solution values for each component, get the mean and variance
    #comps <- colnames(df_list[[1]]) %>% str_extract("X\\d")
    df_mv <- sapply(comps, function(x){
        d <- lapply(df_list, function(y){
            y[, which(str_detect(colnames(y),x))] %>% unlist
        }) %>% unlist %>% as.numeric
        c(mean(d, na.rm = T), sd(d, na.rm = T))
    })
    
    df_zs <- sapply(df_list, function(x){#
        #browser()
        x_names <- colnames(x)
        x <- sapply(x_names, function(y){ #browser()
            centrs <- df_mv[, which(str_detect(y, colnames(df_mv)))] %>% unlist
            (x[, y] - centrs[1])/centrs[2]
        }) %>% data.frame
        if (ncol(x) == 1) x %<>% t %<>% data.frame
        colnames(x) <- x_names
        x
    })
    
    z_files <- solution_filez %>% str_remove("\\.dat") %>% paste0("_nz.dat")
    sapply(1:length(df_zs), function(x){#browser()
        write_delim(df_zs[[x]], path = z_files[x], col_names = T, delim = "\t")
    })
    
    all_sols <- sapply(comps, function(x){
        d <- lapply(df_zs, function(y){
            y[, which(str_detect(colnames(y),x))] %>% unlist
        }) %>% unlist %>% as.numeric
    }) %>% data.frame
    
    colnames(all_sols) <- comps
    
    N_states <- str_extract(names(df_list), "_\\d+\\.") %>% str_remove("_") %>% str_remove("\\.") %>% as.integer()
    N_lengths <- sapply(df_list, nrow)*(N_states)
    all_sols$n_states <- rep(N_states, N_lengths) %>% as.character()
    
    write_delim(all_sols, paste0(folder, "_nfull.txt"), delim = "\t")
    
    setwd(curr)
    #all_sols
}



# the folder must contain 3 sub_folders that are repeats
## Options:
#       folder : the path of the folder in which the triplet results are stored
#       solution : which solution files to chose. Possible options are:
#                   1) a number - functions picks up the systems with these many solutions
#                   2) "all" - applies the log calculator and picks up lz files
#                   3) "normalized" - aplies the normalized log calculator and picks up the nz files
#       n_cluster : number of clusters to divide the data into - passed on the clustering algorithms
#       z_score : T or F. if T, it picks up dat files containing z-scores. else, it picks up normal files.
#       folder_name : name of the folder to be passed on to the resulting PDF files. if blank, extracts the name from folder name.
#       clustering_method : HCA or K-means. Currently only works with HCA
#       scatter : components to be used to derive the scatter plot
#       color_by : 1 or 2 - color the scatter plots by clusters or number of states
#       comp_delete : natural number less than or equal to number of components. Decides if any componenet needs to be eliminated for cluster analysis
# color_by option : if 1, coloring in scatter plot based on clusters. if 2, coloring based on number of states the system had
steady_state_patterns <- function(folder, solution = "all", n_cluster, z_score = T, folder_name = "", clustering_method = 1, 
                                  scatter = c(1,2), color_by = 1, comp_delete  = 0)
{
    
    curr <- getwd()
    setwd(folder)
    
    dir1 <- list.dirs(recursive = F)
    dir1 <- dir1[-which(str_detect(dir1, "RACIPE"))]
    if(solution == "all") sapply(dir1, function(y){
        f <- list.files(path = y, pattern = "lz")
        if(length(f) == 0)
            z_score_calculator_log(y)
    })
    if(solution == "normalized") sapply(dir1, function(y){
        f <- list.files(path = y, pattern = "nz")
        if(length(f) == 0)
            z_score_normalized_calculator(y)
    })
    
    ## Set pattern for the selected files and the names for pdf files
    header <- T
    name <- ""
    n <- as.numeric(solution)
    if(is.na(n)){
        name <- paste0(name, "_full")
        if(z_score){
            if (solution == "all") pattern <- "lfull"
            if (solution == "normalized") pattern <- "nfull"
            name <- paste0(name, "_Z")
            if (solution == "normalized") name <- paste0(name, "_normalized")
        }
        else pattern <- "rfull"
    }
    
    else {
        name <- paste0(name, "_solution_", n)
        if(z_score){
            pattern <- paste0("solution_", n, "_lz")
            
        }
        else{
            pattern <- paste0("solution_", n, ".dat")
            header <- F
        }
    }
    #browser()
    f <- list.files(pattern = pattern, recursive = T)
    if(f == "" || length(f) == 0) {
        message(paste0("Invalid parameters in the folder ", folder,". Please check the parameters and try again."))
        return()
    }
    soln_clustered <- lapply(f, function(x){#browser()
        d <- read.delim(x, header = header)
        if(!is.na(n)){#browser()
            if(!z_score) d <- d[, -(1:2)]
            if(n > 1)
            {
                n_comp <- ncol(d)/n
                if(!is.na(n) && z_score)
                {
                    df <- lapply(1:n, function(y){
                        k <- d[, ((0:(n_comp-1))*n + rep(1, n_comp) + (y-1))]
                        colnames(k) <- 1:n_comp
                        k
                    }) %>% reduce(rbind.data.frame)
                    colnames(df) <- c(paste0("X", 1:n_comp), "n_states")
                }
                else
                {
                    df <- lapply(1:n, function(y){
                        k <- d[, ((y-1)*n_comp + 1:n_comp)]
                        colnames(k) <- 1:n_comp
                        k
                    }) %>% reduce(rbind.data.frame)
                    colnames(df) <- paste0("X", 1:n_comp)
                }
                d <- df
            }
            d$n_states <- str_extract(x, "_\\d+\\.") %>% str_remove("_") %>% str_remove("\\.")
        }
        #browser()
        if (clustering_method == 1)
        {
            if(comp_delete >0 && comp_delete <ncol(d)) d1 <- d[, -comp_delete]
            else d1 <- d
            d$Cluster <- d1[, -which(colnames(d1) == "n_states")] %>% dist %>% hclust %>% cutree(k = n_cluster)
        }
        
        
        d <- d[order(d$Cluster, decreasing = F),]
        d$num <- 1:nrow(d)
        d$Cluster %<>% as.character
        
        d
    })
    
    if(folder_name == "") folder_name <- folder %>% str_remove("\\/$") %>% str_extract("\\/(?:.(?!\\/))+$") %>% str_remove("\\.") %>% str_remove("\\/")
    folder <- folder_name
    
    colr <- ifelse(color_by == 1, "Cluster", "n_states")
    name <- paste0(name, "_", n_cluster, "_clusters")
    if(comp_delete >0 && comp_delete <ncol(soln_clustered[[1]])) folder <- paste0(folder, "_X", comp_delete, "_deleted")
    
    
    heatmaps <- lapply(soln_clustered, function(d){#browser()
        #colnames(d)[1:2] <- c("A", "B")
        gatherd <- d %>% gather(key = "Component", value = "Z_score", -Cluster, -num, -n_states)
        gatherd$Cluster <- as.character(gatherd$Cluster)
        #gatherd$Component <- as.factor(gatherd$Component)
        #levels(gatherd$Component) <- c("CDH1", "VIM")
        #gatherd$Component <- as.character(gatherd$Component)
        ggplot(gatherd, aes(x = Component, y = num, fill = Z_score)) + geom_tile(width = 0.9) + scale_fill_gradientn(colours = rev(rainbow(5))) + 
            theme_stata() + 
            theme(axis.text.y = element_blank(), line = element_blank(), axis.title.y = element_blank(), legend.position = "right", axis.text.x = element_text(size = 29)) + 
            labs(x = "")
    })
    #browser()
    pdf(paste0(folder, "_Heatmaps", name, ".pdf"))
    sapply(heatmaps, print)
    dev.off()
    
    
    
    scatter_plots <- lapply(soln_clustered, function(d){
        d$n_states %<>% as.character
        d$Cluster %<>% as.character
        ggplot(d, aes_string(x = colnames(d)[scatter[1]], y = colnames(d)[scatter[2]], color = colr)) + geom_point() + geom_jitter() +
            theme_stata()
    })
    
    pdf(paste0(folder, "_ScatterPlots", name, "_color_", colr,"_X",scatter[1], "_X",scatter[2], ".pdf"))
    sapply(scatter_plots, print)
    dev.off()
    
    if(curr != getwd()) system(paste0("cp *pdf ", curr))
    
    setwd(curr)
    
}

frequency_plot_maker <- function(folder, n_cluster = 3, primary_comps = c(1,2)){#browser()
    curr <- getwd()
    setwd(folder)
    filez <- list.files(pattern = "lfull", recursive = T)
    prs <- read.delim(list.files(pattern = "prs", recursive = T)[1])
    comps <- prs$Parameter[which(str_detect(prs$Parameter, "Prod"))] %>% str_remove("Prod_of_")
    
    n_factor <- sapply(filez, function(x){
        nrow(read.delim(x))
    }) %>% max
    
    n1 <- primary_comps[1]
    n2 <- primary_comps[2]
    
    parms <- list.files(pattern = "parameters", recursive = T)
    
    NSS_frequency_data <- lapply(parms, function(x){
        d <- read.delim(x)
        d[,2] %>% unlist %>% table
    })
    
    lengths <- sapply(NSS_frequency_data,length)
    n_states <- names(NSS_frequency_data[[which.max(lengths)]])
    
    NSS_frequency_data <- lapply(NSS_frequency_data, function(x){#browser()
        for (j in n_states){
            if (is.na(x[j])) x[j] <- 0
        }
        k <- names(x)[which(names(x) %in% n_states)]
        x <- x[k]
        x <- x[order(as.numeric(names(x)))]
        as.numeric(x)
    }) %>% reduce(cbind.data.frame)
    
    NSS_frequency_data$Mean <- apply(NSS_frequency_data, 1, mean)
    sd <- NSS_frequency_data %>% apply(1, sd)
    NSS_frequency_data$Min <- NSS_frequency_data$Mean - sd
    NSS_frequency_data$Max <- NSS_frequency_data$Mean + sd
    NSS_frequency_data$n_states <- n_states
    colnames(NSS_frequency_data)[1:3] <- c("A", "B", "C")
    
    plots <- lapply(filez, function(x){
        d <- read.delim(x)
        d$Quadrant <- "1"
        d$Quadrant[d[[n1]]*d[[n2]] < 0 & d[[n1]] > 0] <- "4"
        d$Quadrant[d[[n1]]*d[[n2]] < 0 & d[[n1]] < 0] <- "2"
        d$Quadrant[d[[n1]]*d[[n2]] >0 & d[[n1]] < 0] <- "3"
        ggplot(d, aes(x= Quadrant)) + geom_bar(aes(y = ..prop.., group = 1)) + theme_stata()
    })
    
    Quadrant_combined <- lapply(filez, function(x){
        d <- read.delim(x)
        d$Quadrant <- "1"
        d$Quadrant[d[[n1]]*d[[n2]] < 0 & d[[n1]] > 0] <- "4"
        d$Quadrant[d[[n1]]*d[[n2]] < 0 & d[[n1]] < 0] <- "2"
        d$Quadrant[d[[n1]]*d[[n2]] >0 & d[[n1]] < 0] <- "3"
        table(d$Quadrant) %>% as.numeric
    }) %>% reduce(cbind.data.frame)
    Quadrant_combined <- Quadrant_combined/colSums(Quadrant_combined)
    Quadrant_combined$Mean <- apply(Quadrant_combined, 1, mean)
    sd <- Quadrant_combined %>% apply(1, sd)
    Quadrant_combined$Min <- Quadrant_combined$Mean - sd
    Quadrant_combined$Max <- Quadrant_combined$Mean + sd
    Quadrant_combined$Quadrant <- 1:4 %>% as.character
    colnames(Quadrant_combined)[1:3] <- c("A", "B", "C")
    
    Quadrant_states <- lapply(filez, function(x){
        d <- read.delim(x)
        d$Quadrant <- "1"
        d$Quadrant[d[[n1]]*d[[n2]] < 0 & d[[n1]] > 0] <- "4"
        d$Quadrant[d[[n1]]*d[[n2]] < 0 & d[[n1]] < 0] <- "2"
        d$Quadrant[d[[n1]]*d[[n2]] >0 & d[[n1]] < 0] <- "3"
        group_by(d, Quadrant, n_states) %>% summarise(Count = n())
    }) %>% reduce(merge, by = c("Quadrant", "n_states"))
    Quadrant_states[, 3:5] <- Quadrant_states[, 3:5]/colSums(Quadrant_states[, 3:5])
    Quadrant_states$Mean <- apply(Quadrant_states[, 3:5], 1, mean)
    sd <- Quadrant_states[, 3:5] %>% apply(1, sd)
    Quadrant_states$Min <- Quadrant_states$Mean - sd
    Quadrant_states$Max <- Quadrant_states$Mean + sd
    #Quadrant_states$Quadrant <- 1:4 %>% as.character
    colnames(Quadrant_states)[3:5] <- c("A", "B", "C")
    Quadrant_states$n_states %<>% as.character
    
    cluster_plots_data <- lapply(filez, function(x){
        d <- read.delim(x)
        d$Cluster <- d[, -ncol(d)] %>% dist %>% hclust %>% cutree(k = n_cluster)
        d
    })
    
    clusters <- cluster_plots_data %>% lapply(function(x){
        as.numeric(table(x$Cluster))
    }) %>% reduce(cbind.data.frame)
    clusters <- clusters/colSums(clusters)
    clusters$Mean <- apply(clusters, 1, mean)
    sd <- clusters %>% apply(1, sd)
    clusters$Min <- clusters$Mean - sd
    clusters$Max <- clusters$Mean + sd
    clusters$Cluster <- 1:n_cluster %>% as.character
    colnames(clusters)[1:3] <- c("A", "B", "C")
    
    cluster_plots <- lapply(cluster_plots_data, function(d){#browser()
        ggplot(d, aes (x = Cluster)) + geom_bar(aes(y = ..prop.., group = 1)) + theme_stata() +
            labs(y = "Frequency")
    })
    
    #browser()
    cluster_states <- cluster_plots_data %>% lapply(function(x){
        x %>% group_by(n_states, Cluster) %>% summarise(Count = n())
    }) %>% reduce(merge, by = c("Cluster", "n_states"))
    cluster_states[, 3:5] <- cluster_states[, 3:5]/colSums(cluster_states[,3:5])
    cluster_states$Mean <- apply(cluster_states[, 3:5], 1, mean)
    sd <- cluster_states[, 3:5] %>% apply(1, sd)
    cluster_states$Min <- cluster_states$Mean - sd
    cluster_states$Max <- cluster_states$Mean + sd
    #cluster_states$Cluster <- 1:2 %>% as.character
    colnames(cluster_states)[3:5] <- c("A", "B", "C")
    cluster_states$n_states %<>% as.character
    
    
    cluster_plots_state_wise <- lapply(cluster_plots_data, function(d){#browser()
        d$n_states %<>% as.character
        ggplot(d, aes (x = Cluster)) + geom_bar(aes(fill = n_states, y = ..count../sum(..count..)), position = position_dodge()) +
            theme_stata() + labs(y = "Frequency")
    })
    
    corr_plots <- lapply(filez, function(x){
        d <- read.delim(x)
        if (!all(colnames(d) %in% comps)) {
            cs <- which(!colnames(d) %in% comps)
            d <- d[, -cs]
        }
        
        colnames(d) <- comps
        corr <- cor(d) %>% data.frame
        corr$Comps1 <- colnames(corr)
        #browser()
        corr <- gather(corr, key = "Comps2", value = "Cor_Pearson", -Comps1)
        corr$Cor_Pearson %<>% round(3)
        m <- min(corr$Cor_Pearson)
        M <- max(corr$Cor_Pearson)
        ggplot(corr, aes(x = Comps1, y = Comps2)) + 
            geom_tile(aes(fill = Cor_Pearson)) + geom_text(aes(label = Cor_Pearson)) + 
            scale_fill_gradient_tableau(breaks = round(seq(m, M, length.out = 7), 3))
    })
    
    folder_name <- folder %>% str_remove("\\/$") %>% str_extract("\\/(?:.(?!\\/))+$") %>% str_remove("\\.") %>% str_remove("\\/")
    
    setwd(curr)
    
    pdf(paste0(folder_name, "_Quadrant_state_frequency_plots.pdf"))
    print(ggplot(Quadrant_states, aes(x = Quadrant, y = A, fill = n_states)) + geom_bar(stat = "identity", position = position_dodge()) + 
              labs (y = "Frequency") + theme_stata())
    print(ggplot(Quadrant_states, aes(x = Quadrant, y = B, fill = n_states)) + geom_bar(stat = "identity", position = position_dodge()) + 
              labs (y = "Frequency") + theme_stata())
    print(ggplot(Quadrant_states, aes(x = Quadrant, y = B, fill = n_states)) + geom_bar(stat = "identity", position = position_dodge()) + 
              labs (y = "Frequency") + theme_stata())
    print(ggplot(Quadrant_states, aes(x = Quadrant, y = Mean, fill = n_states)) + geom_bar(stat = "identity", position = position_dodge()) + 
              geom_errorbar(aes(ymin = Min, ymax = Max), width=.2,
                            position=position_dodge(.9)) + labs (y = "Frequency") + theme_stata())
    dev.off()
    
    pdf(paste0(folder_name,"_Quadrant_frequency_plots.pdf"))
    sapply(plots, print)
    print(ggplot(Quadrant_combined, aes(x = Quadrant, y = Mean)) + geom_bar(stat = "identity") + 
              geom_errorbar(aes(ymin = Min, ymax = Max), width=.2, position=position_dodge(.9)) +
              labs(y = "Frequency") + theme_stata())
    dev.off()
    
    pdf(paste0(folder_name,"_state_frequency_plots.pdf"))
    print(ggplot(NSS_frequency_data, aes(x = n_states, y = A/10000)) + geom_bar(stat = "identity") + labs (y = "Frequency") + theme_stata())
    print(ggplot(NSS_frequency_data, aes(x = n_states, y = B/10000)) + geom_bar(stat = "identity") + labs (y = "Frequency") + theme_stata())
    print(ggplot(NSS_frequency_data, aes(x = n_states, y = B/10000)) + geom_bar(stat = "identity") + labs (y = "Frequency") + theme_stata())
    print(ggplot(NSS_frequency_data, aes(x = n_states, y = Mean/10000)) + geom_bar(stat = "identity") + 
              geom_errorbar(aes(ymin = Min/10000, ymax = Max/10000), width=.2,
                            position=position_dodge(.9)) + labs (y = "Frequency") + theme_stata())
    dev.off()
    
    pdf(paste0(folder_name, "_Cluster_frequency_plots_", n_cluster, "_clusters.pdf"))
    sapply(cluster_plots, print)
    print(ggplot(clusters, aes(x = Cluster, y = Mean)) + geom_bar(stat = "identity") + 
              geom_errorbar(aes(ymin = Min, ymax = Max), width=.2, position=position_dodge(.9)) +
              labs(y = "Frequency") + theme_stata())
    dev.off()
    #browser()
    pdf(paste0(folder_name, "_Cluster_frequency_with_state_frequency_plots_", n_cluster, "_clusters.pdf"))
    sapply(cluster_plots_state_wise, print)
    print(ggplot(cluster_states, aes(x = Cluster, y = Mean, fill = n_states)) + 
              geom_bar(stat = "identity", position = position_dodge()) + 
              geom_errorbar(aes(ymin = Min, ymax = Max), width=.2, position=position_dodge(.9)) +
              labs(y = "Frequency") + theme_stata())
    dev.off()
    
    pdf(paste0(folder_name, "_Correlation_plot.pdf"))
    sapply(corr_plots, print)
    dev.off()
    
    
}

## Here, primary should be such that the first element indicates zeb and second indicates mir200
Quadrant_and_phase_calc <- function(folder, primary = c(1,2), zs = "log"){#browser()
    curr <- getwd()
    setwd(folder)
    folders <- list.dirs(recursive = F)
    n1 <- primary[1]
    n2 <- primary[2]
    z1 <- ifelse(zs == "log", "lz", "nz")
    dummy <- sapply(folders, function(x){#browser()
        setwd(x)
        soln_files <- list.files(pattern = "solution")
        soln_files <- soln_files[str_detect(soln_files, z1)]
        comps <- paste0("X", 1:(ncol(read.delim(soln_files[1]))))
        n_comp <- length(comps)
        df_zs <- lapply(soln_files, function(y){
            d <- read.delim(y)
            d <- d[complete.cases(d), ]
            n_sols <- str_extract(y, "solution_\\d+_") %>% str_remove("solution_") %>% str_remove("_") %>% as.numeric
            df <- apply(d, 1, function(z){
                #print(z)
                st <- lapply(1:n_sols, function(w){
                    z %<>% as.numeric
                    z[(w + cumsum(c(0,rep(n_sols, n_comp-1))))]
                })
                phase <- sapply(st, function(w){
                    if(w[n1]*w[n2] >0 && w[n1] >0) "1"
                    else if (w[n1]*w[n2] >0 && w[n1] <0) "3"
                    else if (w[n1]*w[n2] <0 && w[n1] <0) "2"
                    else "4"
                })
                if(length(st) == 1) st <- st[[1]] %>% t %>% data.frame
                else st %<>% reduce(rbind.data.frame)
                colnames(st) <- comps
                st$Quadrant <- phase
                Phase_dict <- c("H", "E", "L", "M")
                names(Phase_dict) <- as.character(1:4)
                phase1 <- sapply(names(Phase_dict), function(w){
                    paste0("_", Phase_dict[w],sum(phase == w))
                }) %>% paste(collapse = "") %>% str_remove("_")
                phase2 <- sapply(names(Phase_dict), function(w){
                    g <- sum(phase == w)
                    if(g>0) Phase_dict[w]
                    else NULL
                    #paste0("_", Phase_dict[w],sum(phase == w))
                }) %>% paste(collapse = "") %>% str_remove_all("NULL") 
                
                
                st$Phase <- phase1
                st$short_Phase <- phase2
                st
            }) %>% reduce(rbind.data.frame)
            df$n_states <- n_sols
            df
        }) %>% reduce(rbind.data.frame)
        write_delim(df_zs, paste0(zs,"_full_qpn.txt"), delim = "\t")
        
        setwd("..")
    })
    setwd(curr)
}

Quadrant_plots <- function(folders, primary = c(1,2), redo = F){
    # if(length(folders) < 2) {
    #     message("Minimum 2 folders needed")
    #     return()
    # }
    curr <- getwd()
    Qdf <- lapply(folders, function(x){
        name <- x %>% str_remove("\\/$") %>% str_extract("\\/(?:.(?!\\/))+$") %>% str_remove("\\.") %>% str_remove("\\/")
        setwd(x)
        f <- list.files(pattern = "qpn", recursive = T)
        if(length(f) < 3 || redo)
        {
            if(length(list.files(pattern = "lz", recursive = T)) < 3) sapply(list.dirs(recursive = F), z_score_calculator_log)
            setwd("..")
            Quadrant_and_phase_calc(x, primary)
            setwd(x)
        }
        f <- list.files(pattern = "qpn", recursive = T)
        qdf <- lapply(f, function(y){
            d <- read.delim(y)
            p <- table(d$Quadrant) %>% as.numeric
            p/nrow(d)
        }) %>% reduce(cbind.data.frame)
        Mean <- apply(qdf, 1, mean)
        Sd <- apply(qdf, 1, sd)
        qdf$Mean <- Mean
        qdf$Min <- Mean-Sd
        qdf$Max <- Mean+Sd
        qdf <- qdf[, -(1:length(f))]
        colnames(qdf) %<>% paste0(name, "_", .)
        setwd(curr)
        qdf
    }) %>% reduce(cbind.data.frame)
    #browser()
    Qdf$Quadrants <- as.character(1:4)
    
    Qdf1 <- Qdf %>% 
        select(Quadrants, contains("Mean")) %>%
        gather(key = "Circuit", value = "Mean", -Quadrants)
    Qdf1$Circuit %<>% str_remove("_Mean")
    
    Qdf2 <- Qdf %>% 
        select(Quadrants, contains("Min")) %>%
        gather(key = "Circuit", value = "Min", -Quadrants)
    Qdf2$Circuit %<>% str_remove("_Min")
    
    Qdf3 <- Qdf %>% 
        select(Quadrants, contains("Max")) %>%
        gather(key = "Circuit", value = "Max", -Quadrants)
    Qdf3$Circuit %<>% str_remove("_Max")
    
    Qlist <- list(Qdf1, Qdf2, Qdf3)
    
    Qdf <- reduce(Qlist, merge, by = c("Quadrants", "Circuit"))
    names <- folders %>% str_remove("\\/$") %>% str_extract("\\/(?:.(?!\\/))+$") %>% str_remove("\\.") %>% str_remove("\\/")
    
    if(length(unique(Qdf$Circuit)) > 1)
        p <- ggplot(Qdf, aes(x = Quadrants, y = Mean, fill = Circuit)) + geom_bar(stat = "identity", position = position_dodge()) + 
        geom_errorbar(aes(ymin = Min, ymax = Max), position = position_dodge(0.9)) + theme_stata() + labs(y = "Frequency")
    
    else
        p <- ggplot(Qdf, aes(x = Quadrants, y = Mean)) + geom_bar(stat = "identity", position = position_dodge()) + 
        geom_errorbar(aes(ymin = Min, ymax = Max), position = position_dodge(0.9)) + theme_stata() + labs(y = "Frequency")
    
    setwd(curr)
    
    pdf(paste0(paste(names, collapse = "_"), "_quadrant.pdf"))
    print(p)
    dev.off()
}

Phase_plots <- function(folders, primary = c(1,2), redo = F, min_phase_freq = 0, rank_by_freq = Inf, short_code = F, zs = "log"){
    # if(length(folders) < 2) {
    #     message("Minimum 2 folders needed")
    #     return()
    # }
    column <- "Phase"
    if(short_code) column = "short_Phase"
    curr <- getwd()
    Qdf <- lapply(folders, function(x){
        name <- x %>% str_remove("\\/$") %>% str_extract("\\/(?:.(?!\\/))+$") %>% str_remove("\\.") %>% str_remove("\\/")
        setwd(x)
        f <- list.files(pattern = paste0(zs,"_full_qpn"), recursive = T)
        if(length(f) < 3 || redo)
        {
            if(length(list.files(pattern = "lz", recursive = T)) < 3) sapply(list.dirs(recursive = F), z_score_calculator_log)
            if(length(list.files(pattern = "nz", recursive = T)) < 3) sapply(list.dirs(recursive = F), z_score_normalized_calculator)
            setwd("..")
            Quadrant_and_phase_calc(x, primary, zs = zs)
            setwd(x)
        }
        f <- list.files(pattern = paste0(zs,"_full_qpn"), recursive = T)
        qdf <- lapply(f, function(y){
            d <- read.delim(y)
            p <- table(d[column]) %>% data.frame
            p$Freq <- p$Freq/nrow(d)
            colnames(p) <- c("Phase", str_extract(y, "\\d"))
            p
        }) %>% reduce(merge , by = "Phase", all = T)
        qdf[is.na(qdf)] <- 0
        Mean <- apply(qdf[,-1], 1, mean)
        Sd <- apply(qdf[, -1], 1, sd)
        qdf$Mean <- Mean
        qdf$Min <- Mean-Sd
        qdf$Max <- Mean+Sd
        qdf <- qdf[, -(2:(length(f)+1))]
        colnames(qdf)[-1] %<>% paste0(name, "_", .)
        setwd(curr)
        qdf
    }) %>% reduce(merge, by = "Phase", all = T)
    #browser()
    Qdf[is.na(Qdf)] <- 0
    
    Qdf1 <- Qdf %>% 
        select(Phase, contains("Mean")) %>%
        gather(key = "Circuit", value = "Mean", -Phase)
    Qdf1$Circuit %<>% str_remove("_Mean")
    
    Qdf2 <- Qdf %>% 
        select(Phase, contains("Min")) %>%
        gather(key = "Circuit", value = "Min", -Phase)
    Qdf2$Circuit %<>% str_remove("_Min")
    
    Qdf3 <- Qdf %>% 
        select(Phase, contains("Max")) %>%
        gather(key = "Circuit", value = "Max", -Phase)
    Qdf3$Circuit %<>% str_remove("_Max")
    
    Qlist <- list(Qdf1, Qdf2, Qdf3)
    #browser()
    Qdf <- reduce(Qlist, merge, by = c("Phase", "Circuit"))
    Qdf <- Qdf %>% select(Phase, Mean) %>% group_by(Phase) %>% summarise(total = sum(Mean)) %>% mutate(Rank = dense_rank(desc(total))) %>% merge(Qdf, by = "Phase")
    #browser()
    
    names <- folders %>% str_remove("\\/$") %>% str_extract("\\/(?:.(?!\\/))+$") %>% str_remove("\\.") %>% str_remove("\\/") %>% paste(collapse = "_")
    
    if(short_code) names %<>% paste0("_Short_codes")
    
    if (!is.na(as.numeric(min_phase_freq)) && min_phase_freq > 0 && min_phase_freq < 1) 
    {
        Qdf %<>% filter(total >= min_phase_freq)
        names %<>% paste0("_min_freq_", min_phase_freq)
    }
    
    if (!is.infinite(rank_by_freq) && is.numeric(rank_by_freq))
    {
        Qdf %<>% filter(Rank <= rank_by_freq)
        names %<>% paste0("_max_freq_rank_", rank_by_freq)
    }
    if(short_code)
    {
        l_row <- which(Qdf$Phase == "L")
        lh_rows <- which(str_detect(Qdf$Phase, "L"))
        
        Qdf_noL <- Qdf[-l_row,]
        Qdf_noL$Phase[lh_rows] <- Qdf_noL$Phase[lh_rows] %>% str_remove("L")
    }
    else
    {
        l_row <- which(str_detect(Qdf$Phase, "^L\\d?$"))
        lh_rows <- which(str_detect(Qdf$Phase, "L"))
        
        Qdf_noL <- Qdf[-l_row,]
        Qdf_noL$Phase[lh_rows] <- Qdf_noL$Phase[lh_rows] %>% str_remove("L\\d?")
    }
    
    Qdf_noL <- Qdf_noL %>% group_by(Phase) %>% summarise_if(is.numeric, sum) %>% mutate(Rank = dense_rank(desc(total)))
    Qdf_noL$Circuit <- unique(Qdf$Circuit)
    
    if(length(unique(Qdf$Circuit)) > 1)
    {
        p <- ggplot(Qdf, aes(x = Phase, y = Mean, fill = Circuit)) + geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
            geom_errorbar(aes(ymin = Min, ymax = Max), position = position_dodge(0.9)) + theme_stata() + labs(y = "Frequency") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        p1 <- ggplot(Qdf_noL, aes(x = Phase, y = Mean, fill = Circuit)) + geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
            geom_errorbar(aes(ymin = Min, ymax = Max), position = position_dodge(0.9)) + theme_stata() + labs(y = "Frequency") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
        
    else
    {
        p <- ggplot(Qdf, aes(x = Phase, y = Mean)) + geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
            geom_errorbar(aes(ymin = Min, ymax = Max), position = position_dodge(0.9)) + theme_stata() + labs(y = "Frequency") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        p1 <- ggplot(Qdf_noL, aes(x = Phase, y = Mean)) + geom_bar(stat = "identity", position = position_dodge(), width = 0.8) + 
            geom_errorbar(aes(ymin = Min, ymax = Max), position = position_dodge(0.9)) + theme_stata() + labs(y = "Frequency") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    
    
    
    setwd(curr)
    if(str_length(names) > 50) names <- "All_circuits"
    write_delim(Qdf, paste0(names, "_",zs,"_phase.txt"), delim = "\t")
    
    pdf(paste0(names, "_",zs,"_phase.pdf"))
    print(p)
    dev.off()
    
    pdf(paste0(names, "_",zs,"_phase_noLow.pdf"))
    print(p1)
    dev.off()
}

boolean_converter <- function(folder,n){
    curr <- getwd()
    setwd(folder)
    solution_filez <- list.files(pattern = "nfull", recursive = T)
    result_filez <- solution_filez %>% str_remove("_nfull.txt") %>% paste0("_boolean.txt")
    prs <- read.delim(list.files(pattern = "prs", recursive = T)[1])
    comps <- prs$Parameter[which(str_detect(prs$Parameter, "Prod"))] %>% str_remove("Prod_of_")
    s <- lapply(solution_filez, function(x){#browser()
        sol <- read.delim(x)
        colnames(sol)[1:(ncol(sol)-1)] <- comps
        for (i in comps) 
            sol[[i]] <- ifelse(sol[[i]]>0, 1, 0)
        sol$n_states <- 1/sol$n_states
        sol <- sol %>% group_by_at(1:length(comps)) %>% summarise(Freq = sum(n_states)/n)
        result_file <- x %>% str_remove("_nfull.txt") %>% paste0("_boolean.txt")
        write_delim(sol, path = result_file, delim = "\t")
        sol
    })
    s <- s %>% reduce(merge, by = comps, all = T)
    s[is.na(s)] <- 0
    s$mean <- rowSums(s[, ((ncol(s)-2):ncol(s))])/3
    s <- s[, -((ncol(s)-3):(ncol(s)-1))]
    colnames(s)[ncol(s)] <- folder %>% str_remove("\\./")
    folder <- folder %>% str_remove("\\./") %>% str_remove("/")
    write_delim(s, paste0(folder, "_boolean_m.txt"), delim = "\t")
    setwd(curr)
}


PCA_plots <- function(folder, pattern = "nfull"){
    curr <- getwd()
    setwd(folder)
    
    filez <- list.files(".", pattern = pattern, recursive = T)
    df_list <- lapply(filez, function(f){#browser()
        df <- read.delim(f)
        df_pca <- prcomp(df[, -(ncol(df))])
        pca_comps <- data.frame(t(df_pca$rotation))
        pca_comps$Comps <- colnames(pca_comps)
        pca_comps$SD <- df_pca$sdev
        pca_comps$Percentage = df_pca$sdev*100/sum(df_pca$sdev)
        pca_comps$Cum_per <- cumsum(pca_comps$Percentage)
        pca_comps
    })
    
    lapply(1:length(df_list), function(i){
        df <- df_list[[i]]
        pdf(paste0("PCA_weights_", i, ".pdf"))
        j <- 0
        lapply(data.frame(t(df[, 1:nrow(df)])), function(x){#browser()
            j <<- j + 1
            p <- ggplot(df, aes(x = df$Comps)) + geom_bar( aes_string(y = x), stat = "identity") + theme_stata() + 
                labs(x = "Components", y = "Weight", title = paste0("PC", j, "-" ,round(df$Percentage[j],2),"% SD"))
            print(p)
        })
        dev.off()
        pdf(paste0("PCA_SD_percentage_", i, ".pdf"))
        print(ggplot(rbind.data.frame(0,df), aes(x = c(0,rownames(df)), y = c(0,df$Cum_per), group =1)) + geom_point() + geom_line() + theme_stata() + labs(y = "Percentage deviation", x = "PCA Components"))
        dev.off()
        ## Component wise weight in the PCR's
    })
    
    setwd(curr)
}
