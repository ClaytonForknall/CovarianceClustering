cluster.response.random <- function(data = NULL, kmeans.clust.output.rand = NULL,
                                     data.merge.label = NULL, kmeans.clust.merge.label = "Labels",
                                     asreml.workspace = "15gb", asreml.maxit = 35, output.dir="cluster.responses.rand_outputObjects"){
    ##set up the functions needed
    source("icREML.R")
    library("tidyverse")
    library("asreml")

    ##set up some useful variables and initialise the storage directory
    tot.clust.rand <- ncol(kmeans.clust.output.rand$cluster.labels) - 1

    model.list <- list()
    error.i <- 0
    mod.list <- NULL
    dir.create(output.dir)

    ##set up clust.list.df
    clust.list.df <- data.frame(rand.clust = as.numeric(names(kmeans.clust.output.rand$kmeans.output.list)))
    clust.list.df$unique.clust <- 1:tot.clust.rand

    ##################################
    ##SET UP TO BEGIN MODELLING
    ##################################
    ##Set asreml.options to make icREML easier to use later
    asreml.options(Cfixed = TRUE, gammaPar = FALSE)

    for (i in 1:tot.clust.rand){

        ##define unique cluster name
        clustname <- paste0("k",clust.list.df$rand.clust[i])
        clustnum <- i

        rand.clustnum <- clust.list.df$rand.clust[i]
        rand.clustname <- paste0("Clusterk",rand.clustnum)

        ##append the cluster being considered to the data
        datawclust <- merge(data, kmeans.clust.output.rand$cluster.labels[,c(kmeans.clust.merge.label, rand.clustname)], by.x = data.merge.label, by.y = kmeans.clust.merge.label)
        ##label this cluster variable as "RClust" - this is important for asreml later
        datawclust$RClust <- datawclust[,ncol(datawclust)]
        ##convert the RClust variable to a factor
        datawclust$RClust <- factor(datawclust$RClust)

        ##check for singleton clusters (clusters containing single proteins)
        ##if so, remove these clusters from the cluster variable
        num.prots.rand <- datawclust %>% filter(MaltTime == 0) %>% count(RClust)
        num.single.rand <- 0
        q1.clust.r <- summary(num.prots.rand$n)[1]
        q2.clust.r <- summary(num.prots.rand$n)[2]
        q3.clust.r <- summary(num.prots.rand$n)[3]
        q4.clust.r <- summary(num.prots.rand$n)[5]
        q5.clust.r <- summary(num.prots.rand$n)[6]
        if (any(num.prots.rand$n == 1)){
            clust.with.one.r <- num.prots.rand %>% filter(n==1)
            datawclust$RClust[datawclust$RClust %in% clust.with.one.r$RClust] <- NA
            datawclust$RClust <- factor(datawclust$RClust)
            num.single.rand <- nrow(clust.with.one.r)
        }

            message(paste0("Current Cluster: ", clustname, "; Iteration: i=", i, " of ", tot.clust.rand))
            ##################################
            ##SET STARTING VALUES
            ##################################
             start.fit <- Sys.time()

             ##################################
             ##FIT MODEL
             ##################################
             mod.asr <- tryCatch(
                 {
                     suppressWarnings(asreml(solution ~ Protein + malttime + Protein:malttime,
                                     random =~ spl(malttime) + RClust:spl(malttime) +
                                               Protein:spl(malttime) + MaltTime + Protein:MaltTime,
                                     family = asr_gaussian(dispersion = 1),
                                     weights = weights,
                                     data=datawclust,
                                     workspace = asreml.workspace,
                                     maxit=asreml.maxit,
                                     na.action=na.method(x="include")))
              },
                 error = function(error.message){
                     message(paste0("Error in asreml model, mod=",clustname))
                     message("Here's the error message:")
                     print(error.message)
                     error.i <- i
                     attr(error.i, "class") <- "error"
                     return(error.i)
                 }
             )
            ##check the length of the mod.asr object
            if(attributes(mod.asr)$class == "asreml"){
             ##check for convergence and if not, then give more iterations
             if(!mod.asr$converge){
                 mod.asr <- tryCatch(
                     {
                         suppressWarnings(update(mod.asr, maxit=asreml.maxit))
                     },
                     error = function(error.message){
                         message(paste0("Error in asreml model, mod=",clustname))
                         message("Here's the error message:")
                         print(error.message)
                         error.i <- i
                         attr(error.i, "class") <- "error"
                         return(error.i)
                     }
                     )
             }

             if(attributes(mod.asr)$class == "asreml"){
             ##next, test for change in varcomps of greater than 1%
             mod.asr.varcomps <- summary(mod.asr)$varcomp
             mod.asr.perc.change <- mod.asr.varcomps$`%ch`
             mod.asr.perc.change <- mod.asr.perc.change[!is.na(mod.asr.perc.change)]

             while(any(mod.asr.perc.change > 1)){
                 mod.asr <- tryCatch(
                     {
                         suppressWarnings(update(mod.asr, maxit=asreml.maxit))
                     },
                     error = function(error.message){
                         message(paste0("Error in asreml model, mod=",clustname))
                         message("Here's the error message:")
                         print(error.message)
                         error.i <- i
                         attr(error.i, "class") <- "error"
                         return(error.i)
                     }
                     )
                 if(attributes(mod.asr)$class == "asreml"){
                 mod.asr.varcomps <- summary(mod.asr)$varcomp
                 mod.asr.perc.change <- mod.asr.varcomps$`%ch`
                 mod.asr.perc.change <- mod.asr.perc.change[!is.na(mod.asr.perc.change)]
             }
                 else{
                     error.i <- mod.asr
                     break
                 }
             }
             stop.fit <- Sys.time()

             ##store model
             model.list[[i]] <- mod.asr

             ##use icREML to keep a log of logl, along with AIC, etc.
             mod.sum <- tryCatch(
                 {
                     suppressWarnings(icREML(list(mod.asr)))
                 },
                 error = function(error.message){
                     message(paste0("Error when running icREML, mod=",clustname))
                     message("Here's the error message:")
                     print(error.message)
                     temp.icREML <- data.frame('model'=clustname,
                                               'rclust'=rand.clustnum,
                                               'unique.mod'=clustnum,
                                               'res.loglik'=NA,
                                               'full.loglik'=NA,
                                               'p'=NA,
                                               'q'=NA,
                                               'b'=NA,
                                               'AIC'=NA,
                                               'BIC'=NA)
                     return(temp.icREML)
                 }
                 )
             mod.sum$model <- clustname
             mod.sum$rclust <- rand.clustnum
             mod.sum$num.singleton.r <- num.single.rand
             mod.sum$q1.clust.r <- q1.clust.r
             mod.sum$q2.clust.r <- q2.clust.r
             mod.sum$q3.clust.r <- q3.clust.r
             mod.sum$q4.clust.r <- q4.clust.r
             mod.sum$q5.clust.r <- q5.clust.r
             mod.sum$fit.time <- stop.fit - start.fit
             mod.list <- rbind(mod.list, mod.sum)
             print(mod.list)
             fname <- file.path(output.dir, "ModelList_RandomResponses_icREML.csv")
             write.csv(mod.list, fname, row.names=FALSE)
        }
             else{
                 error.i <- mod.asr
             }
         }
            else{
                error.i <- mod.asr
            }
    }

         ##prepare output
         return.objects <- list("icREML.list" = mod.list,
                                "model.list" = model.list)
         return(return.objects)

}
