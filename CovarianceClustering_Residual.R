cluster.residvar.fixed <- function(data = NULL, kmeans.clust.output = NULL, data.merge.label = NULL, kmeans.clust.merge.label = "Labels",
                             asreml.workspace = "5gb", asreml.maxit = 35, output.dir="cluster.residvar_outputObjects"){

    ##source the functions needed
    source("icREML.R")
    library("asreml")
    library("tidyverse")

    ##set up some useful variables
    tot.clust <- ncol(kmeans.clust.output$cluster.labels) - 1
    model.list <- list()
    error.i <- 0

    ##################################
    ##SET UP TO BEGIN MODELLING
    ##################################
    for (i in 1:tot.clust){

        ##define cluster name
        clustname <- names(kmeans.clust.output$cluster.labels)[i+1]
        clustnum <- unlist(strsplit(clustname, "k"))[2]
        ##append the cluster being considered to the data
        datawclust <- merge(data, kmeans.clust.output$cluster.labels[,c(kmeans.clust.merge.label, clustname)], by.x = data.merge.label, by.y = kmeans.clust.merge.label)
        ##label this cluster variable as "Cluster" - this is important for asreml later
        datawclust$Cluster <- datawclust[,ncol(datawclust)]
        ##convert the cluster variable to a factor
        datawclust$Cluster <- factor(datawclust$Cluster)

        ##check for singleton clusters (clusters containing single proteins)
        ##if so, remove these clusters from the Cluster variable
        num.prots <- datawclust %>% filter(Order==30) %>% count(Cluster)
        num.single <- 0
        q1.clust <- summary(num.prots$n)[1]
        q2.clust <- summary(num.prots$n)[2]
        q3.clust <- summary(num.prots$n)[3]
        q4.clust <- summary(num.prots$n)[5]
        q5.clust <- summary(num.prots$n)[6]
        if (any(num.prots$n == 1)){
            clust.with.one <- num.prots %>% filter(n==1)
            datawclust$Cluster[datawclust$Cluster %in% clust.with.one$Cluster] <- NA
            datawclust$Cluster <- factor(datawclust$Cluster)
            num.single <- nrow(clust.with.one)
        }

        if (i == 1){
            message(paste0("Current Cluster: k=", clustnum, "; Iteration: i=", i, " of ", tot.clust))
            start.fit <- Sys.time()
            ##################################
            ##SET STARTING VALUES
            ##################################
             mod.kinit.sv <- asreml(log(abundance.scOut+1) ~ -1 + Protein:MaltTime,
                                    random =~ MaltRep + MaltTime:MaltRep +
                                              LabRep + LabRep:LabBlock + Order +
                                              Protein:MaltRep + Protein:MaltTime:MaltRep +
                                              Protein:LabRep + Protein:LabRep:LabBlock +
                                              diag(Cluster):diag(MaltTime):MaltRep:LabRep,
                                    residual =~ diag(Protein):diag(MaltTime):MaltRep:LabRep,
                                    data=datawclust,
                                    workspace= asreml.workspace,
                                    start.values=TRUE,
                                    na.action=na.method(x="include"))

             kinit.sv <- mod.kinit.sv$vparameters.table
             ##Need to fix the second time point variance
             kinit.sv$Value[grep("MaltTime_16.5", kinit.sv$Component)] <- 1
             kinit.sv$Constraint[grep("MaltTime_16.5", kinit.sv$Component)] <- "F"

             ##################################
             ##FIT MODEL
             ##################################
             ##Implementing some error catching
             mod.kinit.asr <- tryCatch(
                 {
                     suppressWarnings(asreml(log(abundance.scOut+1) ~ -1 + Protein:MaltTime,
                                     random =~ MaltRep + MaltTime:MaltRep +
                                               LabRep + LabRep:LabBlock + Order +
                                               Protein:MaltRep + Protein:MaltTime:MaltRep +
                                               Protein:LabRep + Protein:LabRep:LabBlock +
                                               diag(Cluster):diag(MaltTime):MaltRep:LabRep,
                                     residual =~ diag(Protein):diag(MaltTime):MaltRep:LabRep,
                                     sparse =~ Protein:MaltTime,
                                     data=datawclust,
                                     workspace = asreml.workspace,
                                     G.param=kinit.sv, R.param=kinit.sv,
                                     maxit=asreml.maxit,
                                     na.action=na.method(x="include")))
              },
                 error = function(error.message){
                     message(paste0("Error in asreml model, k=",clustnum))
                     message("Here's the error message:")
                     print(error.message)
                     error.i <- i
                     attr(error.i, "class") <- "error"
                     return(error.i)
                 }
                 )
             ##check to see if error.trigger hasn't been triggered
             if(attributes(mod.kinit.asr)$class == "asreml"){
             ##check for convergence and if not, then give more iterations
             if(!mod.kinit.asr$converge){
                 mod.kinit.asr <- tryCatch(
                     {
                         suppressWarnings(update(mod.kinit.asr, maxit=asreml.maxit))
                  },
                     error = function(error.message){
                         message(paste0("Error in asreml model, k=",clustnum))
                         message("Here's the error message:")
                         print(error.message)
                         error.i <- i
                         attr(error.i, "class") <- "error"
                         return(error.i)
             }
                     )
             }
             if(attributes(mod.kinit.asr)$class == "asreml"){
             ##next, test for change in varcomps of greater than 1%
             mod.kinit.asr.varcomps <- summary(mod.kinit.asr)$varcomp
             mod.kinit.asr.perc.change <- mod.kinit.asr.varcomps$`%ch`
             mod.kinit.asr.perc.change <- mod.kinit.asr.perc.change[!is.na(mod.kinit.asr.perc.change)]

             while(any(mod.kinit.asr.perc.change > 1)){
                 mod.kinit.asr <- tryCatch(
                     {
                         suppressWarnings(update(mod.kinit.asr, maxit=asreml.maxit))
                     },
                     error = function(error.message){
                         message(paste0("Error in asreml model, k=",clustnum))
                         message("Here's the error message:")
                         print(error.message)
                         error.i <- i
                         attr(error.i, "class") <- "error"
                         return(error.i)
                     }
                     )
                 if(attributes(mod.kinit.asr)$class == "asreml"){
                 mod.kinit.asr.varcomps <- summary(mod.kinit.asr)$varcomp
                 mod.kinit.asr.perc.change <- mod.kinit.asr.varcomps$`%ch`
                 mod.kinit.asr.perc.change <- mod.kinit.asr.perc.change[!is.na(mod.kinit.asr.perc.change)]
             }
                 else{
                     error.i <- mod.kinit.asr
                     break
                 }
             }
            stop.fit <- Sys.time()

             ##store model
             model.list[[i]] <- mod.kinit.asr
             dir.create(output.dir)
             fname <- file.path(output.dir, paste0("mod",i,".rds"))

             ##use icREML to keep a log of logl, along with AIC, etc.
             asr.opt <- asreml.options()
             asreml.options(Cfixed = TRUE, gammaPar=FALSE)
             message("Commencing icREML update")
             mod.kinit.asr <- tryCatch(
                     {
                         suppressWarnings(update(mod.kinit.asr, maxit=1, sparse=~ . - Protein:MaltTime))
                  },
                     error = function(error.message){
                         message(paste0("Error in asreml model, k=",clustnum))
                         message("Here's the error message:")
                         print(error.message)
                         error.i <- i
                         attr(error.i, "class") <- "error"
                         return(error.i)
             }
                     )

              mod.sum <- tryCatch(
                 {
                     suppressWarnings(icREML(list(mod.kinit.asr)))
                 },
                 error = function(error.message){
                     message(paste0("Error when running icREML, k=",clustnum))
                     message("Here's the error message:")
                     print(error.message)
                     temp.icREML <- data.frame('model'=paste0("k=",i),
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
             asreml.options(asr.opt)
             mod.sum$model <- paste0("k=",clustnum)
             mod.sum$num.singleton <- num.single
             mod.sum$q1.clust <- q1.clust
             mod.sum$q2.clust <- q2.clust
             mod.sum$q3.clust <- q3.clust
             mod.sum$q4.clust <- q4.clust
             mod.sum$q5.clust <- q5.clust
             mod.sum$fit.time <- stop.fit - start.fit
             mod.list <- mod.sum
             print(mod.list)
             fname <- file.path(output.dir, "ModelList_icREML.csv")
             write.csv(mod.list, fname, row.names=FALSE)
         }
             else{
                 error.i <- mod.kinit.asr
             }
        }
            else{
                error.i <- mod.kinit.asr
            }
        }
        else {
            message(paste0("Current Cluster: k=", clustnum, "; Iteration: i=", i, " of ", tot.clust))
            ##################################
            ##SET STARTING VALUES
            ##################################
             start.fit <- Sys.time()
             mod.sv <- asreml(log(abundance.scOut+1) ~ -1 + Protein:MaltTime,
                                    random =~ MaltRep + MaltTime:MaltRep +
                                              LabRep + LabRep:LabBlock + Order +
                                              Protein:MaltRep + Protein:MaltTime:MaltRep +
                                              Protein:LabRep + Protein:LabRep:LabBlock +
                                              diag(Cluster):diag(MaltTime):MaltRep:LabRep,
                                    residual =~ diag(Protein):diag(MaltTime):MaltRep:LabRep,
                                    data=datawclust,
                                    workspace= asreml.workspace,
                                    start.values=TRUE,
                                    na.action=na.method(x="include"))

             sv <- mod.sv$vparameters.table
             ##Need to fix the second time point variance
             sv$Value[grep("MaltTime_16.5", sv$Component)] <- 1
             sv$Constraint[grep("MaltTime_16.5", sv$Component)] <- "F"

             ##################################
             ##FIT MODEL
             ##################################
             mod.asr <- tryCatch(
                 {
                     suppressWarnings(asreml(log(abundance.scOut+1) ~ -1 + Protein:MaltTime,
                                     random =~ MaltRep + MaltTime:MaltRep +
                                               LabRep + LabRep:LabBlock + Order +
                                               Protein:MaltRep + Protein:MaltTime:MaltRep +
                                               Protein:LabRep + Protein:LabRep:LabBlock +
                                               diag(Cluster):diag(MaltTime):MaltRep:LabRep,
                                     residual =~ diag(Protein):diag(MaltTime):MaltRep:LabRep,
                                     sparse=~ Protein:MaltTime,
                                     data=datawclust,
                                     workspace = asreml.workspace,
                                     G.param=sv, R.param=sv,
                                     maxit=asreml.maxit,
                                     na.action=na.method(x="include")))
              },
                 error = function(error.message){
                     message(paste0("Error in asreml model, k=",clustnum))
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
                         message(paste0("Error in asreml model, k=",clustnum))
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
                         message(paste0("Error in asreml model, k=",clustnum))
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
             asr.opt <- asreml.options()
             asreml.options(Cfixed = TRUE, gammaPar=FALSE)
             message("Commencing icREML update")
             mod.asr <- tryCatch(
                     {
                         suppressWarnings(update(mod.asr, maxit=1, sparse=~ . - Protein:MaltTime))
                  },
                     error = function(error.message){
                         #error.trigger <<- 1
                         #error.i <<- i
                         message(paste0("Error in asreml model, k=",clustnum))
                         message("Here's the error message:")
                         print(error.message)
                         error.i <- i
                         attr(error.i, "class") <- "error"
                         return(error.i)
             }
                     )
             mod.sum <- tryCatch(
                 {
                     suppressWarnings(icREML(list(mod.asr)))
                 },
                 error = function(error.message){
                     message(paste0("Error when running icREML, k=",clustnum))
                     message("Here's the error message:")
                     print(error.message)
                     temp.icREML <- data.frame('model'=paste0("k=",clustnum),
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
             asreml.options(asr.opt)
             mod.sum$model <- paste0("k=",clustnum)
             mod.sum$num.singleton <- num.single
             mod.sum$q1.clust <- q1.clust
             mod.sum$q2.clust <- q2.clust
             mod.sum$q3.clust <- q3.clust
             mod.sum$q4.clust <- q4.clust
             mod.sum$q5.clust <- q5.clust
             mod.sum$fit.time <- stop.fit - start.fit
             mod.list <- rbind(mod.list, mod.sum)
             print(mod.list)
             fname <- file.path(output.dir, "ModelList_icREML.csv")
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
    }

         ##prepare output
         return.objects <- list("icREML.list" = mod.list,
                                "model.list" = model.list)
         return(return.objects)

}
