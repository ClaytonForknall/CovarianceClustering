kmeans.clust <- function(input.mat=NULL, tot.clust=NULL, seed=5995, cluster.sequence="narrow", growth.rate=0.2){
    #set up a list to store clustering output
    temp.clusters <- list()

    if (cluster.sequence == "narrow"){
    ##begin loop that will perform 2:tot.clust iterations
    for (i in 2:tot.clust){
        ##set the seed
        set.seed(seed)
        ##perform clustering
        temp.kmeans <- kmeans(input.mat, i, iter.max=25)
        ##save a copy of the output
        temp.clusters[[i]] <- temp.kmeans
    }

    ##now want to pull out the clusters and store them in a dataframe
    cluster.labels <- data.frame("Labels" = names(temp.clusters[[2]]$cluster))
    for (i in 2:tot.clust){
        cluster.labels$temp <- temp.clusters[[i]]$cluster
        names(cluster.labels)[ncol(cluster.labels)] <- paste0("Clusterk",i)
    }
}
    else if (cluster.sequence == "broad"){
        ##generate a range of clusters using a geometric growth of growth rate each iteration from 2 to tot.clust
        clust.seq <- list()
        current.clust <- 2
        j <- 1

        while(j <= tot.clust & current.clust < 617){
            clust.seq[[j]] <- current.clust
            current.clust <- ceiling((1 + growth.rate)*current.clust)
            j <- j+1
        }

        ##provide this list of clusters to the loop defined for the narrow search above
        clust.seq <- unlist(clust.seq)
        ##begin loop that will perform 2:length(clust.seq) iterations
        for (i in 1:length(clust.seq)){
            ##set the seed
            set.seed(seed)
            ##perform clustering
            temp.kmeans <- kmeans(input.mat, clust.seq[i], iter.max=25)
            ##save a copy of the output
            temp.clusters[[paste0(clust.seq[i])]] <- temp.kmeans
        }

        ##now want to pull out the clusters and store them in a dataframe
        cluster.labels <- data.frame("Labels" = names(temp.clusters[["2"]]$cluster))
        for (i in 1:length(clust.seq)){
            cluster.labels$temp <- temp.clusters[[i]]$cluster
            names(cluster.labels)[ncol(cluster.labels)] <- paste0("Clusterk",clust.seq[i])
        }
    }
    else {
        stop("Please specify cluster.sequence to be either 'narrow' or 'broad'", call.=FALSE)
    }

    ##prepare output
    return.objects <- list("cluster.labels" = cluster.labels,
                           "kmeans.output.list" = temp.clusters)
    return(return.objects)
}
