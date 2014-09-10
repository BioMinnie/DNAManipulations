#
#    Script: motifOccurrence.R
#    Author: Benjamin Tovar
#    Date: 21/April/2012
#
################################################################################

#                       ############
#                        FUNCTIONS:
#                       ############

##############################################################################
iterateComputeDistance <- function(multinomialDNAmodel,
                                   lengthDNAseq,
                                   motif,
                                   numberOfIterations){

    # This function returns the mean distance 
    # of a given motif given X number of DNA sequences given a multinomial model
    # (probability distribution of each base).
    
    # So, it will generate X number of DNA sequences using a given 
    # probability distribution and then it will compute the distance among 
    # that mofit within the total set of sequences to finally returns 
    # the average distance of the motif.
    
    result <- rep(NA,numberOfIterations)
    for(i in 1:numberOfIterations){
        currentGenome <- NA
        currentCoordinatesOfMotif <- NA        
        currentSequence <- sample(c("A","C","G","T"),
                                    lengthDNAseq,rep=T,
                                    prob=multinomialDNAmodel)
        currentCoordinatesOfMotif <- coordMotif(currentSequence,motif)
        result[i] <- computeDistance(currentCoordinatesOfMotif)
        cat(" *** Iteration number: ",i," completed *** | average distance = "
            ,result[i],"\n")  
    }
    result <- trunc(mean(result))
    cat(" \n*** Computation status: DONE ***\n\n")
    return(result)  
}
##############################################################################
coordMotif <- function(targetSequence,motif){

    # This function returns the coordinates of the motif of study in a target 
    # DNA sequence. In other words, if I found the motif, tell me exactly in
    # which position of the DNA sequence is.
    
    lengthMotif <- length(motif)
    lengthTargetSeq <- (length(targetSequence)-lengthMotif)
    motif <- toString(motif)
    motif <- gsub(", ","",motif)
    res <- 1    
    for(i in 1:lengthTargetSeq){
        currentTargetSeq <- targetSequence[i:(i+(lengthMotif)-1)]
        currentTargetSeq <- toString(currentTargetSeq)
        currentTargetSeq <- gsub(", ","",currentTargetSeq)
        if(currentTargetSeq == motif){
            res[(length(res)+1)] <- i
        }
    }
    return(res)
}
##############################################################################
computeDistance <- function(coordinatesOfMotif){
    
    # This function returns the mean distance 
    # of a given motif given its coordinates within a target DNA sequence.
    # In other words, If I already got a list with the coordinates where the 
    # motif is inside a DNA sequence, tell me the average distance between
    # this coordinates to get the expected distance of that motif.

    currentDistance <- rep(NA,(length(coordinatesOfMotif)-1))
    lengthCoord <- length(currentDistance)
    for(i in 1:lengthCoord){
        currentDistance[i] <- coordinatesOfMotif[i+1]-coordinatesOfMotif[i]
    }
    res <- trunc(mean(currentDistance))
    return(res)   
}
##############################################################################
computeExpectedDistance <- function(multinomialModel,
                                    lengthDNAseq,
                                    motif){
                                    
    # This function computes the expected distance of a given motif in a DNA
    # sequence given its multinomial model (probability distribution of 
    # each base)
    
    # Convert the motif into an index                                       
    motifIndex <- gsub("A",1,motif); motifIndex <- gsub("C",2,motifIndex)
    motifIndex <- gsub("G",3,motifIndex); motifIndex <- gsub("T",4,motifIndex)
    motifIndex <- as.numeric(motifIndex)
    # Compute p value of the motif given the multinomial model
    p <- rep(NA,length(motif))
    for(i in 1:length(motifIndex)){
        p[i] <- multinomialModel[motifIndex[i]]
    }    
    p <- prod(p)
    result <- trunc(lengthDNAseq/(lengthDNAseq*p))
    return(result)
}                                   
##############################################################################

# Benjamin
