#########################################################################################
# License Agreement
#
# THIS WORK IS PROVIDED UNDER THE TERMS OF THIS CREATIVE COMMONS PUBLIC LICENSE
# ("CCPL" OR "LICENSE"). THE WORK IS PROTECTED BY COPYRIGHT AND/OR OTHER
# APPLICABLE LAW. ANY USE OF THE WORK OTHER THAN AS AUTHORIZED UNDER THIS LICENSE
# OR COPYRIGHT LAW IS PROHIBITED.
#
# BY EXERCISING ANY RIGHTS TO THE WORK PROVIDED HERE, YOU ACCEPT AND AGREE TO BE
# BOUND BY THE TERMS OF THIS LICENSE. TO THE EXTENT THIS LICENSE MAY BE CONSIDERED
# TO BE A CONTRACT, THE LICENSOR GRANTS YOU THE RIGHTS CONTAINED HERE IN
# CONSIDERATION OF YOUR ACCEPTANCE OF SUCH TERMS AND CONDITIONS.
#
# BASELIne: Bayesian Estimation of Antigen-Driven Selection in Immunoglobulin Sequences
# Coded by: Mohamed Uduman & Gur Yaari
# Copyright 2012 Kleinstein Lab
# Version: 1.3 (01/23/2014)
#########################################################################################

op <- options();
options(showWarnCalls=FALSE, showErrorCalls=FALSE, warn=-1)
library('seqinr')
#if( Sys.info()[1]=="Linux"){
#  library("multicore")
#}


# Load functions and initialize global variables
source("Baseline_Functions_Version1.3.r")
# Initialize parameters with user provided arguments
  arg <- commandArgs(TRUE)
  #arg = c(2,1,5,5,0,1,"1:26:38:55:65:104:116", "test.fasta","","sample")
  #arg = c(1,1,5,5,0,1,"1:38:55:65:104:116:200", "test.fasta","","sample")
  #arg = c(1,1,5,5,1,1,"1:26:38:55:65:104:116", "/home/mu37/Wu/Wu_Cloned_gapped_sequences_D-masked.fasta","/home/mu37/Wu/","Wu")
  testID <- as.numeric(arg[1])                    # 1 = Focused, 2 = Local
  species <- as.numeric(arg[2])                   # 1 = Human. 2 = Mouse
  substitutionModel <- as.numeric(arg[3])         # 0 = Uniform substitution, 1 = Smith DS et al. 1996, 5 = FiveS
  mutabilityModel <- as.numeric(arg[4])           # 0 = Uniform mutablity, 1 = Tri-nucleotide (Shapiro GS et al. 2002)  , 5 = FiveS
  clonal <- as.numeric(arg[5])                    # 0 = Independent sequences, 1 = Clonally related, 2 = Clonally related & only non-terminal mutations
  fixIndels <- as.numeric(arg[6])                 # 0 = Do nothing, 1 = Try and fix Indels
  region <- as.numeric(strsplit(arg[7],":")[[1]]) # StartPos:LastNucleotideF1:C1:F2:C2:F3:C3
  inputFilePath <- arg[8]                         # Full path to input file
  outputPath <- arg[9]                            # Full path to location of output files
  outputID <- arg[10]                             # ID for session output

  if(testID==5){
    traitChangeModel <- 1
    if( !is.na(any(arg[11])) ) traitChangeModel <- as.numeric(arg[11])    # 1 <- Chothia 1998
    initializeTraitChange(traitChangeModel)
  }

# Initialize other parameters/variables

  # Initialzie the codon table ( definitions of R/S )
  computeCodonTable(testID)

  # Initialize
  # Test Name
  testName<-"Focused"
  if(testID==2) testName<-"Local"
  if(testID==3) testName<-"Imbalanced"
  if(testID==4) testName<-"ImbalancedSilent"

  # Indel placeholders initialization
  indelPos <- NULL
  delPos <- NULL
  insPos <- NULL

  # Initialize in Tranistion & Mutability matrixes
  substitution <- initializeSubstitutionMatrix(substitutionModel,species)
  mutability <- initializeMutabilityMatrix(mutabilityModel,species)

  # FWR/CDR boundaries
  flagTrim <- F
  if( is.na(region[7])){
    flagTrim <- T
    region[7]<-region[6]
  }
  readStart = min(region,na.rm=T)
  readEnd = max(region,na.rm=T)
  if(readStart>1){
    region = region - (readStart - 1)
  }
  region_Nuc = c( (region[1]*3-2) , (region[2:7]*3) )
  region_Cod = region

  readStart = (readStart*3)-2
  readEnd = (readEnd*3)

    FWR_Nuc <- c( rep(TRUE,(region_Nuc[2])),
                  rep(FALSE,(region_Nuc[3]-region_Nuc[2])),
                  rep(TRUE,(region_Nuc[4]-region_Nuc[3])),
                  rep(FALSE,(region_Nuc[5]-region_Nuc[4])),
                  rep(TRUE,(region_Nuc[6]-region_Nuc[5])),
                  rep(FALSE,(region_Nuc[7]-region_Nuc[6]))
                )
    CDR_Nuc <- (1-FWR_Nuc)
    CDR_Nuc <- as.logical(CDR_Nuc)
    FWR_Nuc_Mat <- matrix( rep(FWR_Nuc,4), ncol=length(FWR_Nuc), nrow=4, byrow=T)
    CDR_Nuc_Mat <- matrix( rep(CDR_Nuc,4), ncol=length(CDR_Nuc), nrow=4, byrow=T)

    FWR_Codon <- c( rep(TRUE,(region[2])),
                  rep(FALSE,(region[3]-region[2])),
                  rep(TRUE,(region[4]-region[3])),
                  rep(FALSE,(region[5]-region[4])),
                  rep(TRUE,(region[6]-region[5])),
                  rep(FALSE,(region[7]-region[6]))
                )
    CDR_Codon <- (1-FWR_Codon)
    CDR_Codon <- as.logical(CDR_Codon)


# Read input FASTA file
  tryCatch(
    inputFASTA <- baseline.read.fasta(inputFilePath, seqtype="DNA",as.string=T,set.attributes=F,forceDNAtolower=F)
    , error = function(ex){
      cat("Error|Error reading input. Please enter or upload a valid FASTA file.\n")
      q()
    }
  )

  if (length(inputFASTA)==1) {
    cat("Error|Error reading input. Please enter or upload a valid FASTA file.\n")
    q()
  }

  # Process sequence IDs/names
  names(inputFASTA) <- sapply(names(inputFASTA),function(x){trim(x)})

  # Convert non nucleotide characters to N
  inputFASTA[length(inputFASTA)] = gsub("\t","",inputFASTA[length(inputFASTA)])
  inputFASTA <- lapply(inputFASTA,replaceNonFASTAChars)

  # Process the FASTA file and conver to Matrix[inputSequence, germlineSequence]
  processedInput <- processInputAdvanced(inputFASTA)
  matInput <- processedInput[[1]]
  germlines <- processedInput[[2]]
  lenGermlines = length(unique(germlines))
  groups <- processedInput[[3]]
  lenGroups = length(unique(groups))
  rm(processedInput)
  rm(inputFASTA)

#   # remove clones with less than 2 seqeunces
#   tableGL <- table(germlines)
#   singletons <- which(tableGL<8)
#   rowsToRemove <- match(singletons,germlines)
#   if(any(rowsToRemove)){
#     matInput <- matInput[-rowsToRemove,]
#     germlines <- germlines[-rowsToRemove]
#     groups <- groups[-rowsToRemove]
#   }
#
#   # remove unproductive seqs
#   nonFuctionalSeqs <- sapply(rownames(matInput),function(x){any(grep("unproductive",x))})
#   if(any(nonFuctionalSeqs)){
#     if(sum(nonFuctionalSeqs)==length(germlines)){
#       write.table("Unproductive",file=paste(outputPath,outputID,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
#       q()
#     }
#     matInput <- matInput[-which(nonFuctionalSeqs),]
#     germlines <- germlines[-which(nonFuctionalSeqs)]
#     germlines[1:length(germlines)] <- 1:length(germlines)
#     groups <- groups[-which(nonFuctionalSeqs)]
#   }
#
#   if(class(matInput)=="character"){
#     write.table("All unproductive seqs",file=paste(outputPath,outputID,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
#     q()
#   }
#
#   if(nrow(matInput)<10 | is.null(nrow(matInput))){
#     write.table(paste(nrow(matInput), "seqs only",sep=""),file=paste(outputPath,outputID,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
#     q()
#   }

# replace leading & trailing "-" with "N:
  matInput <- t(apply(matInput,1,replaceLeadingTrailingDashes,readEnd))

  # Trim (nucleotide) input sequences to the last codon
  #matInput[,1] <- apply(matrix(matInput[,1]),1,trimToLastCodon)

#   # Check for Indels
#   if(fixIndels){
#     delPos <- fixDeletions(matInput)
#     insPos <- fixInsertions(matInput)
#   }else{
#     # Check for indels
#     indelPos <- checkForInDels(matInput)
#     indelPos <- apply(cbind(indelPos[[1]],indelPos[[2]]),1,function(x){(x[1]==T & x[2]==T)})
#   }

  # If indels are present, remove mutations in the seqeunce & throw warning at end
  #matInput[indelPos,] <- apply(matrix(matInput[indelPos,],nrow=sum(indelPos),ncol=2),1,function(x){x[1]=x[2]; return(x) })

  colnames(matInput)=c("Input","Germline")

  # If seqeunces are clonal, create effective sequence for each clone & modify germline/group definitions
  germlinesOriginal = NULL
  if(clonal){
    germlinesOriginal <- germlines
    collapseCloneResults <- tapply(1:nrow(matInput),germlines,function(i){
                                                                collapseClone(matInput[i,1],matInput[i[1],2],readEnd,nonTerminalOnly=(clonal-1))
                                                              })
    matInput = t(sapply(collapseCloneResults,function(x){return(x[[1]])}))
    names_groups = tapply(groups,germlines,function(x){names(x[1])})
    groups = tapply(groups,germlines,function(x){array(x[1],dimnames=names(x[1]))})
    names(groups) = names_groups

    names_germlines =  tapply(germlines,germlines,function(x){names(x[1])})
    germlines = tapply(   germlines,germlines,function(x){array(x[1],dimnames=names(x[1]))}   )
    names(germlines) = names_germlines
    matInputErrors = sapply(collapseCloneResults,function(x){return(x[[2]])})
  }


# Selection Analysis


#  if (length(germlines)>sequenceLimit) {
#    # Code to parallelize processing goes here
#    stop( paste("Error: Cannot process more than ", Upper_limit," sequences",sep="") )
#  }

#  if (length(germlines)<sequenceLimit) {}

    # Compute expected mutation frequencies
    matExpected <- getExpectedIndividual(matInput)

    # Count observed number of mutations in the different regions
    mutations <- lapply( 1:nrow(matInput),  function(i){
                                              #cat(i,"\n")
                                              seqI = s2c(matInput[i,1])
                                              seqG = s2c(matInput[i,2])
                                              matIGL = matrix(c(seqI,seqG),ncol=length(seqI),nrow=2,byrow=T)
                                              retVal <- NA
                                              tryCatch(
                                                retVal <- analyzeMutations2NucUri(matIGL)
                                                , error = function(ex){
                                                  retVal <- NA
                                                }
                                              )


                                              return( retVal )
                                            })

    matObserved <- t(sapply( mutations, processNucMutations2 ))
    numberOfSeqsWithMutations <- numberOfSeqsWithMutations(matObserved, testID)

    #if(sum(numberOfSeqsWithMutations)==0){
    #  write.table("No mutated sequences",file=paste(outputPath,outputID,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
    #  q()
    #}

    matMutationInfo <- cbind(matObserved,matExpected)
    rm(matObserved,matExpected)


    #Bayesian  PDFs
    bayes_pdf = computeBayesianScore(matMutationInfo, test=testName, max_sigma=20,length_sigma=4001)
    bayesPDF_cdr = bayes_pdf[[1]]
    bayesPDF_fwr = bayes_pdf[[2]]
    rm(bayes_pdf)

    bayesPDF_germlines_cdr = tapply(bayesPDF_cdr,germlines,function(x) groupPosteriors(x,length_sigma=4001))
    bayesPDF_germlines_fwr = tapply(bayesPDF_fwr,germlines,function(x) groupPosteriors(x,length_sigma=4001))

    bayesPDF_groups_cdr = tapply(bayesPDF_cdr,groups,function(x) groupPosteriors(x,length_sigma=4001))
    bayesPDF_groups_fwr = tapply(bayesPDF_fwr,groups,function(x) groupPosteriors(x,length_sigma=4001))

    if(lenGroups>1){
      groups <- c(groups,lenGroups+1)
      names(groups)[length(groups)] = "All sequences combined"
      bayesPDF_groups_cdr[[lenGroups+1]] =   groupPosteriors(bayesPDF_groups_cdr,length_sigma=4001)
      bayesPDF_groups_fwr[[lenGroups+1]] =   groupPosteriors(bayesPDF_groups_fwr,length_sigma=4001)
    }

    #Bayesian  Outputs
    bayes_cdr =  t(sapply(bayesPDF_cdr,calcBayesOutputInfo))
    bayes_fwr =  t(sapply(bayesPDF_fwr,calcBayesOutputInfo))
    bayes_germlines_cdr =  t(sapply(bayesPDF_germlines_cdr,calcBayesOutputInfo))
    bayes_germlines_fwr =  t(sapply(bayesPDF_germlines_fwr,calcBayesOutputInfo))
    bayes_groups_cdr =  t(sapply(bayesPDF_groups_cdr,calcBayesOutputInfo))
    bayes_groups_fwr =  t(sapply(bayesPDF_groups_fwr,calcBayesOutputInfo))

    #P-values
    simgaP_cdr = sapply(bayesPDF_cdr,computeSigmaP)
    simgaP_fwr = sapply(bayesPDF_fwr,computeSigmaP)

    simgaP_germlines_cdr = sapply(bayesPDF_germlines_cdr,computeSigmaP)
    simgaP_germlines_fwr = sapply(bayesPDF_germlines_fwr,computeSigmaP)

    simgaP_groups_cdr = sapply(bayesPDF_groups_cdr,computeSigmaP)
    simgaP_groups_fwr = sapply(bayesPDF_groups_fwr,computeSigmaP)


    #Format output

    # Round expected mutation frequencies to 3 decimal places
    matMutationInfo[germlinesOriginal[indelPos],] = NA
    if(nrow(matMutationInfo)==1){
      matMutationInfo[5:8] = round(matMutationInfo[,5:8]/sum(matMutationInfo[,5:8],na.rm=T),3)
    }else{
      matMutationInfo[,5:8] = t(round(apply(matMutationInfo[,5:8],1,function(x){ return(x/sum(x,na.rm=T)) }),3))
    }

    listPDFs = list()
    nRows = length(unique(groups)) + length(unique(germlines)) + length(groups)

    matOutput = matrix(NA,ncol=18,nrow=nRows)
    rowNumb = 1
    for(G in unique(groups)){
      #print(G)
      matOutput[rowNumb,c(1,2,11:18)] = c("Group",names(groups)[groups==G][1],bayes_groups_cdr[G,],bayes_groups_fwr[G,],simgaP_groups_cdr[G],simgaP_groups_fwr[G])
      listPDFs[[rowNumb]] = list("CDR"=bayesPDF_groups_cdr[[G]],"FWR"=bayesPDF_groups_fwr[[G]])
      names(listPDFs)[rowNumb] = names(groups[groups==paste(G)])[1]
      #if(names(groups)[which(groups==G)[1]]!="All sequences combined"){
      gs = unique(germlines[groups==G])
      rowNumb = rowNumb+1
      if( !is.na(gs) ){
        for( g in gs ){
          matOutput[rowNumb,c(1,2,11:18)] = c("Germline",names(germlines)[germlines==g][1],bayes_germlines_cdr[g,],bayes_germlines_fwr[g,],simgaP_germlines_cdr[g],simgaP_germlines_fwr[g])
          listPDFs[[rowNumb]] = list("CDR"=bayesPDF_germlines_cdr[[g]],"FWR"=bayesPDF_germlines_fwr[[g]])
          names(listPDFs)[rowNumb] = names(germlines[germlines==paste(g)])[1]
          rowNumb = rowNumb+1
          indexesOfInterest = which(germlines==g)
          numbSeqsOfInterest =  length(indexesOfInterest)
          rowNumb = seq(rowNumb,rowNumb+(numbSeqsOfInterest-1))
          matOutput[rowNumb,] = matrix(   c(  rep("Sequence",numbSeqsOfInterest),
                                              rownames(matInput)[indexesOfInterest],
                                              c(matMutationInfo[indexesOfInterest,1:4]),
                                              c(matMutationInfo[indexesOfInterest,5:8]),
                                              c(bayes_cdr[indexesOfInterest,]),
                                              c(bayes_fwr[indexesOfInterest,]),
                                              c(simgaP_cdr[indexesOfInterest]),
                                              c(simgaP_fwr[indexesOfInterest])
          ), ncol=18, nrow=numbSeqsOfInterest,byrow=F)
          increment=0
          for( ioi in indexesOfInterest){
            listPDFs[[min(rowNumb)+increment]] =  list("CDR"=bayesPDF_cdr[[ioi]] , "FWR"=bayesPDF_fwr[[ioi]])
            names(listPDFs)[min(rowNumb)+increment] = rownames(matInput)[ioi]
            increment = increment + 1
          }
          rowNumb=max(rowNumb)+1

        }
      }
    }
    colsToFormat = 11:18
    matOutput[,colsToFormat] = formatC(  matrix(as.numeric(matOutput[,colsToFormat]), nrow=nrow(matOutput), ncol=length(colsToFormat)) ,  digits=3)
    matOutput[matOutput== " NaN"] = NA



    colnames(matOutput) = c("Type", "ID", "Observed_CDR_R", "Observed_CDR_S", "Observed_FWR_R", "Observed_FWR_S",
                            "Expected_CDR_R", "Expected_CDR_S", "Expected_FWR_R", "Expected_FWR_S",
                            paste( rep(testName,6), rep(c("Sigma","CIlower","CIupper"),2),rep(c("CDR","FWR"),each=3), sep="_"),
                            paste( rep(testName,2), rep("P",2),c("CDR","FWR"), sep="_")
    )
    fileName = paste(outputPath,outputID,".txt",sep="")
    write.table(matOutput,file=fileName,quote=F,sep="\t",row.names=F,col.names=T)
    fileName = paste(outputPath,outputID,".RData",sep="")
    save(listPDFs,file=fileName)

indelWarning = FALSE
if(sum(indelPos)>0){
  indelWarning = "<P>Warning: The following sequences have either gaps and/or deletions, and have been ommited from the analysis.";
  indelWarning = paste( indelWarning , "<UL>", sep="" )
  for(indels in names(indelPos)[indelPos]){
    indelWarning = paste( indelWarning , "<LI>", indels, "</LI>", sep="" )
  }
  indelWarning = paste( indelWarning , "</UL></P>", sep="" )
}

cloneWarning = FALSE
if(clonal==1){
  if(sum(matInputErrors)>0){
    cloneWarning = "<P>Warning: The following clones have sequences of unequal length.";
    cloneWarning = paste( cloneWarning , "<UL>", sep="" )
    for(clone in names(matInputErrors)[matInputErrors]){
      cloneWarning = paste( cloneWarning , "<LI>", names(germlines)[as.numeric(clone)], "</LI>", sep="" )
    }
    cloneWarning = paste( cloneWarning , "</UL></P>", sep="" )
  }
}
cat(paste("Success",outputID,indelWarning,cloneWarning,sep="|"))
