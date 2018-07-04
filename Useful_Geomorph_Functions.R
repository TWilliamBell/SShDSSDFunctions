### Will's Useful Functions for Geomorph ###
source('./R_Default_Functions/Default_Package_of_Functions.R')
library(geomorph)

cumFact <- function(data,obs,excludeSex = FALSE) {
  # Make a factor that is a combination of the factors it already has to allow each treatment to be
  # its own category.
  # Exclude sex argument requires that the category for sex be named "Sex".
  factorlogical <- unname(unlist(lapply(data,is.factor))) # Tells me whether a given row is a factor
  if (isTRUE(excludeSex)) {
    a <- which(names(data)=="Sex")
    factorlogical[a] <- FALSE
  }
  factors <- array(dim=c(obs,length(data))) # Make an empty array we shall fill
  for (i in 1:length(factorlogical)) {
    if (isTRUE(factorlogical[i])) {
      factors[,i] <- as.character(data[[i]])
      # Fill spots in the array with the columns in the data that are factors
    }
  }
  cumFactor <- rep(NA_character_,obs) # Make empty vector for our end product to fill
  for (i in 1:obs) {
    for (j in 1:length(factorlogical)) {
      if (!is.na(factors[i,j])) {
        cumFactor[i] <- paste(cumFactor[i],as.character(factors[i,j])) # Add to every empty spot in the name of all associated factors via repeated concactenations
      }
    }
  }
  cumFactor <- as.factor(cumFactor)
}

Cumulative_Factors <- function(data, excludeSex = FALSE) { # Takes in a data frame or list and 
  # returns all extant combinations of factors as individual factors.
  factorlogical <- unname(unlist(lapply(data,is.factor)))
  obsnum <- length(data[[which(factorlogical)[1]]])
  A <- cumFact(data,obsnum,excludeSex)
}

# Generic SShD calculation

SShDFunc <- function(Coords, SexVec, Zeroed = F) {

    FData <- Coords[as.character(SexVec)=='f', , drop = FALSE] ## Matrices we will use for overall analysis of sexual dimorphism
    MData <- Coords[as.character(SexVec)=='m', , drop = FALSE]

    ConsensusF <- apply(FData,2,mean) ## Average male and female coordinates
    ConsensusM <- apply(MData,2,mean)

    DiffFM <- ConsensusF-ConsensusM

    SShD <- euclidean(DiffFM)
    
    if (isTRUE(Zeroed)) {
      Zero <- replicate(1000, zeroSShD(Coords, SexVec))
      return(data.frame(SShD = SShD, ZeroTerm = mean(Zero)))
    }
    else {return(SShD)}
    }

# The SSD for the entire data set.

SSDFunc <- function(SizeVec, SexVec) {
  SexVec <- as.character(SexVec)
  MSize <- mean(SizeVec[SexVec == 'm'])
  FSize <- mean(SizeVec[SexVec == 'f'])
  if (FSize >= MSize) {
    SSD <- (FSize/MSize)-1
  }
  else {
    SSD <- -((MSize/FSize)-1)
  }
  SSD
}

# Bootstrapped SShD calculator with optional SSD calculator

bootSShD <- function(Coords, SexVec, SizeVec = NULL) { ## Bootstrapping function.
  n <- nrow(Coords)
  Male <- Coords[as.character(SexVec) == 'm', , drop = F] ## Requires Sex to be given in terms of 'm' and 'f'
  nM <- nrow(Male)
  Female <- Coords[as.character(SexVec) == 'f', , drop = F]
  if (!is.null(SizeVec)) {
    FSize <- SizeVec[as.character(SexVec) == 'f']
    MSize <- SizeVec[as.character(SexVec) == 'm']
  }
  if (length(Male) == 0 & length(Female) == 0) {
    return(NA)
    warning("No data found for either sex, unused level?\n")
  }
  else if (length(Male) == 0 | length(Female) == 0) {
    return(NA)
    warning("No data found for one sex, missing data?\n")
  }
  else {
    nF <- nrow(Female)
    BootF <- sample(1:nF, size = nF, replace = T)
    FBoot <- Female[BootF, , drop = F]
    BootM <- sample(1:nM, size = nM, replace = T)
    MBoot <- Male[BootM, , drop = F]
    Sex <- c(rep('f', nF), rep('m', nM))
    if (!is.null(SizeVec)) {
      FSizeBoot <- FSize[BootF]
      MSizeBoot <- MSize[BootM]
      Size <- c(FSizeBoot,MSizeBoot)
      SSDdat <- SSDFunc(Size,Sex)
    }
    BootFrame <- data.frame(Sex)
    Boot <- matrix(nrow = n, ncol = ncol(Coords))
    Boot[1:nF, ] <- FBoot
    Boot[(nF+1):n, ] <- MBoot
    BootFrame$coords <- Boot
    LM <- lm(coords ~ Sex + 1, data = BootFrame)
    Vector <- LM$coef[2,]
    SShD <- euclidean(Vector)
    if (!is.null(SizeVec)) {
      return(list(SShD = SShD, SSD = SSDdat))
    }
    SShD
  }
}

## SShD error estimation

SShDwSEFunc <- function(Coords, SexVec, SizeVec = NULL, rep = 1000) {
    Boot <- replicate(rep,bootSShD(Coords, SexVec, SizeVec))
    if (!is.null(SizeVec)) {
      SShD <- mean(unlist(Boot[1, ]))
      StandardErrorforSShD <- sd(unlist(Boot[1, ]))
      SSD <- mean(unlist(Boot[2, ]))
      StandardErrorforSSD <- sd(unlist(Boot[2, ]))
      Results <- Boot
    }
    else {
      SShD <- mean(Boot)
      StandardErrorforSShD <- sd(Boot)
      Results <- Boot
    }
    SShD_Results <- list(SShD = SShD, StandardErrorforSShD = StandardErrorforSShD, Results = Results)
    if (!is.null(SizeVec)) {
      SShD_Results$SSD <- SSD
      SShD_Results$StandardErrorforSSD <- StandardErrorforSSD
    }
    SShD_Results
}

## Let's make a stratified SShD calculator.

StratSShDFunc <- function(Coords, NonSexFactorVec, SexVec) {
    NonSexFactorLevelsVec <- levels(NonSexFactorVec)
    StratSShD <- rep(NA_real_,length(NonSexFactorLevelsVec))
    zeroSShD <- rep(NA_real_, length(NonSexFactorLevelsVec))

    for (i in 1:length(NonSexFactorLevelsVec)) { 
      ## Level-based SShD calculations w/o bootstrap
      level <- NonSexFactorLevelsVec[i]
      StratCoord <- Coords[NonSexFactorVec==level, , drop = FALSE]
      if (length(StratCoord) == 0) {
        StratSShD[i] <- NA
        warning(paste("No coordinate data found for th level,", level, ", unused level?\n", sep = ""))
      }
      else {
        TempDataFrame1 <- data.frame(Sex = SexVec[NonSexFactorVec==level])
        TempDataFrame1$coords <- StratCoord[ , , drop = F]
        if (!is.element('m', TempDataFrame1$Sex) | 
            !is.element('f', TempDataFrame1$Sex)) {
          StratSShD[i] <- NA
          warning(paste("No data found for at least one sex for the level ", level, ", missing data?\n", sep = ""))
        }
        else {
          LinearModel1 <- lm(coords ~ Sex+1, data = TempDataFrame1)
          StratSShD[i] <- euclidean(LinearModel1$coefficients[2,])
          zeroSShD[i] <- mean(replicate(1000, zeroSShD(TempDataFrame1$coords, TempDataFrame1$Sex)))
        }
      }
    }

    SShDStratifiedResults <- list()  ## Start collecting the information
    SShDStratifiedResults$Factors <- NonSexFactorLevelsVec
    SShDStratifiedResults$SShD <- StratSShD
    SShDStratifiedResults$ZeroSShDBiasTerm <- zeroSShD
    SShDStratifiedResults
}

## Stratified SShD Calculator with Standard Error, optional SSD calculator.

StratSShDwSEFunc <- function(Coords, NonSexFactorVec, SexVec, SizeVec = NULL, print.progress = T, rep = 1000) {
    StratSShDmean <- rep(NA_real_, nlevels(NonSexFactorVec))
    StratSShDse <- rep(NA_real_, nlevels(NonSexFactorVec))
    StratSSDmean <- rep(NA_real_, nlevels(NonSexFactorVec))
    StratSSDse <- rep(NA_real_, nlevels(NonSexFactorVec))
    nvec <- unname(table(NonSexFactorVec))
    Levels <- levels(NonSexFactorVec)
    for (i in 1:length(Levels)) {
      ## Determine SShD for each group and standard error using bootstrap analysis.
      dat <- Coords[NonSexFactorVec == Levels[i], , drop = F]
      SexVector <- SexVec[NonSexFactorVec == Levels[i]]
      a <- replicate(rep,bootSShD(
        Coords = dat, SexVec = SexVector, SizeVec = SizeVec))
      if (!is.null(SizeVec) & !(typeof(a)=="logical")) {
          StratSSDmean[i] <- mean(unlist(a[2, ]))
          StratSSDse[i] <- sd(unlist(a[2, ]))
          StratSShDmean[i] <- mean(unlist(a[1, ]))
          StratSShDse[i] <- sd(unlist(a[1, ]))
      }
      else {
        StratSShDmean[i] <- mean(a)
        StratSShDse[i] <- sd(a)
      }
      if (isTRUE(print.progress)) {
        if (i==1) {
          cat("Finished the ", i,"st level of ", length(Levels), " levels total for SShD bootstrap.\n", sep="")
        }
        else if (i==2) {
          cat("Finished the ", i,"nd level of ", length(Levels), " levels total for SShD bootstrap.\n", sep="")
        }
        else if (i==3) {
          cat("Finished the ", i,"rd level of ", length(Levels), " levels total for SShD bootstrap.\n", sep="")
        }
        else {cat("Finished the ", i,"th level of ", length(Levels), " levels total for SShD bootstrap.\n", sep="")}
      }
    }
  SShDwSEStratifiedResults <- list()
  if (!is.null(SexVec)) {
    SShDwSEStratifiedResults$BootstrappedSSD <- StratSSDmean
    SShDwSEStratifiedResults$SSD_Standard_Error <- StratSSDse
  }
  SShDwSEStratifiedResults$BootstrappedSShD <- StratSShDmean
  SShDwSEStratifiedResults$SShD_Standard_Error <- StratSShDse
  SShDwSEStratifiedResults
}

## Let's make a function that determines SSD

StratSSDFunc <- function(SizeVec, NonSexFactorVec, SexVec) {
    Levels <- levels(NonSexFactorVec) 
    StratSSD <- rep(NA_real_,length(Levels)) 
    for (i in 1:length(Levels)) {
      level <- Levels[i]
      StratFSize <- SizeVec[as.character(SexVec)=='f' & NonSexFactorVec==level] 
      StratMSize <- SizeVec[as.character(SexVec)=='m' & NonSexFactorVec==level]
      if (length(StratFSize) == 0 | length(StratMSize) == 0) {
        StratSSD[i] <- NA
        warning(paste("No data found for at least one sex for ", level, ", missing data?\n", sep = ""))
      }
      else {
        StratSSD[i] <- SSDFunc(SizeVec[NonSexFactorVec == level], SexVec[NonSexFactorVec == level])
      }
    }
    StratSSD
}

zeroSShD <- function(Coords, SexVec) {
  Sex <- sample(SexVec, length(SexVec))
  NullSShD <- SShDFunc(Coords, Sex)
  NullSShD
}

read.xl.allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

clear()