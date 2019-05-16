#  Copyright (c) 2018 - RWTH Aachen University
#
#  File author(s): Enio Gjerga
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  email: enio.gjerga@gmai.com
#
##############################################################################
# 13:22 23/03/2018
# This script is used to produce the PHONEMeS Outputs by calling the functions
# which write the ILP formulation

start.time <- Sys.time()
#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)
library(dplyr)
library(XML)

##
# Loading the prepared PHONEMeS inputs and background network

#load("../Input-Data/dataGMM_UACC257.RData") # for UACC cell line
load("../Input-Data/dataGMM_A2058.RData") # for A2058 cell line
load("../Background-Network/interactionsALL.RData")

##
# Loading all functions needed for the analysis

source("Public/buildDataMatrix.R")
source("Public/ilpFunctions.R")
source("Public/buildDataObject.R")
source("Public/build_Nw.R")
source("Public/build_PKN.R")
source("Public/runPHONEMeS_dt.R")
source("Public/assignAttributes.R")

##
# Preparing background network as a PHONEMeS input

bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)

##
# Choose the conditions for each time-point

#conditions <- list(c("tp_2min"), c("tp_10min"), c("tp_30min"), c("tp_60min"), c("tp_90min")) # for uacc
conditions <- list(c("tp_2min"), c("tp_10min"), c("tp_20min"), c("tp_30min"), c("tp_60min"), c("tp_90min")) # for a2058

#names(conditions) <- c("tp_2min", "tp_10min", "tp_30min", "tp_60min", "tp_90min") # for uacc
names(conditions) <- c("tp_2min", "tp_10min", "tp_20min", "tp_30min", "tp_60min", "tp_90min") # for a2058

##
# Choose the targets for each time-point (EDNRB - the same)

#targets.P <- list(tp_2min=c("EDNRB_HUMAN"), tp_10min=c("EDNRB_HUMAN"), tp_30min=c("EDNRB_HUMAN"), tp_60min=c("EDNRB_HUMAN"), tp_90min=c("EDNRB_HUMAN")) # for uacc
targets.P <- list(tp_2min=c("EDNRB_HUMAN"), tp_10min=c("EDNRB_HUMAN"), tp_20min=c("EDNRB_HUMAN"), tp_30min=c("EDNRB_HUMAN"), tp_60min=c("EDNRB_HUMAN"), tp_90min=c("EDNRB_HUMAN")) # for a2058

##
# Assign to each time-point a set of experiments accordingly

#experiments <- list(tp_2min=c(1), tp_10min=c(2), tp_30min=c(3), tp_60min=c(4), tp_90min=c(5)) # for uacc
experiments <- list(tp_2min=c(1), tp_10min=c(2), tp_20min=c(3), tp_30min=c(4), tp_60min=c(5), tp_90min=c(6)) # for a2058


#set.seed(383789) # for UACC cell line
set.seed(102384) # for A2058 cell line

resList = runPHONEMeS_dt(targets.P = targets.P, conditions = conditions, dataGMM = dataGMM, experiments = experiments, bg = bg, nIter = 100)

#save(resList, file = "resList_UACC257.RData")
save(resList, file = "resList_A2058.RData")
