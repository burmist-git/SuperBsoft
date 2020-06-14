#-----------------------------------------------------
# PacPidNtupleSetup.tcl
#
# Main setup script for selecting BPC - samples in CM2
#   Moving from $env to FwkCfgVar, but not there yet!
#
# $Id: Exp $
#
# ----------------------------------------------------

# ----------------------------------------------
# Enable Filtering on Tag-DB
# ----------------------------------------------
###unset TagBitFiltering 

# if its available, execute the Patches.tcl file from the current directory
catch { source Patches.tcl }

# ----------------------------------------------
# PART TWO : Init PID - Stuff 
# ----------------------------------------------

# disable some selectors and notify SampleBPCAnalysis
sourceFoundFile PacPidCalib/PacPidNtupleSelectors.tcl

# ----------------------------------------------
# PART THREE : Init BtaPidCalib - Stuff 
# ----------------------------------------------

sequence create PacPidPCNtupleSequence
sequence append PacPidPCNtupleSequence PacPidExampleCalibProcessor
sequence append PacPidPCNtupleSequence PacPidPCNtupleWriter
mod disable PacPidPCNtupleWriterMC

if [ info exists MC ] {

 mod talk PacPidPCNtupleWriter
    writeMcTruth set t
 exit
}

# append MyAnalysis to path
sequence append MyAnalysisSequence PacPidPCNtupleSequence

# -----------------------------------------------------
# PART FOUR : Configure Modules
# -----------------------------------------------------

# add some tag bits
sourceFoundFile PacPidCalib/PacPidPCNtupleTagBits.tcl

# add other tag bits 
if [ info exists TagBitList ] {
  foreach s $TagBitList {
    module talk PacPidPCNtupleWriter 
      TagBits set $s
     exit
  }
}


if [ expr [info exists env(BtaPCNtupleTight)] && {$env(BtaPCNtupleTight) == "yes"} ] {
    module talk PacPidExampleCalibProcessor
      pidKsFillNtp set f
      pidDstarFillNtp set f
    exit
} else {
    module talk PacPidExampleCalibProcessor
      # NA: commenting out
      #sourceFoundFile BetaPidCalib/BtaPidKs.tcl

      sourceFoundFile BetaPidCalib/BtaPidDstarNotSoTight.tcl
      #sourceFoundFile BetaPidCalib/BtaPidDstar.tcl
      # loose proton setup

      # NA: commenting out
      #pidLamPr_mode set loose
      #pidLamPr_flight set 0.5
      #pidLamPr_deltaVtxIPCut set 0.3   
   
      #LambdaTrkOnlySamplePostVtxMassWindow set 0.1
      #LambdaTrkOnlySamplePointBackMax set 0.1
    
      # NA: commenting out
      #pidKsFillNtp set f

      pidDstarFillNtp set f
    exit
}


if [ info exists SampleList ] {
  foreach s $SampleList {
     module talk PacPidExampleCalibProcessor 
      sampleNames set $s
     exit
  }
}

# setup conversion selector 
#sourceFoundFile  BetaPidCalib/BtaPidGammaConv.tcl

#silly BPC - everything is different in opr and in the ntuple
#maker.  So disable the opr sequence to make ntuples.
#sequence disable BtaPidCalibOprSequence
