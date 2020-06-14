
# root output filename
set HistFileName pacPidNtupleKaon.root

# number of events to be simulated
#set NEVENTS 10
set NEVENTS 1000

######################
### PacMC Sequence ###
######################

#
#  create the PacMC main sequence
#
sourceFoundFile PacMC/PacMC.tcl

#############
### Decay ###
#############

## Customize generator
## Note on DECAY.DEC: DECAY.DEC in PacMC was modified to use the Ali-Greub model
## instead of Kagan-Neubert for b->Xs gamma. This saves ~2min initialization time.
## Edit the lines containing "XSGAMMA" or use PARENT/EvtGen/DECAY.DEC to switch to Kagan-Neubert.
disableGenerators 0
module enable GfiEvtGen
talkto GfiEvtGen {
    GENERATE set Upsilon(4S)
    DECAY    set RELEASE/PacMC/DECAY.DEC
    UDECAY   set PARENT/ProdDecayFiles/B+B-_DstarD0_D0p_Kpi.dec
}

#######################
### PID ntuple part ###
#######################

# From BetaPidCalibNtuple/BtaPCNtupleMiniKaon.tcl

## choose the flavor of ntuple to write (hbook or root) and the file name
##
set BetaMiniTuple "root"

# list of samples
# This list is read in BetaPidCalibNtuple/BtaPCNtupleMiniSetup.tcl
set SampleList [ list PacPidDstarSample ]
set TagBitList [ list PacPidDstarSample ]

# The following statement enables the Tag-filtermode
set TagBitFiltering 0

# setup PacPidPCNtupleWriter
module talk PacPidPCNtupleWriter
  doMuha        set false
  doDstar       set true
  doLambda      set false
  doTrkLambda   set false
  doConversions set false
exit

# setup everything else and run
sourceFoundFile PacPidCalib/PacPidNtupleSetup.tcl

#action enable NameAction

#sourceFoundFile BetaPidCalibNtuple/BtaPCNtupleMiniRun.tcl

echo
echo "*******************"
echo "*** action list ***"
echo "*******************"
echo

action list

echo
echo "*****************"
echo "*** path list ***"
echo "*****************"
echo

path list

echo
echo "****************"
echo "*** seq list ***"
echo "****************"
echo

seq list

echo
echo "****************"
echo "*** mod list ***"
echo "****************"
echo

mod list

echo

if [info exists NEVENTS] {
    ev begin -nev $NEVENTS
    ErrMsg trace "completed OK"
    exit
}

