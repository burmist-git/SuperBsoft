# $Id: Exp $
# --------------------------------------------------------
# setup selectors for BPC - ntuples
# Look whether YOUR selector is
#  a) enabled
#  b) added to 'SelectorList'
# --------------------------------------------------------


# Good Tracks
mod enable GoodTrackVeryLooseSelection
mod enable GoodTrackLooseSelection
mod enable GoodTrackTightSelection
mod enable GoodPhotonLooseSelection
mod enable GoodTrackAccSelection
mod enable GoodTrackAccLooseSelection
mod enable GoodNeutralLooseAccSelection
mod enable GoodPhotonDefaultSelection
mod disable TrkMicroDispatch

# Muons
sequence enable PacPidMuonSequence
module enable TruthBasedMuonSelection

# Electrons
sequence enable PacPidElectronSequence
module enable  NoDeDxFirstElectronSelection

# Kaons
sequence enable  PacPidKaonSequence
module enable NotApionFirstKaonSelection
module enable VeryLooseFirstKaonSelection
module enable LooseFirstKaonSelection
module enable TightFirstKaonSelection
module enable VeryTightFirstKaonSelection
module enable TruthBasedKaonSelection

# pions
sequence enable  PacPidPionSequence
module enable VeryLooseFirstPionSelection
module enable LooseFirstPionSelection
module enable TightFirstPionSelection
module enable VeryTightFirstPionSelection
module enable TruthBasedPionSelection

#protons

exit


set SelectorList [ list ChargedTracks GoodTracksVeryLoose GoodTracksLoose GoodTracksTight \
		       NoDeDxFirstElectronSelection \
		       NotApionFirstKaonSelection VeryLooseFirstKaonSelection LooseFirstKaonSelection TightFirstKaonSelection VeryTightFirstKaonSelection TruthBasedKaonSelection \
		       VeryLooseFirstPionSelection LooseFirstPionSelection TightFirstPionSelection VeryTightFirstPionSelection TruthBasedPionSelection
		        ]

# BEWARE The maximum size of the PidLists arrays is fixed
# BEWARE in BtaCand2Ntuple.cc -- now set to 80 because...
#              dyaeb 20060605 -- we have just exceeded 60!
# BEWARE in BtaCand2Ntuple.cc -- now set to 100 in anticipation of new selector families, 17 nov 06 (kflood)
#   19 Aug 2007: with muBDT and without Kalanand's new selectors, we already have 80! (avtelnov)
#   21 Aug 2007: !Ay, caramba! it's exactly 100!!! (kalanand, avtelnov)

# aparently you can not put comments in the middle of a list...
# next time please test your tcl changes
#                        pMicroLoose pMicroDefault pMicroTight \
# below removed for release 18   dyaeb 20050217
#            muLikeVeryLoose muLikeLoose muLikeTight \
#                        KMicroNotPionGTL
#                piRoyLoose piRoyNotKaon \
#  K, pi & p GLHTight added for release 18
#  K, pi & p GLH expanded to usual 4 levels for release 20 & analysis-31
#  p ELH added for release 20 & analysis-31

if [ info exists ThisisBPCWriterMC ] {
    foreach s $SelectorList {
	mod talk  BtaPCNtupleWriterMC
	PidLists set $s
	exit
    }
} else {
    if [ info exists ThisisBPCWriter_dEdxKp ] {
	foreach s $SelectorList {
	    mod talk  BtaPCNtupleWriter_dEdxKp
	    PidLists set $s
	    exit
	}
    } else {
	foreach s $SelectorList {
	    mod talk  PacPidPCNtupleWriter
	    PidLists set $s
	    exit
	}
    }
}
