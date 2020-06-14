# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile 66 2002-07-14 02:00:42Z gautier $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := 
LIBDDLORDERED := 

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := 

#++ Extra binary rules
$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']

BINS := PacPidNtupleApp
BINCCFILES := $(BINS:=.cc) $(EXTRABINS:=.cc)

#++ Binary rules		 [in production, build it with 'bin']
$(addprefix $(bindir),$(BINS)): $(bindir)% : %.o


#++ shell script products.. 			[build it with 'bin']
#BINSCRIPTS := testscript

#################################################################
#++ regression test scripts			[build it with 'test']

#$(testdir)mytest.T : mytest.tcl mytesttemp2

#########################################
#++ make sure -fno-automatic flag is used

###ifneq (,$(findstring Linux2,$(BFARCH)))
noautomatic_pkg += $(PACKAGE)
##endif ## Linux2

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk

