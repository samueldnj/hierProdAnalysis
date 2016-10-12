# A makefile for setting up the coastwide sim-est procedure 
# (currently only for a mac - will improve later)

# First, let's make sure we know where all the bits are
ADMB="/usr/local/bin/admb"
PROJECT="project"
FORENSICS="badfits"

install: 
	make ADMB
	make dirs


ADMB:
	$(ADMB) msProdCV.tpl
	$(ADMB) ssProdCV.tpl

dirs:
	mkdir $(PROJECT)
	mkdir $(FORENSICS)

# clean:
