# A makefile for setting up the coastwide sim-est procedure 
# (currently only for a mac - will improve later)

# First, let's make sure we know where all the bits are
ADMB="/usr/local/bin/admb"
PROJECT="project"
FORENSICS="badfits"

install: 
	make ADMB
	make dirs

# Run ADMB to set up executables
ADMB:
	$(ADMB) msProdCV.tpl

# Create directories to hold project files and forensics
dirs:
	mkdir $(PROJECT)
	mkdir $(FORENSICS)

# Dry run of cleaning the directory down to the git index. Will show files to be
# deleted
clean:
	git clean -fdn

# Real version of directory cleaner
cleanReal:
	git clean -fd
