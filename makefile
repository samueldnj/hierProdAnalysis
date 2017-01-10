# A makefile for setting up the coastwide sim-est procedure 
# (currently only for a mac - will improve later)

# First, let's make sure we know where all the bits are
RSCRIPT="Rscript"
PROJECT="project"
FORENSICS="badfits"
INIT="Rscript ./init.R"

install: 
	make TMB
	make dirs

# Run ADMB to set up executables
TMB:
	$(RSCRIPT) init.R

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
