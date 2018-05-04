######################################################################
# Project: TMVA - Toolkit for Multivariate Data Analysis             #
# Code   : Example usage for TMVA analyses                           #
###################################################################### 
CXX = g++
CCFLAGS = -O0 -std=c++1y #-std=c++11 #-ggdb -g


LD = g++
LDFLAGS = 

LIBS = $(shell root-config --libdir) 
INCS = $(shell root-config --incdir)  


default : TMVAClassification

TMVAClassification : ./TMVAClassification.C
	@echo -n "Building $@ ... "
	$(CXX) $(CCFLAGS) $<  $(addprefix -I, $(INCS))  $(shell root-config --libs) -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO  -lMLP -lm  -o $@
	@echo "Done"
clean: 
	rm -f *.o TMVAClassification

