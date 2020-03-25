SHELL = /bin/sh

#==========================
# Files
#==========================
EXE = test_snap.exe
SRC = *.cpp

#==========================
# Machine specific info
#==========================
# compilers and options
CXX = clang++
CXXFLAGS = -g -O2
DEFINE = -DREFDATA_TWOJ=$(ref_data) # 8 or 14 or 2
CXXFLAGS = -g -O2 $(DEFINE) -std=c++11
CXXFLAGS += -mfma
ifeq ($(OPENMP),y)
	DEFINE += -Dopenmp_version
	DEFINE += -ffp-contract=fast -fstrict-aliasing -Wno-openmp-target -Wall -Wno-unused-variable -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda
endif

#==========================
# Compiler commands
#==========================
CXXOBJ       = $(CXX) $(CXXFLAGS) -c 
CXXLD         = $(CXX) $(CXXFLAGS)

#==========================
# Make the executable
#==========================

$(EXE): $(SRC) $(INC)
	echo $(SRC)
	$(CXXLD) $(SRC) -o $(EXE)

#==========================
#remove all objs
#==========================
clean:
	/bin/rm -f *.o $(EXE)
