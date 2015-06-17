SYSTEM     = x86_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
# Cplex & Concert Directorys
#------------------------------------------------------------

#CPLEX         = ./../../../../../../../../opt/ibm/ILOG/CPLEX_Studio126
CPLEX         = /opt/ibm/ILOG/CPLEX_Studio126
CPLEXDIR      = $(CPLEX)/cplex
CONCERTDIR    = $(CPLEX)/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m32 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -m32 -lm -lpthread

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
INCLUDE     = ./../include

EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
DATA        = /data

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(INCLUDE)

#------------------------------------------------------------

OUT=./../../bin
OBJS:=$(patsubst %.cpp,$(OUT)/%.o,$(wildcard *.cpp))
CPP_EX = $(OUT)/PMCLAP

all:
	make all_cpp

all_cpp: $(CPP_EX)

# ------------------------------------------------------------

$(CPP_EX): $(OBJS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) $(OBJS) -o $@ $(CCLNFLAGS)

$(OUT)/%.o: %.cpp $(INCLUDE)/%.h
	$(CCC) $(CCFLAGS) -c $< -o $@

# ------------------------------------------------------------
.PHONY: clean
clean:
	rm $(OUT)/*.o
