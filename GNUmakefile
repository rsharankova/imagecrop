
# compiler
CXX = g++

CXXFLAGS = -g -fPIC

# include options for this package
INCFLAGS  = -I.
INCFLAGS += $(shell root-config --cflags)
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larlite-config --includes)/../UserDev
INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool
INCFLAGS += $(shell larlite-config --includes)/../UserDev/SelectionTool
INCFLAGS += -I$(LARLITE_USERDEVDIR)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell larlitecv-config --includes)
##INCFLAGS += -I$(LAROPENCV_BASEDIR)
INCFLAGS += -I$(LARCV_INCDIR)/LArOpenCVHandle/

LDFLAGS += $(shell root-config --ldflags --libs)
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_FhiclLite
LDFLAGS += $(shell larlite-config --libs) -lSelectionTool_OpT0FinderAna -lSelectionTool_OpT0FinderApp \
	-lSelectionTool_OpT0PhotonLibrary -lSelectionTool_OpT0FinderAlgorithms -lSelectionTool_OpT0FinderBase 
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)

#LDFLAGS += -lLArOpenCV_Core
#LDFLAGS += -lLArOpenCV_ImageClusterBase
#LDFLAGS += -lLArOpenCV_ImageClusterAlgoFunction
#LDFLAGS += -lLArOpenCV_ImageClusterAlgoData
#LDFLAGS += -lLArOpenCV_ImageClusterAlgoClass

CXXFLAGS += $(INCFLAGS)

# note: llcvprocessor headers and libraries are in larlitecv/build/include and lib
#LDFLAGS += -l

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

PROGRAMS = larflow_imgcrop larflow_imgcrop_data
PROGSRC = $(addsuffix .cxx, $(PROGRAMS) )
SRCS = $(filter-out $(PROGSRC),$(wildcard *.cxx))
OBJS = $(addprefix obj/,$(patsubst %.cxx, %.o, $(SRCS)))
BINS  = $(addprefix bin/,$(PROGRAMS))

all: $(OBJS) $(BINS)


clean:
	rm $(BINS)

bin/%: %.cxx $(OBJS)
	$(CXX) $(CXXFLAGS) $*.cxx -o $@ $(OBJS) $(LDFLAGS)

obj/%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $*.cxx -o $@
