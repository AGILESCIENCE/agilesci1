#############################################################################
# Use make variable_name=' options ' to override the variables or make -e to
# override the file variables with the environment variables
#	make CXXFLAGS='-g'
#	make prefix='/usr'
#	make CXX=g++
# External environment variable: CFISIO, ROOTSYS, CTARTA, ICEDIR
# Instructions:
# - modify the section 1)
# - in section 10), modify the following action:
#	* all: and or remove exe and lib prerequisite
#	* lib: and or remove staticlib and dynamiclib prerequisite
#	* clean: add or remove the files and directories that should be cleaned
#	* install: add or remove the files and directories that should be installed
#	* uninstall: add or remove the files and directories that should be uninstalled
#############################################################################

PROJECT= agilesci1
SHELL = /bin/sh

####### 1) Project names and system

SYSTEM= $(shell gcc -dumpmachine)
#ice, ctarta, mpi, cfitsio
LINKERENV= cfitsio, pil, wcs, root, agile

# Applications
AG_ADD_DIFF = AG_add_diff
AG_CHECK_MAP_VALUE = AG_checkMapValue
AG_CIRCLE = AG_circle
AG_MAP2CSV = AG_map2csv
AG_SUMMAPGEN = AG_summapgen
AG_DIFF_CONV = AG_diff_conv
AG_EXPMAPGEN = AG_expmapgen
AG_CTSMAPGEN = AG_ctsmapgen
AG_GASMAPGEN = AG_gasmapgen
AG_INTMAPGEN = AG_intmapgen
AG_DIFMAPGEN = AG_difmapgen
TESTEDP = testedp
AG_AP = AG_ap
AG_LM = AG_lm5
AG_SELECT = AG_select
AG_ITERATIVEGENSRCLIST = AG_iterativeGenSrcList
AG_THETAMAPGEN = AG_thetamapgen
#AG_MULTI2 = AG_multi2
AG_MULTI5 = AG_multi
AG_MULTI5EXT = AG_multi5ext
AG_MULTISIM5 = AG_multisim
AG_DIFFSIM5 = AG_diffsim
AG_MULTITERATIVE5 = AG_multiterative
AG_PASTEMAP5 = AG_pasteMap
AG_PIXEXTRACT= AG_pixextract
AG_SPOTFINDER = AG_spotfinder
AG_GENAPP = AG_genapp
AG_FITPSFARRAY = AG_fitpsfarray
AG_FITPSFARRAY3 = AG_fitpsfarray3
AG_FITPSFARRAY3_H = AG_fitpsfarray3_H
AG_CONVERTTOSKYMAP5 = AG_converttoSkyMap
AG_PADMAP5 = AG_padMap
AG_CREATEPSD = AG_createpsd
AG_CREATEPSD3 = AG_createpsd3
AG_NORM = AG_norm
AG_INDEXGEN = AG_indexgen
AG_EXPRATIO = AG_expratio
AG_LM6 = AG_lm6
AG_ADDRINGTOAITOFFMAP = AG_addringtoaitoffmap

# Libraries
AGILE_MAP = AgileMap
ALIKE_DATA5 = AlikeData5
ALIKE_PSF = AlikePsf
CALIBUTILS = CalibUtils
CONNECTED_REGION = ConnectedRegion
ELLIPSE = Ellipse
FITSUTILS = FitsUtils
GENMAPPARAMS = GenmapParams
INTERVALS = Intervals
LABELING = Labeling
MATH_UTILS = MathUtils
PIL_PARAMS = PilParams
PLOT_CTS2D = PlotCts2D3
POLYGON = Polygon
ROI_MULTI5 = RoiMulti5
ROTATION = Rotation
SKYMAP = SkyMap

VER_FILE_NAME = version.h
#the name of the directory where the conf file are copied (into $(datadir))
CONF_DEST_DIR =
#the name of the icon for the installation
ICON_NAME=

####### 2) Directories for the installation

# Prefix for each installed program. Only ABSOLUTE PATH
prefix=$(AGILE)
exec_prefix=$(prefix)
# The directory to install the binary files in.
bindir=$(exec_prefix)/bin
# The directory to install the local configuration file.
datadir=$(exec_prefix)/share
# The directory to install the libraries in.
libdir=$(exec_prefix)/lib
# The directory to install the info files in.
infodir=$(exec_prefix)/info
# The directory to install the include files in.
includedir=$(exec_prefix)/include
# The directory to install the icon
icondir=$(HOME)/.local/share/applications/

####### 3) Directories for the compiler

OBJECTS_DIR = obj
SOURCE_DIR = code
INCLUDE_DIR = code
DOC_DIR = ref
DOXY_SOURCE_DIR = code_filtered
EXE_DESTDIR  = bin
LIB_DESTDIR = lib
CONF_DIR=conf
ICON_DIR = ui

####### 4) Compiler, tools and options

ifneq (, $(findstring mpi, $(LINKERENV)))
CXX = mpic++
else
CXX = g++
endif

CC = gcc

CXXFLAGS = -g -O2 -pipe -I $(INCLUDE_DIR)

ifneq (, $(findstring agile, $(LINKERENV)))
    ifeq (, $(findstring -I $(AGILE)/include, $(CXXFLAGS)))
        CXXFLAGS += -I $(AGILE)/include
    endif
    LIBS += -L$(AGILE)/lib -lagilesci
endif
ifneq (, $(findstring wcs, $(LINKERENV)))
    ifeq (,$(findstring -I $(AGILE)/include, $(CXXFLAGS)))
        CXXFLAGS += -I $(AGILE)/include
    endif
    LIBS += -L$(AGILE)/lib -lagilewcs
endif
ifneq (, $(findstring pil, $(LINKERENV)))
    ifeq (,$(findstring -I $(AGILE)/include, $(CXXFLAGS)))
        CXXFLAGS += -I $(AGILE)/include
    endif
    LIBS += -L$(AGILE)/lib -lagilepil
endif
ifneq (, $(findstring root, $(LINKERENV)))
    CXXFLAGS += -W -fPIC -D_REENTRANT $(shell root-config --cflags)
    LIBS += $(shell root-config --glibs) -lMinuit
endif
ifneq (, $(findstring cfitsio, $(LINKERENV)))
    CXXFLAGS += -I$(CFITSIO)/include
    LIBS += -L$(CFITSIO)/lib -lcfitsio
endif

LINK     = $(CXX)
#for link
LFLAGS = -shared -Wl,-soname,$(TARGET1) -Wl,-rpath,$(DESTDIR)
AR       = ar cqs
TAR      = tar -cf
GZIP     = gzip -9f
COPY     = cp -f -r
COPY_FILE= $(COPY) -p
COPY_DIR = $(COPY) -pR
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rm -rf
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p

####### 5) VPATH

VPATH=$(SOURCE_DIR):$(INCLUDE_DIR):
vpath %.o $(OBJECTS_DIR)

####### 6) Files of the project

INCLUDE=$(foreach dir,$(INCLUDE_DIR), $(wildcard $(dir)/*.h))
SOURCE=$(foreach dir,$(SOURCE_DIR), $(wildcard $(dir)/*.cpp))
SOURCE+=$(foreach dir,$(SOURCE_DIR), $(wildcard $(dir)/*.c))
#Objects to build
OBJECTS=$(addsuffix .o, $(basename $(notdir $(SOURCE))))
#only for documentation generation
DOC_INCLUDE= $(addprefix $(DOXY_SOURCE_DIR)/, $(notdir $(INCLUDE)))
DOC_SOURCE= $(addprefix $(DOXY_SOURCE_DIR)/, $(notdir $(SOURCE)))

####### 7) Only for library generation

TARGET  = $(LIB_NAME).so.$(shell cat version)
TARGETA = $(LIB_NAME).a
TARGETD = $(LIB_NAME).so.$(shell cat version)
TARGET0 = $(LIB_NAME).so
TARGET1 = $(LIB_NAME).so.$(shell cut version -f 1 -d '.')
TARGET2 = $(LIB_NAME).so.$(shell cut version -f 1 -d '.').$(shell cut version -f 2 -d '.')

####### 8) Preliminar operations

$(shell  cut $(INCLUDE_DIR)/$(VER_FILE_NAME) -f 3 > version)
#WARNING: use -d ' ' if in the version.h the separator is a space

####### 9) Pattern rules

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $(OBJECTS_DIR)/$@

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $(OBJECTS_DIR)/$@

#only for documentation generation
$(DOXY_SOURCE_DIR)/%.h : %.h
	cp  $<  $@

$(DOXY_SOURCE_DIR)/%.cpp : %.cpp
	cp  $<  $@

####### 10) Build rules

#all: compile the entire program.
all: exe
	#only if conf directory is present:
	#$(SYMLINK) $(CONF_DIR) $(CONF_DEST_DIR)

lib: staticlib

exe: makeobjdir $(OBJECTS)
	test -d $(EXE_DESTDIR) || mkdir -p $(EXE_DESTDIR)

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_MULTI5) $(OBJECTS_DIR)/AG_multi.o $(LIBS)

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_MULTISIM5) $(OBJECTS_DIR)/AG_multisim.o $(LIBS)

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_DIFFSIM5) $(OBJECTS_DIR)/AG_diffsim.o $(LIBS)

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_MULTI5EXT) $(OBJECTS_DIR)/AG_multiext.o  $(LIBS)

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(TESTEDP) $(OBJECTS_DIR)/testedp.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_MULTITERATIVE5) $(OBJECTS_DIR)/AG_multiterative.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_CHECK_MAP_VALUE) $(OBJECTS_DIR)/AG_checkMapValue.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_CIRCLE) $(OBJECTS_DIR)/AG_circle.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_MAP2CSV) $(OBJECTS_DIR)/AG_map2csv.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_SUMMAPGEN) $(OBJECTS_DIR)/AG_summapgen.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_DIFF_CONV) $(OBJECTS_DIR)/AG_diff_conv.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_PASTEMAP5) $(OBJECTS_DIR)/AG_pasteMap.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_GENAPP) $(OBJECTS_DIR)/AG_genapp.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_ADD_DIFF) $(OBJECTS_DIR)/AG_add_diff.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_EXPMAPGEN) $(OBJECTS_DIR)/AG_expmapgen.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_AP) $(OBJECTS_DIR)/AG_ap.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_LM) $(OBJECTS_DIR)/AG_lm5.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_SELECT) $(OBJECTS_DIR)/AG_select.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_CTSMAPGEN) $(OBJECTS_DIR)/AG_ctsmapgen.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_GASMAPGEN) $(OBJECTS_DIR)/AG_gasmapgen.o  $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_INTMAPGEN) $(OBJECTS_DIR)/AG_intmapgen.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_DIFMAPGEN) $(OBJECTS_DIR)/AG_difmapgen.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_ITERATIVEGENSRCLIST) $(OBJECTS_DIR)/AG_iterativeGenSrcList.o $(LIBS)

	$(CXX) -g $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_SPOTFINDER) $(OBJECTS_DIR)/AG_spotfinder.o $(LIBS)

	#$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_THETAMAPGEN) $(OBJECTS_DIR)/AG_ $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_FITPSFARRAY) $(OBJECTS_DIR)/AG_fitpsfarray.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_FITPSFARRAY3) $(OBJECTS_DIR)/AG_fitpsfarray3.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_FITPSFARRAY3_H) $(OBJECTS_DIR)/AG_fitpsfarray3_H.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_CONVERTTOSKYMAP5) $(OBJECTS_DIR)/AG_converttoSkyMap.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_PADMAP5) $(OBJECTS_DIR)/AG_padMap.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_CREATEPSD) $(OBJECTS_DIR)/AG_createpsd.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_CREATEPSD3) $(OBJECTS_DIR)/AG_createpsd3.o $(LIBS)

	$(CXX)  $(ALL_CFLAGS_NO_ROOT) -o $(EXE_DESTDIR)/$(AG_NORM) $(OBJECTS_DIR)/AG_norm.o $(LIBS)

	$(CXX) -o $(EXE_DESTDIR)/$(AG_INDEXGEN) $(OBJECTS_DIR)/AG_indexgen.o -lz

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_EXPRATIO) $(OBJECTS_DIR)/AG_expratio.o $(LIBS)

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_LM6) $(OBJECTS_DIR)/AG_lm6.o $(LIBS)

	$(CXX) $(CPPFLAGS) $(ALL_CFLAGS) -o $(EXE_DESTDIR)/$(AG_ADDRINGTOAITOFFMAP) $(OBJECTS_DIR)/AG_addringtoaitoffmap.o $(LIBS)


staticlib: makelibdir makeobjdir $(OBJECTS)
	test -d $(LIB_DESTDIR) || mkdir -p $(LIB_DESTDIR)
	$(DEL_FILE) $(LIB_DESTDIR)/$(TARGETA)
	$(AR) $(LIB_DESTDIR)/$(TARGETA) $(OBJECTS_DIR)/*.o

dynamiclib: makelibdir makeobjdir $(OBJECTS)
	$(DEL_FILE) $(TARGET) $(TARGET0) $(TARGET1) $(TARGET2)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS_DIR)/*.o $(LIBS)
	$(SYMLINK) $(TARGET) $(TARGET0)
	$(SYMLINK) $(TARGET) $(TARGET1)
	$(SYMLINK) $(TARGET) $(TARGET2)
	test $(LIB_DESTDIR) = . || $(DEL_FILE) $(LIB_DESTDIR)/$(TARGET)
	test $(LIB_DESTDIR) = . || $(DEL_FILE) $(LIB_DESTDIR)/$(TARGET0)
	test $(LIB_DESTDIR) = . || $(DEL_FILE) $(LIB_DESTDIR)/$(TARGET1)
	test $(LIB_DESTDIR) = . || $(DEL_FILE) $(LIB_DESTDIR)/$(TARGET2)
	test $(LIB_DESTDIR) = . || $(MOVE) $(TARGET) $(TARGET0) $(TARGET1) $(TARGET2) $(LIB_DESTDIR)

makeobjdir:
	test -d $(OBJECTS_DIR) || mkdir -p $(OBJECTS_DIR)

makelibdir:
	test -d $(LIB_DESTDIR) || mkdir -p $(LIB_DESTDIR)

#clean: delete all files from the current directory that are normally created by building the program.
clean:
	$(DEL_FILE) $(OBJECTS_DIR)/*.o
	$(DEL_FILE) *~ core *.core
	$(DEL_FILE) $(LIB_DESTDIR)/*.a
	$(DEL_FILE) $(LIB_DESTDIR)/*.so*
	$(DEL_FILE) $(EXE_DESTDIR)/AG_*
	$(DEL_FILE) version
	$(DEL_FILE) prefix
	$(DEL_FILE) $(PROJECT).dvi
	$(DEL_FILE) $(PROJECT).pdf
	$(DEL_FILE) $(PROJECT).ps
	test $(OBJECTS_DIR) = . || $(DEL_DIR) $(OBJECTS_DIR)
	test $(EXE_DESTDIR) = . || $(DEL_DIR) $(EXE_DESTDIR)
	test $(LIB_DESTDIR) = . || $(DEL_DIR) $(LIB_DESTDIR)
	test $(DOXY_SOURCE_DIR) = . || $(DEL_DIR) $(DOXY_SOURCE_DIR)
	test $(DOC_DIR) = . || $(DEL_DIR) $(DOC_DIR)

#Delete all files from the current directory that are created by configuring or building the program.
distclean: clean

#install: compile the program and copy the executables, libraries,
#and so on to the file names where they should reside for actual use.
install: all
	$(shell echo $(prefix) > prefix)
	#test -d $(datadir)/$(CONF_DEST_DIR) || mkdir -p $(datadir)/$(CONF_DEST_DIR)
	#test -d $(infodir) || mkdir -p $(infodir)


	# For exe installation
	test -d $(bindir) || mkdir -p $(bindir)
	$(COPY_FILE) $(EXE_DESTDIR)/* $(bindir)
	#copy icon
	#test -d $(icondir) || mkdir -p $(icondir)
	#$(COPY_FILE) $(ICON_DIR)/$(ICON_NAME) $(icondir)

	# For conf files installation
	test -d $(datadir) || mkdir -p $(datadir)
	$(COPY_FILE) $(CONF_DIR)/* $(datadir)/$(CONF_DEST_DIR)

#uninstall: delete all the installed files--the copies that the `install' target creates.
uninstall:
	#For library uninstall
	$(DEL_FILE) $(libdir)/$(TARGETA)
	$(DEL_FILE) $(libdir)/$(TARGETD)
	$(DEL_FILE) $(libdir)/$(TARGET0)
	$(DEL_FILE) $(libdir)/$(TARGET1)
	$(DEL_FILE) $(libdir)/$(TARGET2)
	$(DEL_FILE) $(addprefix $(includedir)/, $(notdir $(INCLUDE)))

	# For exe uninstall
	$(DEL_FILE) $(bindir)/$(EXE_NAME)
	#$(DEL_FILE) $(icondir)/$(ICON_NAME)

#dist: create a distribution tar file for this program
dist: all

# dvi, pdf, ps, for documentation generation
dvi: info
	cd $(DOC_DIR)/latex && $(MAKE)
	$(SYMLINK) $(DOC_DIR)/latex/refman.dvi $(PROJECT).dvi

pdf: info
	cd $(DOC_DIR)/latex && $(MAKE) pdf
	$(SYMLINK) $(DOC_DIR)/latex/refman.pdf $(PROJECT).pdf

ps: info
	cd $(DOC_DIR)/latex && $(MAKE) ps
	$(SYMLINK) $(DOC_DIR)/latex/refman.ps $(PROJECT).ps

#info: generate any Info files needed.
info: makedoxdir $(DOC_INCLUDE) $(DOC_SOURCE)
	test -d $(DOC_DIR) || mkdir -p $(DOC_DIR)
	doxygen Doxyfile

makedoxdir:
	test -d $(DOXY_SOURCE_DIR) || mkdir -p $(DOXY_SOURCE_DIR)
