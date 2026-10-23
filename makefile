.PHONY: all all_sources test_sources RHSTEST GOMTEST SWSHTEST SFMTEST MSTREGRESSION COMPACT_TEST clean help
all: all_sources

CC := g++
CXXFLAGS := -g -O3 -std=c++17 -Wall -Wno-reorder
LDFLAGS :=

SRC := ./src
MAIN_SRC := ./main_src
INC := ./include
OBJ := ./obj
EXE := ./exe

# Set these paths to your Boost installation.
INCLUDES := -I$(INC) -I/YOUR/PATH/TO/boost_VERSION_NUMBER
BOOST_LIB_DIR := /YOUR/PATH/TO/boost_VERSION_NUMBER/libs
LIBS := -L$(BOOST_LIB_DIR)

CLASSES := RadialHomogeneousSolution GeodesicOrbitalMotion SpinWeightedSpheroidalHarmonics special_functions config ScalarFluxMode CompactifiedCharacteristicRadialSolution ScalarRadialBackend
OBJECTS := $(CLASSES:%=$(OBJ)/%.o)

RHSTEST: $(EXE)/RHS_test
GOMTEST: $(EXE)/GOM_test
SWSHTEST: $(EXE)/SWSH_test
SFMTEST: $(EXE)/SFM_test
MSTREGRESSION: $(EXE)/MST_regression_test
COMPACT_TEST: $(EXE)/compactified_test
all_sources: test_sources
test_sources: RHSTEST GOMTEST SWSHTEST SFMTEST

$(OBJECTS): $(wildcard $(INC)/*.h)
$(OBJ)/RHS_test.o $(OBJ)/GOM_test.o $(OBJ)/SWSH_test.o $(OBJ)/SFM_test.o $(OBJ)/MST_regression_test.o $(OBJ)/compactified_test.o: $(wildcard $(INC)/*.h)

$(OBJ)/%.o: $(SRC)/%.cpp $(INC)/%.h
	@mkdir -p $(OBJ)
	$(CC) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/%_test.o: $(MAIN_SRC)/%_test.cpp
	@mkdir -p $(OBJ)
	$(CC) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/compactified_test.o: tests/compactified_test.cpp
	@mkdir -p $(OBJ)
	$(CC) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(EXE)/%_test: $(OBJECTS) $(OBJ)/%_test.o
	@mkdir -p $(EXE)
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@

clean:
	rm -rf $(OBJ) $(EXE)

help:
	@echo "make RHSTEST GOMTEST SWSHTEST SFMTEST MSTREGRESSION COMPACT_TEST all clean"
