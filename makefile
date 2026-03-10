.PHONY : all RHSTEST GOMTEST SWSHTEST SFMTEST

# compilation flags and options

CC := g++
#uncomment the 2 lines below if you are using a mac
#CXXFLAGS := -g -O3 -std=c++17 -Wall -Xclang -fopenmp
#LDFLAGS := -Xclang -fopenmp -lomp

#uncomment 2 two linew below if you are using linux
CXXFLAGS := -g -O3 -std=c++17 -Wall  -Wno-reorder -fopenmp 
LDFLAGS := -fopenmp

# Project directory structure
SRC := ./src
MAIN_SRC := ./main_src
INC := ./include
OBJ := ./obj
EXE := ./exe

# Include directories (first row cluster sara, second row local sara)
#INCLUDES := -I./$(INC) -I/home/sgliorio/boost_1_85_0
INCLUDES := -I./$(INC) -I/usr/local/boost_1_85_0

# Boost libraries directory (first row cluster sara, second row local sara)
BOOST_LIB_DIR := /home/sgliorio/boost_1_85_0/libs
#BOOST_LIB_DIR := /usr/local/boost_1_85_0/libs

# Boost libraries
BOOST_LIBS := -L$(BOOST_LIB_DIR) 

# Combined Libraries
LIBS := $(BOOST_LIBS) 

# List of classes and namespaces
CLASSES := RadialHomogeneousSolution GeodesicOrbitalMotion SpinWeightedSpheroidalHarmonics special_functions config ScalarFluxMode

# Define objects from classes
OBJECTS := $(CLASSES:%=$(OBJ)/%.o)

# Define the executables
RHSTEST : $(EXE)/RHS_test
GOMTEST : $(EXE)/GOM_test
SWSHTEST : $(EXE)/SWSH_test
SFMTEST : $(EXE)/SFM_test


#define rule to compile objects from classes
$(OBJ)/%.o: $(SRC)/%.cpp $(INC)/%.h
	@if [ ! -d $(OBJ) ]; then mkdir -pv $(OBJ); fi
	g++ $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

#define rules to compile objects from main sources
$(OBJ)/RHS_test.o : $(MAIN_SRC)/RHS_test.cpp
	@if [ ! -d $(OBJ) ]; then mkdir -pv $(OBJ); fi
	g++ $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/GOM_test.o : $(MAIN_SRC)/GOM_test.cpp
	@if [ ! -d $(OBJ) ]; then mkdir -pv $(OBJ); fi
	g++ $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/SWSH_test.o : $(MAIN_SRC)/SWSH_test.cpp
	@if [ ! -d $(OBJ) ]; then mkdir -pv $(OBJ); fi
	g++ $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/SFM_test.o : $(MAIN_SRC)/SFM_test.cpp
	@if [ ! -d $(OBJ) ]; then mkdir -pv $(OBJ); fi
	g++ $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

# defining the variable all
all : all_sources

test_sources : RHSTEST GOMTEST SWSHTEST SFMTEST
 
all_sources : test_sources

#dependencies and compilation rules for executables

$(EXE)/RHS_test : $(OBJECTS) $(OBJ)/RHS_test.o
	@if [ ! -d $(EXE) ]; then mkdir -pv $(EXE); fi
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@

$(EXE)/GOM_test : $(OBJECTS) $(OBJ)/GOM_test.o
	@if [ ! -d $(EXE) ]; then mkdir -pv $(EXE); fi
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@

$(EXE)/SWSH_test : $(OBJECTS) $(OBJ)/SWSH_test.o
	@if [ ! -d $(EXE) ]; then mkdir -pv $(EXE); fi
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@

$(EXE)/SFM_test : $(OBJECTS) $(OBJ)/SFM_test.o
	@if [ ! -d $(EXE) ]; then mkdir -pv $(EXE); fi
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@

# Clean targets to remove generated files and directories
clean:
	rm -rfv $(OBJ) $(EXE)

# Help target to display available targets and their descriptions
help:
	@echo "  make RHSTEST               - Build RHSTEST"
	@echo "  make GOMTEST                -Build GOMTEST"
	@echo "  make SWSHTEST              - Build SWSHTEST"
	@echo "  make SFMTEST				- Build SFMTEST"
	@echo "  make clean                 - Remove object files and executables"
	@echo "  make all                   - Build all"