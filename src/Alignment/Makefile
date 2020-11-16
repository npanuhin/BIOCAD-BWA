CXX = g++
CXXFLAGS =
RM = del /Q /F /S
RM_DIR = rd

OBJ_DIR = obj
BIN_DIR = bin

DIRS = $(OBJ_DIR) $(BIN_DIR)
SRC = $(wildcard *.cpp)
OBJ = $(SRC:%.cpp=$(OBJ_DIR)/%.o)
# OBJS = $(addprefix $(OBJ_DIR)/,$(addsuffix .0,$(basename $(pathsubst %,%,$(SOURCES)))))



all: compile

compile: $(DIRS) $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(BIN_DIR)/run

run:
	$(BIN_DIR)/run

clean: $(DIRS)
	$(RM) *.o
	$(RM) *.exe
	$(RM_DIR) bin
	$(RM_DIR) obj

# debug:
# 	echo $(OBJ)



$(DIRS):
	mkdir $@

$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@



K_MERS = k-mers.cpp utils.cpp main.cpp
K_MERS_OBJ = $(K_MERS:%.cpp=$(OBJ_DIR)/%.o)
k-mers: $(DIRS) $(K_MERS_OBJ)
	$(CXX) $(CXXFLAGS) $(K_MERS_OBJ) -o $(BIN_DIR)/k-mers
	$(BIN_DIR)/k-mers
