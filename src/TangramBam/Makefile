SOURCES:=$(shell ls *.cpp)

OBJS:=$(addprefix $(OBJ_DIR)/,$(SOURCES:.cpp=.o))

PROGRAM:=$(BIN_DIR)/tangram_bam

REQUIRED_OBJS = $(OBJ_DIR)/Fasta.o \
		$(OBJ_DIR)/split.o \
		$(OBJ_DIR)/ssw_cpp.o \
		$(OBJ_DIR)/ssw.o

$(PROGRAM): $(OBJS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(CXXFLAGS) $(PTHREAD) -o $@ $^ $(INCLUDES) $(REQUIRED_OBJS) -lbamtools -lz

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(INCLUDES)

