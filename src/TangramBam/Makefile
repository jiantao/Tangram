SOURCES:=$(shell ls *.cpp)
CSOURCES:=$(shell ls *.c)

OBJS:=$(addprefix $(OBJ_DIR)/,$(SOURCES:.cpp=.o))
COBJS:=$(addprefix $(OBJ_DIR)/,$(CSOURCES:.c=.o))

PROGRAM:=$(BIN_DIR)/tangram_bam

REQUIRED_OBJS = $(OBJ_DIR)/Fasta.o \
		$(OBJ_DIR)/split.o \
		$(OBJ_DIR)/ssw_cpp.o \
		$(OBJ_DIR)/ssw.o \
		$(OBJ_DIR)/md5.o 

$(PROGRAM): $(OBJS) $(COBJS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(CXXFLAGS) $(PTHREAD) -o $@ $^ $(INCLUDES) $(REQUIRED_OBJS) -lbamtools -lbam -lz

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) -std=c++0x $(INCLUDES)

$(COBJS): $(CSOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c $(CFLAGS) $(INCLUDES)

