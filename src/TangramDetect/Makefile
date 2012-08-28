SOURCES:=$(shell ls *.cpp)
CSOURCES:=$(shell ls *.c)

OBJS:=$(addprefix $(OBJ_DIR)/,$(SOURCES:.cpp=.o))
COBJS:=$(addprefix $(OBJ_DIR)/,$(CSOURCES:.c=.o))

PROGRAM:=$(BIN_DIR)/tangram_detect

$(PROGRAM): $(OBJS) $(COBJS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(CXXFLAGS) $(PTHREAD) -o $@ $^ $(INCLUDES) -lbamtools -lz

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(INCLUDES)

$(COBJS): $(CSOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c $(CFLAGS) $(INCLUDES)
