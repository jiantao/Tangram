SOURCES:=$(shell ls *.cpp)
OUT:=TGM_Reference.o TGM_Error.o TGM_GetOpt.o md5.o 
OBJS:=$(addprefix $(OBJ_DIR)/,$(SOURCES:.cpp=.o))
OUT_OBJS:=$(addprefix $(OBJ_DIR)/,$(OUT))

LOCAL_INCLUDES:= -I../TangramDetect

PROGRAM:=$(BIN_DIR)/tangram_index

$(PROGRAM): $(OBJS) $(OUT_OBJS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LOCAL_INCLUDES) -lz

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(LOCAL_INCLUDES)


.PHONY: $(OUT_OBJS)
