SOURCES  := TGM_MergeLibGetOpt.c TGM_MergeLibMain.c TGM_LibInfo.c TGM_FragLenHist.c TGM_GetOpt.c TGM_Utilities.c
OBJS = $(SOURCES:.c=.o)

REQUIRED_OBJS = $(OBJ_DIR)/TGM_Error.o \
                $(OBJ_DIR)/TGM_BamHeader.o

OBJS = $(SOURCES:.c=.o)
PROGRAM:=$(BIN_DIR)/tangram_merge

$(PROGRAM): $(OBJS)
	@echo "  * linking $(PROGRAM)"
	@$(CC) $(CFLAGS) -o $(@) $(^) $(REQUIRED_OBJS) $(INCLUDES) -lbam -lz -lm

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c $(CFLAGS) $(INCLUDES) -I ../OutSources/util/

clean:
	@rm -rf *.o

.PHONY: clean
