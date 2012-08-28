SOURCES  := TGM_MergeLibGetOpt.c TGM_MergeLibMain.c TGM_LibInfo.c TGM_FragLenHist.c TGM_GetOpt.c TGM_Error.c TGM_Utilities.c
OBJS = $(SOURCES:.c=.o)

OBJS = $(SOURCES:.c=.o)
PROGRAM:=$(BIN_DIR)/tangram_merge

$(PROGRAM): $(OBJS)
	@echo "  * linking $(PROGRAM)"
	@$(CC) $(CFLAGS) -o $(@) $(^) $(INCLUDES) -lbam -lz -lm

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c $(CFLAGS) $(INCLUDES)

clean:
	@rm -rf *.o

.PHONY: clean
