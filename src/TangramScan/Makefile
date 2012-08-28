SOURCES  := TGM_ReadPairScan.c TGM_ReadPairScanGetOpt.c TGM_ReadPairScanMain.c TGM_LibInfo.c TGM_FragLenHist.c TGM_GetOpt.c TGM_Error.c TGM_Utilities.c TGM_BamInStream.c \
	    TGM_BamHeader.c TGM_BamPairAux.c TGM_BamMemPool.c

OBJS = $(SOURCES:.c=.o)
PROGRAM:=$(BIN_DIR)/tangram_scan

$(PROGRAM): $(OBJS)
	@echo "  * linking $(PROGRAM)"
	@$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(INCLUDES) -lbam -lz -lm

$(OBJS): $(SOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c $(CFLAGS) $(INCLUDES)

clean:
	@rm -rf *.o

.PHONY: clean
