BAM_LIB:=$(LIB_DIR)/libbamtools.a
SAM_LIB:=$(LIB_DIR)/libbam.a
SSW:=$(OBJ_DIR)/ssw.o

#SUB_DIRS = stripedSW

libs: $(SSW) $(BAM_LIB) $(SAM_LIB)

$(SSW): stripedSW/ssw.c stripedSW/ssw_cpp.cpp
	@$(MAKE) all --no-print-directory --directory=stripedSW; 

$(BAM_LIB): 
	@echo "  * generating libbamtools.a"
	@cd bamtools/ && mkdir -p build && cd build && cmake .. && $(MAKE)
	@echo "  * copying libbamtools.a to the library directory"
	@cp bamtools/lib/libbamtools.a $(LIB_DIR)
	@echo "  * copying include files to the library directory"
	@cp -r bamtools/include/* $(INCLUDE_DIR)
	@echo " "

$(SAM_LIB):
	@echo "  * generating libbam.a"
	@$(MAKE) lib --no-print-directory -C samtools
	@echo "  * copying libbam.a to the library directory"
	@cp samtools/libbam.a $(LIB_DIR)
	@echo "  * copying include files to the library directory"
	@cp samtools/bgzf.h samtools/bam.h samtools/khash.h samtools/kseq.h $(INCLUDE_DIR)

clean:
	@rm -rf bamtools/bin bamtools/build bamtools/lib/ bamtools/include
	@$(MAKE) clean --no-print-directory -C samtools

.PHONY: all libs clean
