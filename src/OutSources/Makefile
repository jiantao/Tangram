BAM_LIB:=$(LIB_DIR)/libbamtools.a
SAM_LIB:=$(LIB_DIR)/libbam.a

libs: $(BAM_LIB) $(SAM_LIB)

$(BAM_LIB): 
	@cd bamtools/ && mkdir -p build && cd build && cmake .. && $(MAKE)
	@cp bamtools/lib/libbamtools.a $(LIB_DIR)
	@cp -r bamtools/include/* $(INCLUDE_DIR)

$(SAM_LIB):
	@$(MAKE) lib --no-print-directory -C samtools
	@cp samtools/libbam.a $(LIB_DIR)
	@cp samtools/bgzf.h samtools/bam.h samtools/khash.h samtools/kseq.h $(INCLUDE_DIR)

clean:
	@rm -rf bamtools/bin bamtools/build bamtools/lib/ bamtools/include
	@$(MAKE) clean --no-print-directory -C samtools

.PHONY: libs clean
