$(shell mkdir -p ../bin)
$(shell mkdir -p ../obj)
$(shell mkdir -p ./lib)
$(shell mkdir -p ./include)

export SRC_DIR:=$(shell pwd)
export OBJ_DIR:=$(dir $(SRC_DIR))obj
export BIN_DIR:=$(dir $(SRC_DIR))bin
export LIB_DIR:=$(SRC_DIR)/lib
export INCLUDE_DIR:=$(SRC_DIR)/include

export CC := gcc
export CXX := g++

### Please choose NORMAL or VERBOSE_DEBUG mode ###
##############################
##### NORMAL mode ############
##############################
export CFLAGS := --std=gnu99 -Wall -g -O2
export CXXFLAGS := -Wall -g -O2
##############################
##### VERBOSE_DEBUG mode #####
##############################
#export CFLAGS := --std=gnu99 -Wall -g -O2 -D TD_VERBOSE_DEBUG
#export CXXFLAGS := -Wall -g -O2 -D TD_VERBOSE_DEBUG
##############################

export PTHREAD := -pthread
export INCLUDES := -I$(INCLUDE_DIR) -L$(LIB_DIR)

SUBDIRS := OutSources \
		TangramBam \
		TangramScan \
		TangramMerge \
		TangramDetect \
		TangramIndex \
		TangramFilter

all:
	@echo "Building Tangram"
	@echo "========================================================="
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

clean:
	@rm -rf $(BIN_DIR) $(OBJ_DIR)
	@$(MAKE) clean --no-print-directory -C TangramScan
	@$(MAKE) clean --no-print-directory -C TangramMerge

clean_all:
	@rm -rf $(BIN_DIR) $(OBJ_DIR)
	@rm -rf $(LIB_DIR)
	@rm -rf $(INCLUDE_DIR)
	@$(MAKE) clean --no-print-directory -C TangramScan
	@$(MAKE) clean --no-print-directory -C TangramMerge
	@$(MAKE) clean --no-print-directory -C OutSources

.PHONY: all clean clean_all
