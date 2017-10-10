RM = rm -f
GCLIB := ../gclib
INC := -I$(GCLIB)
CC := g++
LDFLAGS =
BASEFLAGS := -std=c++11 $(INC) -Wall -Wextra -D_REENTRANT -fno-exceptions -fno-rtti
ifeq ($(findstring release,$(MAKECMDGOALS)),)
  CFLAGS := -g -DDEBUG $(BASEFLAGS)
  LDFLAGS += -g
else
  CFLAGS := -O2 -DNDEBUG $(BASEFLAGS)
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER  := g++
LIBS := 
OBJS := $(GCLIB)/GBase.o $(GCLIB)/gdna.o $(GCLIB)/codons.o \
 $(GCLIB)/GFaSeqGet.o $(GCLIB)/gff.o GIntervalTree.o $(GCLIB)/GArgs.o

.PHONY : all
all:    trmap

trmap: trmap.o $(OBJS)
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

.PHONY : clean
clean:: 
	@${RM} $(OBJS) trmap.o trmap
	@${RM} core.*


