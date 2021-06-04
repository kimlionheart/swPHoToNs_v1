# this file is for general settings such as file list, etc.
# machine-specific settings such as include paths and #defines are in Makefile.local

##include Makefile.sw


CC   = mpicc
FC   = mpif90
SC  = sw5cc
#ARC  = ar
#LINK = $(CXX)

##CFLAGS +=  -fPIC -Wall -O2 -I$(SRCDIR)

#LFLAGS   +=   -fPIC -fopenmp

#INCLUDES += -I/home/qwang/develop/photoNs-2.1/inc
#LIBS     += -L/home/qwang/local/lib -lgsl -lgslcblas

INC2DFFT += -I./2decomp_fft/include
LIB2DFFT += -L./2decomp_fft/lib
#INC2DFFT += -I/home/export/online1/qwang/develop/2decomp_fft/include
#LIB2DFFT += -L/home/export/online1/qwang/develop/2decomp_fft/lib
SWLU      = -L/home/export/online3/swmore/opensource/swlu/lib -lswlu_mpi 
SWSPC     = -L/home/export/online3/swmore/release/lib -lspc

SRCDIR    = src
INCDIR	  = inc
OBJDIR    = obj
LIBDIR    = lib
#EXEDIR	  = ./run/	#bak
EXEDIR	  = ./
##EXEDIR	  = /home/export/online1/qwang
EXE       = photoNs-lcdm
##LIBS	 +=  -l2decomp_fft  
LIBS	 +=  -l2decomp_fft

OPTS	 += -OPT:IEEE_arith=2 -O2 -g
##OPTS	 += -DTABLE_GRAVITY
OPTS	 += -DLONGSHORT
OPTS	 += -DPERIODIC_CONDITION
OPTS	 += -DMYALLTOALL
#OPTS	 += -DPMTHREAD
##OPTS	 += -DMEMINIT
#OPTS	 += -DADAPTIVE

##OnAthrd	 = -DSWATHREAD
##OPTS	 += -DSWATHREAD
#OPTS	 += -lm
#OPTS	 += -DSWATHREAD
#OPTS	 += -lm_slave
#OPTS	 += -hybrid
#OPTS	 += -DPMONLY
##OPTS	 += -DP2POFF
##OPTS	 += -DDEMONSTRATION
##OPTS    += -DCHECKING
##OPTS    += -fopenmp

#OPTS	 += -Wno-unused-but-set-variable -Wno-unused-variable \
	    -Wno-unused-function -Wno-strict-aliasing -fopenmp

SOURCES	  = photoNs.c \
	    domains.c \
	    initial.c \
	    kernels.c \
	    remotes_task.c \
	    toptree.c \
	    utility.c \
	    adaptive.c \
	    operator.c \
	    snapshot.c \
	    partmesh.c \
	    icreater.c \
	    fmm.c 
#	    task.c 
KSLAVE    = slave_kernels.c get_reply.c 
	    
CONV2D	  = conv.f90

##TASKAT	  = task.c

OBJECTS	  = $(patsubst %.c, $(OBJDIR)/%.o,$(SOURCES))
OBJFORT	  = $(patsubst %.f90, $(OBJDIR)/%.o,$(CONV2D))
##OBJACC	  = $(patsubst %.c, $(OBJDIR)/%.o,$(TASKACC))
#OBJATH	  = $(OBJDIR)/task.o
#OBJATS 	  = $(OBJDIR)/task_slave.o

OBJSLAVE =  $(patsubst %.c, $(OBJDIR)/%.o,$(KSLAVE))

LINK_OBJ =  slave_kernels.o get_reply.o

LINK_ONJS = $(patsubst %.o, $(OBJDIR)/%.o,$(LINK_OBJ))


exe: $(OBJECTS) $(OBJFORT) $(LINK_ONJS) 
	@mkdir -p $(EXEDIR)
#	$(FC) -g -OPT:IEEE_arith=2 -O2 -hybrid  -Wl,--whole-archive,-wrap,athread_init,-wrap,__expt_handler,-wrap,__real_athread_spawn /home/export/online3/swmore/release/lib/libspc.a -Wl,--no-whole-archive $(OBJDIR)/*.o $(LIBS) -lm_slave  libgptl.a $(LIB2DFFT) $(SWLU) $(SWSPC) -o $(EXEDIR)/$(EXE)
	$(FC) -g -OPT:IEEE_arith=2 -O2 -hybrid  -Wl,--whole-archive,-wrap,athread_init,-wrap,__expt_handler,-wrap,__real_athread_spawn /home/export/online3/swmore/release/lib/libspc.a -Wl,--no-whole-archive $(OBJDIR)/*.o $(LIBS) -lm_slave  libgptl.a $(LIB2DFFT) $(SWLU) $(SWSPC) -o ./run/$(EXE)

$(OBJSLAVE): $(OBJDIR)/%.o: $(SRCDIR)/%.c
	sw5cc -c -slave -O2 -msimd -c -o "$@" "$<"
	
$(OBJDIR)/%.o:  $(SRCDIR)/%.c 
	@mkdir -p $(OBJDIR)
	$(CC) -c -g   $(CFLAGS) $(OPTS) $(INCLUDES) $(OnAthrd)   -I$(INCDIR) $(SWLU) $(SWSPC)  -o "$@" "$<"
$(OBJDIR)/%.o:  $(SRCDIR)/%.f90 
##	$(FC) -c  $(INC2DFFT) -o "$@" "$<"
	$(FC) -c -g -OPT:IEEE_arith=2   $(INC2DFFT) -o "$@" "$<"

#$(OBJDIR)/task.o:  $(SRCDIR)/task.c Makefile.sw
#	$(SC) -c -host  -O2 -OPT:fast_math_library  -OPT:IEEE_arith=2 -lm $(OnAthrd) $(CFLAGS) $(INCLUDES) -I$(INCDIR) -o "$@" "$<"

#$(OBJDIR)/task_slave.o:  $(SRCDIR)/task_slave.c Makefile.sw
#	$(SC) -c -slave -O2  -OPT:IEEE_arith=2  $(OnAthrd) $(CFLAGS) -lm_slave $(INCLUDES) -I$(INCDIR) -o "$@" "$<"

.PHONY: demo
demo:
	cd $(EXEDIR); \
	bsub -np 1 ./$(EXE) ../demo/lcdm.run

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(EXEDIR)/$(EXE)
	cp ./math_obj/erf_slave.o ./$(OBJDIR)
	cp ./math_obj/exp_slave.o ./$(OBJDIR) 
	
	



