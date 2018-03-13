# Automatic Makefile

FCM_BLD_CFG = /home/amber/NEMOGCM/TOOLS/COMPILE/bld_tools.cfg
export FCM_VERBOSE ?= 1

export FCM_ROOTDIR := /home/amber/NEMOGCM/TOOLS/REBUILD_NEMO/BLD
export FCM_CACHEDIR := $(FCM_ROOTDIR)/.cache/.bld
export FCM_CFGDIR := $(FCM_ROOTDIR)/cfg
export FCM_SRCDIR := $(FCM_ROOTDIR)/src
export FCM_BINDIR := $(FCM_ROOTDIR)/bin
export FCM_ETCDIR := $(FCM_ROOTDIR)/etc
export FCM_LIBDIR := $(FCM_ROOTDIR)/lib
export FCM_TMPDIR := $(FCM_ROOTDIR)/tmp
export FCM_DONEDIR := $(FCM_ROOTDIR)/done
export FCM_FLAGSDIR := $(FCM_ROOTDIR)/flags
export FCM_INCDIR := $(FCM_ROOTDIR)/inc
export FCM_PPSRCDIR := $(FCM_ROOTDIR)/ppsrc
export FCM_OBJDIR := $(FCM_ROOTDIR)/obj
export FCM_ROOTPATH := $(FCM_ROOTDIR)
export FCM_CACHEPATH := $(FCM_CACHEDIR)
export FCM_CFGPATH := $(FCM_CFGDIR)
export FCM_SRCPATH := $(FCM_SRCDIR)
export FCM_BINPATH := $(FCM_BINDIR)
export FCM_ETCPATH := $(FCM_ETCDIR)
export FCM_LIBPATH := $(FCM_LIBDIR)
export FCM_TMPPATH := $(FCM_TMPDIR)
export FCM_DONEPATH := $(FCM_DONEDIR)
export FCM_FLAGSPATH := $(FCM_FLAGSDIR)
export FCM_INCPATH := $(FCM_INCDIR)
export FCM_PPSRCPATH := $(FCM_PPSRCDIR)
export FCM_OBJPATH := $(FCM_OBJDIR)

export PERL5LIB := /home/amber/NEMOGCM/EXTERNAL/fcm/bin/../lib:$(PERL5LIB)

export AR = ar
export ARFLAGS = rcs
export CC = cc
export CC_COMPILE = -c
export CC_DEFINE = -D
export CC_INCLUDE = -I
export CC_OUTPUT = -o
export CFLAGS = 
export CPPKEYS = 
export DIFF3FLAGS = -E -m
export FC = mpif90
export FC_COMPILE = -c
export FC_DEFINE = -D
export FC_INCLUDE = -I
export FC_MODSEARCH = 
export FC_OUTPUT = -o
export FFLAGS = -c -fpp -r8 -O3 -assume byterecl -heap-arrays -diag-disable 10010 -xSSE4.1 -axSSE4.2 -I/home/amber/XIOS/inc
export FPPKEYS = 
export LD = mpif90
export LDFLAGS = -lstdc++ -shared-intel -diag-disable 10010 -L/home/amber/XIOS/lib -lxios -lnetcdf -lnetcdff -lpnetcdf -lhdf5 -lhdf5_hl
export LD_LIBLINK = -l
export LD_LIBSEARCH = -L
export LD_OUTPUT = -o
export SHELL = /bin/sh

vpath % $(FCM_BINPATH)
vpath %.done $(FCM_DONEPATH)
vpath %.idone $(FCM_DONEPATH)
vpath %.etc $(FCM_ETCPATH)
vpath %.flags $(FCM_FLAGSPATH)
vpath % $(FCM_INCPATH)
vpath %.a $(FCM_LIBPATH)
vpath %.o $(FCM_OBJPATH)

FCM_BLD_TARGETS = rebuild_nemo.exe

.PHONY : all

all : $(FCM_BLD_TARGETS)

$(FCM_DONEDIR)/FCM_CP.dummy:
	touch $@

CFLAGS__nemo.flags: CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__nemo.flags: CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__nemo.flags: FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__nemo.flags: FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__nemo.flags: LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

LD__nemo.flags: LD.flags
	touch $(FCM_FLAGSDIR)/$@

CC.flags:
	touch $(FCM_FLAGSDIR)/$@

CFLAGS.flags: CC.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS.flags:
	touch $(FCM_FLAGSDIR)/$@

FC.flags:
	touch $(FCM_FLAGSDIR)/$@

FFLAGS.flags: FC.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS.flags:
	touch $(FCM_FLAGSDIR)/$@

LD.flags:
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS.flags:
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__nemo__rebuild_nemo.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__nemo__rebuild_nemo.flags: LDFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

LD__nemo__rebuild_nemo.flags: LD__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

rebuild_nemo.exe: rebuild_nemo.o LD__nemo__rebuild_nemo.flags LDFLAGS__nemo__rebuild_nemo.flags
	fcm_internal load:F nemo $< $@

rebuild_nemo.o: /home/amber/NEMOGCM/TOOLS/REBUILD_NEMO/src/rebuild_nemo.f90 FFLAGS__nemo__rebuild_nemo.flags
	fcm_internal compile:F nemo $< $@

# EOF
