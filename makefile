FF = mpif90
FFFlags = -ffree-form -ffree-line-length-none -O3  -fdefault-real-8 -cpp -I$(PETSC_DIR)/include
LDFlags = -lpetsc -L$(PETSC_DIR)/lib

OS := $(shell uname)

ifeq ($(OS), Linux)
   LDFlags += -Wl,-rpath=$(PETSC_DIR)/lib
endif

TARGET = euler

all : euler

srcfiles = $(wildcard *.F90)    #Find all the files with F90 as extension
OBJ  = $(patsubst %.F90,%.o,$(srcfiles))

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

euler : $(OBJ)
	@$(FF) $(FFFlags) $(OBJ) -o $(TARGET) $(LDFlags)

comvar.mod : comvar.F90
	@$(FF) $(FFFlags) -c comvar.F90 $(LDFlags)

%.o: %.F90 comvar.mod
	@$(FF) $(FFFlags) -c $*.F90 $(LDFlags)

cleano::
	@rm -f $(OBJ) *.mod
clean::	
	@rm -f $(TARGET) *.plt *.xml Output/*.xml Output/*.plt
