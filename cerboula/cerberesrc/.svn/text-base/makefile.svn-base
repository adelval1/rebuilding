CERBEREHOME=$(shell pwd)
OBJ=$(CERBEREHOME)/../cerbereobj
SRC=$(CERBEREHOME)/../cerberesrc
OBJ=$(CERBEREHOME)/../cerbereobj
EXE=$(CERBEREHOME)/../exe/cerbere.exe
PEGDIR = $(CERBEREHOME)/../../libraries/tempegase45/exec/lib


# C COMPILER GIVEN BELOW
CC=cc

# FORTRAN-77 compiler given below
FC=f77
# FORTRAN-90 compiler given below
FC=ifort

# LIBRARY NAMES GIVEN BELOW
LIBRARIES=-lm
f90lib=
LIBRARIES=  -L$(PEGDIR)  
Peglib = -lpegase -lm 



# LINKER GIVEN BELOW
LD=ifort

# LINKER FLAGS GIVEN BELOW
LDFLAGS= -O
#LDFLAGS= -g

# DEFINES
AR=ar
ARFLAGS=rv
CP=cp
RM=rm -f

# INCLUDE SEARCH PATH GIVEN BELOW
INCLUDES= $(CERBEREHOME)/../../libraries/tempegase45/src/*mod_* $(CERBEREHOME)/../../libraries/tempegase45/src/*interf*  $(CERBEREHOME)/../../libraries/tempegase45/exec/lib


# COMPILER FLAGS GIVEN BELOW
COMPILERFLAGS= -O
#COMPILERFLAGS= -g
CFLAGS=$(COMPILERFLAGS) 
#-ggdb 
FFLAGS=$(COMPILERFLAGS) $(INCLUDES) #-B108 -YEXT_NAMES=LCS -YCFRL=1 -trap=INVALID,DIVBYZERO -Rc -lg2c -ladebug -g -C
PFLAGS=$(COMPILERFLAGS)
COBFLAGS=$(COMPILERFLAGS) -xref
CXXFLAGS=$(COMPILERFLAGS) -xref


.SUFFIXES: .o .f90

# SOURCES GIVEN BELOW
SRCS=  $(SRC)/modules.f90 \
 $(SRC)/interf.f90\
 $(SRC)/CERBERE.f90\
 $(SRC)/input_reading.f90\
 $(SRC)/enth_rebuilding.f90\
 $(SRC)/chemar.f90\
 $(SRC)/heatflux.f90\
 $(SRC)/initialcond.f90\
 $(SRC)/pegasemodule.90\
 $(SRC)/VKI_BLinput.f90\
 $(SRC)/cat_rebuilding.f90\
 $(SRC)/heat_identif.f90\
 $(SRC)/Kp_rebuilding.f90\
 $(SRC)/write_sol.f90\
 $(SRC)/set_neboula_input.f90\
 $(SRC)/equilcondition.f90\
 $(SRC)/initialdistrib.f90\
 $(SRC)/propequil.f90\
 $(SRC)/propfroz.f90\
 $(SRC)/abacus.f90\
 $(SRC)/plot.f90\
 $(SRC)/init_dist.f90\
 $(SRC)/data_storage.f90\
 $(SRC)/richardson_extrapolation.f90\
 $(SRC)/dwdrho.f90

# OBJECTS GIVEN BELOW
OBJS=  $(OBJ)/modules.o \
 $(OBJ)/interf.o\
 $(OBJ)/CERBERE.o\
 $(OBJ)/input_reading.o\
 $(OBJ)/enth_rebuilding.o\
 $(OBJ)/chemar.o\
 $(OBJ)/heatflux.o\
 $(OBJ)/initialcond.o\
 $(OBJ)/pegasemodule.o\
 $(OBJ)/VKI_BLinput.o\
 $(OBJ)/cat_rebuilding.o\
 $(OBJ)/heat_identif.o\
 $(OBJ)/Kp_rebuilding.o\
 $(OBJ)/write_sol.o\
 $(OBJ)/set_neboula_input.o\
 $(OBJ)/equilcondition.o\
 $(OBJ)/initialdistrib.o\
 $(OBJ)/propequil.o\
 $(OBJ)/propfroz.o\
 $(OBJ)/abacus.o\
 $(OBJ)/plot.o\
 $(OBJ)/init_dist.o\
 $(OBJ)/data_storage.o\
 $(OBJ)/richardson_extrapolation.o\
 $(OBJ)/dwdrho.o

./cerbere.exe : $(OBJS)
	@echo Linking... 
	$(LD) -o  $(EXE) $(LDFLAGS) $(OBJS) $(LIBRARIES) $(Peglib) $(f90lib)
	@echo Success: $(EXE) has been created!
	@echo the rebuilding procedure is all set and ready to go!
	@echo '\07'


.f90.o:
	$(FC) -g -c   -I$(INCLUDES) $< $(FFLAGS) 

.f90.a:
	$(FC90) -c $(FFLAGS) $(INCLUDES) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o
 
$(OBJ)/modules.o : $(SRC)/modules.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/modules.f90
	
$(OBJ)/interf.o : $(SRC)/interf.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/interf.f90

$(OBJ)/CERBERE.o : $(SRC)/CERBERE.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/CERBERE.f90

$(OBJ)/input_reading.o : $(SRC)/input_reading.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/input_reading.f90

$(OBJ)/enth_rebuilding.o : $(SRC)/enth_rebuilding.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/enth_rebuilding.f90

$(OBJ)/chemar.o : $(SRC)/chemar.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/chemar.f90

$(OBJ)/heatflux.o : $(SRC)/heatflux.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/heatflux.f90

$(OBJ)/initialcond.o : $(SRC)/initialcond.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/initialcond.f90

$(OBJ)/pegasemodule.o : $(SRC)/pegasemodule.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/pegasemodule.f90

$(OBJ)/VKI_BLinput.o : $(SRC)/VKI_BLinput.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/VKI_BLinput.f90

$(OBJ)/cat_rebuilding.o : $(SRC)/cat_rebuilding.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/cat_rebuilding.f90

$(OBJ)/heat_identif.o : $(SRC)/heat_identif.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/heat_identif.f90

$(OBJ)/Kp_rebuilding.o : $(SRC)/Kp_rebuilding.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/Kp_rebuilding.f90

$(OBJ)/write_sol.o : $(SRC)/write_sol.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/write_sol.f90

$(OBJ)/set_neboula_input.o : $(SRC)/set_neboula_input.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/set_neboula_input.f90

$(OBJ)/equilcondition.o : $(SRC)/equilcondition.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/equilcondition.f90

$(OBJ)/initialdistrib.o : $(SRC)/initialdistrib.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/initialdistrib.f90

$(OBJ)/propequil.o : $(SRC)/propequil.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/propequil.f90

$(OBJ)/propfroz.o : $(SRC)/propfroz.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/propfroz.f90

$(OBJ)/abacus.o : $(SRC)/abacus.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/abacus.f90

$(OBJ)/plot.o : $(SRC)/plot.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/plot.f90

$(OBJ)/init_dist.o : $(SRC)/init_dist.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/init_dist.f90

$(OBJ)/data_storage.o : $(SRC)/data_storage.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/data_storage.f90

$(OBJ)/richardson_extrapolation.o : $(SRC)/richardson_extrapolation.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/richardson_extrapolation.f90

$(OBJ)/dwdrho.o : $(SRC)/dwdrho.f90
	@echo compiling object
	@cd $(OBJ); $(FC) $(FFLAGS)  -c $(SRC)/dwdrho.f90

clean:
	$(RM) $(OBJS) core $(EXE) ../cerbereobj/*

touchsrcs:
	@touch $(SRCS)

variables:
	@echo "FORTRAN FLAGS" $(FC) $(FFLAGS) 
