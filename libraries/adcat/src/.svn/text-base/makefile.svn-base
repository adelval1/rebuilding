CATCHOME=$(shell pwd)
DIR=$(CATCHOME)/..
SRC=$(DIR)/src
OBJ=$(DIR)/obj
LIB=$(DIR)/exe
EXE=$(DIR)/exe/libCATC.a

# FORTRAN compiler
CC=gcc

#PATH TO INDLUDE FILES
INCLUDES=$(SRC)/firstlayer/secondlayer/thirdlayer/

# COMPILER FLAGS GIVEN BELOW
#Compilerflags=-O4 -I$(INCLUDES)
Compilerflags=-g -O0 -I$(INCLUDES)

# DEFINITIONS
AR=ar
ARFLAGS=ruv
RM=rm -f
RNLIB = touch

OBJS=\
 $(OBJ)/IsoCatWall.o\
 $(OBJ)/AdiaCatWall.o\
 $(OBJ)/REWCatWall.o\
 $(OBJ)/Catalycity.o\
 $(OBJ)/AirBetaCristobaliteBalat.o\
 $(OBJ)/AirQuartzBalat.o\
 $(OBJ)/AirPM1000VKI.o\
 $(OBJ)/AirCSiCVKI.o\
 $(OBJ)/AirPM1000Preci.o\
 $(OBJ)/AirSiO2Fertig.o\
 $(OBJ)/AirSiCFertig.o\
 $(OBJ)/AirSiO2DetailedFertig.o\
 $(OBJ)/AirSiCDetailedFertig.o\
 $(OBJ)/AirPtBarbato.o\
 $(OBJ)/AirNiBarbato.o\
 $(OBJ)/AirCuBarbato.o\
 $(OBJ)/GenericReactionRateBased.o\
 $(OBJ)/DRW.o\
 $(OBJ)/ZGB.o\
 $(OBJ)/ZGBplus.o\
 $(OBJ)/ZGBplusL.o\
 $(OBJ)/EmissivityPM1000Buursink.o\
 $(OBJ)/CatalycityMoleProduction.o\
 $(OBJ)/numerics.o\
 $(OBJ)/les.o\

default: CATClib

$(OBJ)/IsoCatWall.o : $(SRC)/firstlayer/IsoCatWall.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/IsoCatWall.c
$(OBJ)/AdiaCatWall.o : $(SRC)/firstlayer/AdiaCatWall.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/AdiaCatWall.c
$(OBJ)/REWCatWall.o : $(SRC)/firstlayer/REWCatWall.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/REWCatWall.c

$(OBJ)/Catalycity.o : $(SRC)/firstlayer/secondlayer/Catalycity.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/Catalycity.c
$(OBJ)/AirBetaCristobaliteBalat.o : $(SRC)/firstlayer/secondlayer/AirBetaCristobaliteBalat.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirBetaCristobaliteBalat.c
$(OBJ)/AirQuartzBalat.o : $(SRC)/firstlayer/secondlayer/AirQuartzBalat.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirQuartzBalat.c
$(OBJ)/AirPM1000VKI.o : $(SRC)/firstlayer/secondlayer/AirPM1000VKI.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirPM1000VKI.c
$(OBJ)/AirCSiCVKI.o : $(SRC)/firstlayer/secondlayer/AirCSiCVKI.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirCSiCVKI.c
$(OBJ)/AirPM1000Preci.o : $(SRC)/firstlayer/secondlayer/AirPM1000Preci.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirPM1000Preci.c
$(OBJ)/AirSiO2Fertig.o : $(SRC)/firstlayer/secondlayer/AirSiO2Fertig.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirSiO2Fertig.c
$(OBJ)/AirSiCFertig.o : $(SRC)/firstlayer/secondlayer/AirSiCFertig.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirSiCFertig.c
$(OBJ)/AirSiO2DetailedFertig.o : $(SRC)/firstlayer/secondlayer/AirSiO2DetailedFertig.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirSiO2DetailedFertig.c
$(OBJ)/AirPtBarbato.o : $(SRC)/firstlayer/secondlayer/AirPtBarbato.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirPtBarbato.c
$(OBJ)/AirNiBarbato.o : $(SRC)/firstlayer/secondlayer/AirNiBarbato.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirNiBarbato.c
$(OBJ)/AirCuBarbato.o : $(SRC)/firstlayer/secondlayer/AirCuBarbato.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirCuBarbato.c
$(OBJ)/AirSiCDetailedFertig.o : $(SRC)/firstlayer/secondlayer/AirSiCDetailedFertig.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/AirSiCDetailedFertig.c
$(OBJ)/GenericReactionRateBased.o : $(SRC)/firstlayer/secondlayer/GenericReactionRateBased.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/GenericReactionRateBased.c
$(OBJ)/DRW.o : $(SRC)/firstlayer/secondlayer/DRW.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/DRW.c
$(OBJ)/ZGB.o : $(SRC)/firstlayer/secondlayer/ZGB.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/ZGB.c
$(OBJ)/ZGBplus.o : $(SRC)/firstlayer/secondlayer/ZGBplus.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/ZGBplus.c
$(OBJ)/ZGBplusL.o : $(SRC)/firstlayer/secondlayer/ZGBplusL.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/ZGBplusL.c
$(OBJ)/EmissivityPM1000Buursink.o : $(SRC)/firstlayer/secondlayer/EmissivityPM1000Buursink.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/EmissivityPM1000Buursink.c


$(OBJ)/CatalycityMoleProduction.o : $(SRC)/firstlayer/secondlayer/thirdlayer/CatalycityMoleProduction.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/thirdlayer/CatalycityMoleProduction.c
$(OBJ)/numerics.o : $(SRC)/firstlayer/secondlayer/thirdlayer/numerics.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/thirdlayer/numerics.c
$(OBJ)/les.o : $(SRC)/firstlayer/secondlayer/thirdlayer/les.c
	@cd $(OBJ);$(CC) $(Compilerflags) -c $(SRC)/firstlayer/secondlayer/thirdlayer/les.c


CATClib: $(OBJS)
	@echo "Building library CATC"
	@$(RM) $(EXE)
	@cd $(LIB);
	@$(AR) $(ARFLAGS) $(EXE) $(OBJS)
	@cd $(LIB);
	@$(RNLIB) $(EXE)

clean:	
	@echo deleting objects
	@cd $(OBJ);
	@$(RM)  $(OBJS) 
	@echo deleting library
	@$(RM) $(EXE)
