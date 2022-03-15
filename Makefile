#Compilador Tupa
#CC=ftn 

#Compilador 
FC = mpif90

# Opcoes de compilacao
FFLAGS = -c -O3 -g -w

# Opcoes de otimizacao
FFLAGSOPT = -O3 -g -w

VPATH = src
MODDIR = mod
BUILDDIR = build

# Arquivos objeto do annMPCA
SRCMPCA := $(BUILDDIR)/foul.o \
$(BUILDDIR)/uniformR8.o \
$(BUILDDIR)/newTypes.o \
$(BUILDDIR)/normalR8.o \
$(BUILDDIR)/rnaFunctions.o \
$(BUILDDIR)/annGeneralization.o \
$(BUILDDIR)/annTraining.o \
$(BUILDDIR)/mpcaFunctions.o \
$(BUILDDIR)/mpca.o

# Arquivos objeto do annTest
SRCTEST := $(BUILDDIR)/foul.o \
$(BUILDDIR)/rnaFunctions.o \
$(BUILDDIR)/newTypes.o \
$(BUILDDIR)/annTest.o \
$(BUILDDIR)/main_test.o

all: 	$(BUILDDIR)/foul.o \
	$(BUILDDIR)/newTypes.o \
	$(BUILDDIR)/uniformR8.o \
	$(BUILDDIR)/normalR8.o \
	$(BUILDDIR)/rnaFunctions.o \
	$(BUILDDIR)/annGeneralization.o \
	$(BUILDDIR)/annTraining.o \
	$(BUILDDIR)/mpcaFunctions.o \
	$(BUILDDIR)/mpca.o \
	$(BUILDDIR)/annTest.o \
	$(BUILDDIR)/main_generalization.o \
	$(BUILDDIR)/main_test.o \
	annMPCA \
	annTest

annMPCA:
	$(FC) $(FFLAGSOPT) -o annMPCA $(SRCMPCA)

annTest:
	$(FC) $(FFLAGSOPT) -o annTest $(SRCTEST)

$(BUILDDIR)/%.o: $(VPATH)/%.f90
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) $< -o $@

clean:	removemod
	rm -rf *.*~ Makefile~ build/*.o *.mod annTest annTest annMPCA

removemod:
	rm -f build/*.o *.mod
