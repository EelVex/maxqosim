WARN=-Wall -Wextra -fmax-errors=1 #-g -fcheck=all -fbacktrace
#-fstack-arrays -fPIC 
OPTIMIZATIONS=-Ofast -ffast-math -fno-range-check -funroll-loops  -fstack-arrays
OPTIONS=-ffree-line-length-512
FORTRAN=gfortran

OBJDIR=_build/
BUILDDIR=_build/

%.mod: %.F90
	@mkdir -p $(OBJDIR)
	@$(FORTRAN) -c $(WARN) $(OPTIMIZATIONS) $< -o $(OBJDIR)$*.o

main:  main.F90 integrators.F90 parameters.F90 equations.F90
	@echo "> compiling $@"
	@# @$(FORTRAN) -J$(MODDIR) $(WARN) $(OPTIMIZATIONS) -fwhole-program $(OPENMP) $< -o $@  $(LIB) -llapack -lblas
	@$(FORTRAN) $(WARN) $(OPTIONS) $(OPTIMIZATIONS) -flto main.F90 integrators.F90 -o $@ 

clean:
	rm -rf _build/
	rm *.mod

# EOF
