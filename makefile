FC = mpifort

OPT = -O3
FCFLAGS = -g -ffree-line-length-0 $(OPT)

SOURCEDIR = src
OBJDIR = obj

SRCS = comm_variable.f90  read_setup.f90 \
	 MPI_setup.f90 wavenumber.f90 DefineGrid.f90 \
	 boundaries.f90 RHS.f90  analytical_sol.f90 \
	 Operators.f90 Smoother.f90 \
     CSLP_Solver.f90 deflation_setup.f90 \
	 user_module.f90 idrs_module.f90 \
     solvers.f90 Write_data.f90 \
     main.f90

OBJS = $(patsubst %.f90,$(OBJDIR)/%.o,$(SRCS))

EXECUTABLE = helmholtz_2d_solver.out

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SOURCEDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf *.mod $(OBJDIR) $(EXECUTABLE)
