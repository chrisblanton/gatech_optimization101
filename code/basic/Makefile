
# Compiler Options
FC=gfortran
PROFILE=
#PROFILE=-g -pg
OPT=-O1

# Folders
SRC=src
BUILD=build
BIN=bin

# Variables
VPATH=src:build


# Rules

mod_param.o: mod_param.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/mod_param.f90 -J$(BUILD) -o $(BUILD)/$@
mod_io.o: mod_param.o mod_io.f90
	$(FC) $(PROFILE)  $(OPT) -c $(SRC)/mod_io.f90 -J$(BUILD) -o $(BUILD)/$@
mod_random.o: mod_param.o mod_random.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/mod_random.f90 -J$(BUILD) -o $(BUILD)/$@
mod_matops.o: mod_param.o mod_random.o mod_matops.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/mod_matops.f90 -J$(BUILD) -o $(BUILD)/$@
mod_solver.o: mod_param.o mod_random.o mod_matops.o mod_solver.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/mod_solver.f90 -J$(BUILD) -o $(BUILD)/$@
mod_system_generator.o: mod_param.o mod_random.o mod_matops.o mod_system_generator.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/mod_system_generator.f90 -J$(BUILD) -o $(BUILD)/$@
matmul.o: mod_param.o mod_random.o mod_matops.o matmul.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/matmul.f90 -J$(BUILD) -o $(BUILD)/$@
matmul.out: mod_param.o mod_random.o mod_matops.o matmul.o
	$(FC) $(PROFILE) $(OPT) -o $(BIN)/matmul.out $(addprefix $(BUILD)/, mod_param.o  mod_random.o mod_matops.o matmul.o)
generate_systems.o: mod_param.o mod_random.o mod_io.o mod_matops.o mod_system_generator.o generate_systems.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/generate_systems.f90 -J$(BUILD) -o $(BUILD)/$@
generate_systems.out: mod_param.o mod_random.o mod_io.o mod_matops.o mod_system_generator.o generate_systems.o
	$(FC) $(PROFILE) $(OPT) -o $(BIN)/$@ $(addprefix $(BUILD)/, mod_param.o mod_random.o mod_io.o mod_matops.o mod_system_generator.o generate_systems.o)
test_solver.o: mod_param.o mod_random.o mod_io.o mod_matops.o mod_solver.o test_solver.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/test_solver.f90 -J$(BUILD) -o $(BUILD)/$@
Test_solver.out: mod_param.o mod_random.o mod_io.o mod_matops.o mod_solver.o test_solver.o
	$(FC) $(PROFILE) $(OPT) -o $(BIN)/$@ $(addprefix $(BUILD)/, mod_param.o mod_random.o mod_io.o mod_matops.o mod_solver.o test_solver.o)
benchmark_solver.o: mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver.f90
	$(FC) $(PROFILE) $(OPT) -c $(SRC)/benchmark_solver.f90 -J$(BUILD) -o $(BUILD)/$@
benchmark_solver.out: mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver.o
	$(FC) $(PROFILE) $(OPT) -o $(BIN)/$@ $(addprefix $(BUILD)/, mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver.o)

clean:
	rm -f $(BUILD)/*.o $(BUILD)/*.mod $(BIN)/*.out
