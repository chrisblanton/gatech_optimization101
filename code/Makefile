FC=gfortran

mod_param.o: mod_param.f90
	$(FC) -c mod_param.f90
mod_io.o: mod_param.o mod_io.f90
	$(FC) -c mod_io.f90
mod_random.o: mod_param.o mod_random.f90
	$(FC) -c mod_random.f90
mod_matops.o: mod_param.o mod_random.o mod_matops.f90
	$(FC) -c mod_matops.f90
mod_solver.o: mod_param.o mod_random.o mod_matops.o mod_solver.f90
	$(FC) -c mod_solver.f90
mod_system_generator.o: mod_param.o mod_random.o mod_matops.o mod_system_generator.f90
	$(FC) -c mod_system_generator.f90
generate_systems.o: mod_param.o mod_random.o mod_io.o mod_matops.o mod_system_generator.o generate_systems.f90
	$(FC) -c generate_systems.f90
generate_systems.out: mod_param.o mod_random.o mod_io.o mod_matops.o mod_system_generator.o generate_systems.o
	$(FC) -o generate_systems.out  mod_param.o mod_random.o mod_io.o mod_matops.o mod_system_generator.o generate_systems.o
matmul.o: mod_param.o mod_matops.o matmul.f90
	$(FC) -c matmul.f90
matmul.out:  mod_param.o mod_random.o mod_matops.o matmul.o
	$(FC) -o matmul.out mod_param.o mod_random.o mod_matops.o matmul.o
test_solver.o: mod_param.o mod_random.o mod_io.o mod_matops.o mod_solver.o test_solver.f90
	$(FC) -c test_solver.f90
test_solver.out: mod_param.o mod_random.o mod_io.o mod_matops.o mod_solver.o test_solver.o
	$(FC) -o test_solver.out mod_param.o mod_random.o mod_io.o mod_matops.o mod_solver.o test_solver.o
benchmark_solver.o: mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver.f90
	$(FC) -c benchmark_solver.f90
benchmark_solver.out: mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver.o
	$(FC) -o benchmark_solver.out mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver.o
benchmark_solver_tilesize.o: mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver_tilesize.f90
	$(FC) -c benchmark_solver_tilesize.f90
benchmark_solver_tilesize.out: mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver_tilesize.o
	$(FC) -o benchmark_solver_tilesize.out mod_param.o mod_random.o mod_matops.o mod_solver.o mod_system_generator.o benchmark_solver_tilesize.o


clean:
	rm -f *.o *.mod *.out
