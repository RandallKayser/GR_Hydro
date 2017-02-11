SRCS_MAIN = c_files/main.c c_files/write_to_file.c c_files/boundary.c c_files/timestep.c c_files/misc.c c_files/update.c 
SRCS_MINK = c_files/init_files/minkowski_coords.c c_files/init_files/minkowski_sod_shock.c

HEAD = headers/write_to_file.h headers/boundary.h headers/constants.h headers/geometry.h headers/initfuncs.h headers/misc.h headers/timestep.h headers/update.h

default: hydro_mink_1d_sod

hydro_mink_1d_sod: $(SRCS_MAIN) $(SRCS_MINK) $(HEAD) 
	gcc $(SRCS_MAIN) $(SRCS_MINK) -o hydro_mink_1d_sod

clean:
	rm -f *.hydro