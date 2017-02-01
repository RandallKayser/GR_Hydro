SRCS = c_files/main.c c_files/write_to_file.c c_files/boundary.c c_files/timestep.c c_files/misc.c c_files/update.c c_files/init_files/minkowski_coords.c c_files/init_files/minkowski_static_init.c

HEAD = headers/write_to_file.h headers/boundary.h headers/constants.h headers/geometry.h headers/initfuncs.h headers/misc.h headers/timestep.h headers/update.h

default: HYDRO

HYDRO: $(SRCS) $(HEAD) 
	gcc $(SRCS) -o HYDRO

clean:
	rm -f HYDRO 