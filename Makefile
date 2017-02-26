SRCS_MAIN = c_files/main.c c_files/write_to_file.c c_files/boundary.c c_files/timestep.c c_files/misc.c c_files/update.c 

SRCS_MINK_CART = c_files/init_files/minkowski_coords.c
MINK_SOD_INIT_1D = c_files/init_files/minkowski_sod_shock_1d.c
MINK_SKEW_SOD_INIT_2D = c_files/init_files/minkowski_skew_sod_2d.c

HEAD = headers/write_to_file.h headers/boundary.h headers/constants.h headers/geometry.h headers/initfuncs.h headers/misc.h headers/timestep.h headers/update.h

default: hydro_mink_2d_skew_sod

hydro_mink_1d_sod: $(SRCS_MAIN) $(SRCS_MINK_CART) $(MINK_SOD_INIT_1D) $(HEAD) 
	gcc -Wall -O3 $(SRCS_MAIN) $(SRCS_MINK_CART) $(MINK_SOD_INIT_1D) -o hydro_mink_1d_sod

hydro_mink_2d_skew_sod: $(SRC_MAIN) $(SRCS_MINK_CART) $(MINK_SKEW_SOD_INIT_2D) $(HEAD)
	gcc $(SRCS_MAIN) $(SRCS_MINK_CART) $(MINK_SKEW_SOD_INIT_2D) -o hydro_mink_2d_skew_sod

hydro_mink_1d_sod_gdb: $(SRCS_MAIN) $(SRCS_MINK_CART) $(MINK_SOD_INIT_1D) $(HEAD) 
	gcc -g $(SRCS_MAIN) $(SRCS_MINK_CART) $(MINK_SOD_INIT_1D) -o hydro_mink_1d_sod

hydro_mink_2d_skew_sod_gdb: $(SRC_MAIN) $(SRCS_MINK_CART) $(MINK_SKEW_SOD_INIT_2D) $(HEAD)
	gcc -g $(SRCS_MAIN) $(SRCS_MINK_CART) $(MINK_SKEW_SOD_INIT_2D) -o hydro_mink_2d_skew_sod

clean:
	rm -f hydro_mink_1d_sod hydro_mink_2d_skew_sod