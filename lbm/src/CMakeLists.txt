cmake_minimum_required(VERSION 3.5)

set(PATH_MATH "math")
set(PATH_AREA "modeling_area")
set(PATH_PHYS "phys_values")
set(PATH_SOLV "solver")
set(PATH_BC "solver/bc")

set(
	source_list
	"math/array_func_impl.h"
	"math/2d/my_matrix_2d.h"
	"math/2d/my_matrix_2d_impl.h"
	"math/3d/my_matrix_3d.h"
	"math/3d/my_matrix_3d_impl.h"
	"math/my_matrix_interface.h"
	"modeling_area/fluid.h"
	"modeling_area/medium.h"
	"phys_values/2d/distribution_func_2d.h"
	"phys_values/2d/distribution_func_2d_impl.h"
	"phys_values/2d/macroscopic_param_2d.h"
	"phys_values/2d/macroscopic_param_2d_impl.h"
	"phys_values/3d/distribution_func_3d.h"
	"phys_values/3d/distribution_func_3d_impl.h"
	"phys_values/3d/macroscopic_param_3d.h"
	"phys_values/3d/macroscopic_param_3d_impl.h"
	"phys_values/distribution_function_interface.h"
	"solver/solver.h"
	"solver/srt.h"
	"solver/bc/bc.h"
	"modeling_area/fluid.cpp"
	"modeling_area/medium.cpp"
	"solver/im_body/immersed_body.h"
	"solver/im_body/immersed_body.cpp"
	"solver/ib_srt.h"
	"solver/ib_srt.cpp"
	"solver/srt.cpp"
	"solver/bc/bc.cpp"
	"solver/mrt.h"
	"solver/mrt.cpp"
	"main.cpp"
)

add_executable(${PROJECT_NAME} ${source_list})

foreach(source IN LISTS source_list)
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" source_path_msvc "${source_path}")
    source_group("${source_path_msvc}" FILES "${source}")
endforeach()
