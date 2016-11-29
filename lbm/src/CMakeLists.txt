cmake_minimum_required(VERSION 3.5)

set(PATH_MATH "math")
set(PATH_AREA "modeling_area")
set(PATH_PHYS "phys_values")
set(PATH_SOLV "solver")
set(PATH_BC "solver/bc")

set(
	source_list
	"math/array_func_impl.h"
	"math/my_matrix_2d.h"
	"math/my_matrix_2d_impl.h"
	"math/my_matrix_interface.h"
	"modeling_area/fluid.h"
	"modeling_area/medium.h"
	"phys_values/distribution_func.h"
	"phys_values/distribution_func_impl.h"
	"phys_values/macroscopic_param.h"
	"phys_values/macroscopic_param_impl.h"
	"solver/solver.h"
	"solver/srt.h"
	"solver/bc/bc.h"
	"modeling_area/fluid.cpp"
	"modeling_area/medium.cpp"
	"solver/srt.cpp"
	"solver/bc/bc.cpp"
	"main.cpp"
)

add_executable(${PROJECT_NAME} ${source_list})

foreach(source IN LISTS source_list)
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" source_path_msvc "${source_path}")
    source_group("${source_path_msvc}" FILES "${source}")
endforeach()
