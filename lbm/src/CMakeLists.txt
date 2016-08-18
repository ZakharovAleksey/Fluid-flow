cmake_minimum_required(VERSION 3.5)

set(PATH_MATH "math")
set(PATH_AREA "modeling_area")
set(PATH_PHYS "phys_values")
set(PATH_SOLV "solver")
set(PATH_BC "solver/bc")


set(HEADER_LIST 
	${PATH_MATH}/array_func_impl.h
	${PATH_MATH}/my_matrix.h
	${PATH_MATH}/my_matrix_impl.h

	${PATH_AREA}/fluid.h
	${PATH_AREA}/medium.h

	${PATH_PHYS}/distribution_func.h
	${PATH_PHYS}/distribution_func_impl.h
	${PATH_PHYS}/macroscopic_param.h
	${PATH_PHYS}/macroscopic_param_impl.h

	${PATH_SOLV}/solver.h
	${PATH_SOLV}/srt.h

	${PATH_BC}/bc.h	
	)

set(SOURCE_LIST
	${PATH_AREA}/fluid.cpp
	${PATH_AREA}/medium.cpp

	${PATH_SOLV}/srt.cpp

	${PATH_BC}/bc.cpp

	main.cpp
	)


add_executable(lbm ${HEADER_LIST} ${SOURCE_LIST})
