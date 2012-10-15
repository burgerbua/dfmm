######################################################################
### add here your include directories
include_directories (
	${BASE_PATH}
	${BOOST_INCLUDE_PATH}
)


# create test programs
######################################################################
### add here your executables

# performs hadamard products of vector
add_library (had_lib
  ${BASE_PATH}/sources/zhad.f
  ${BASE_PATH}/sources/dhad.f)





# points
add_executable (points ./main_points.C)

# symmetries
add_executable (symmetries ./main_symmetries.C)

# octree
add_executable (octree ./main_octree.C)

# dfmm
add_executable (dfmm ./main_dfmm.C)
target_link_libraries (dfmm
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	had_lib
	)

# bbfmm
add_executable (bbfmm ./main_bbfmm.C)
target_link_libraries (bbfmm
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	had_lib
	)

# aca
add_executable (aca ./main_aca.C)
target_link_libraries (aca
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	)

# aca316
add_executable (aca316 ./main_aca316.C)
target_link_libraries (aca316
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	)

# interpolator
add_executable (interpolator ./main_interpolator.C)
target_link_libraries (interpolator
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
	)




######################################################################
### set here output path for your executables

set_target_properties (
  points
	symmetries
	octree			
	dfmm
	bbfmm
	aca
	aca316
	interpolator
	PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_PATH})

#set_target_properties (dfmm
#	PROPERTIES OUTPUT_NAME
#	dfmmL${DFMM_INTERPOLATIONORDER}E${DFMM_LOWRANKEPSILON})
#set_target_properties (bbfmm
#	PROPERTIES OUTPUT_NAME
#	bbfmmL${DFMM_INTERPOLATIONORDER}E${DFMM_LOWRANKEPSILON})


#set_target_properties (dfmm
#	PROPERTIES OUTPUT_NAME dfmm_${DFMM_ACCURACY})
#set_target_properties (bbfmm
#	PROPERTIES OUTPUT_NAME bbfmm_${DFMM_ACCURACY})


######################################################################
### testing
add_test(
	NAME points
	COMMAND points
	)

add_test(
	NAME octree
	WORKING_DIRECTORY ${BASE_PATH}/cmk
	COMMAND octree ./grid/e2_64_119698.inp 119698 7
	)

add_test(
	NAME interpolator
	COMMAND interpolator
	)

add_test(
	NAME aca
	COMMAND aca .3
	)

add_test(
	NAME aca316
	COMMAND aca316 .3
	)

add_test(
	NAME dfmm
	WORKING_DIRECTORY ${BASE_PATH}/cmk
	COMMAND dfmm ./grid/e2_16_7667.inp 7667 3 .3 1
	)

add_test(
	NAME bbfmm
	WORKING_DIRECTORY ${BASE_PATH}/cmk
	COMMAND bbfmm ./grid/e2_16_7667.inp 7667 3 1
	)
