cmake_minimum_required(VERSION 2.8.1)
project(pcapack)
enable_language(C)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wfatal-errors -O3 -fPIC -march=native")


### External dependencies
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

#find_package(Check)
#find_package(Cobra)

find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



### Source files
set(RAND
	src/rand/rand_normal.c
	src/rand/rand_tab.c
	src/rand/rand_uniform.c
	src/rand/rng_interface.c
	src/rand/samplers.c
	src/rand/utils.c
	src/rand/mrg/rand_mrg.c
	src/rand/mt/rand_mt.c
)

set(SUMSTATS
	src/sumstats/covariance.c
	src/sumstats/crossprod.c
	src/sumstats/means.c
	src/sumstats/scale.c
	src/sumstats/sweeps.c
	src/sumstats/variances.c
)

set(UTILS
	src/utils/distances.c
	src/utils/dnrm2.c
	src/utils/inverse.c
	src/utils/rank.c
	src/utils/xpose.c
)

SET(PCA
	src/fastmap.c
	src/svd.c
	src/pca.c
)

SET(DEBUG
	src/utils/debug/matprinter.c
)


set(PCAPACK_FILES ${RAND} ${SUMSTATS} ${UTILS} ${PCA} ${DEBUG})



### Make and install
add_library(pcapack_objs OBJECT ${PCAPACK_FILES})

add_library(pcapack SHARED $<TARGET_OBJECTS:pcapack_objs>)
add_library(pcapack_static STATIC $<TARGET_OBJECTS:pcapack_objs>)


target_link_libraries (pcapack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_target_properties(pcapack PROPERTIES PREFIX "")
if(NOT WIN32)
  set_target_properties(pcapack PROPERTIES SUFFIX ".so")
  install(TARGETS pcapack LIBRARY DESTINATION ${CMAKE_SOURCE_DIR})
else(NOT WIN32)
  install(TARGETS pcapack RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
endif(NOT WIN32)



#	if(CHECK_FOUND)
#		add_subdirectory(test)
#	endif()

##	if(COBRA_FOUND)
##		add_subdirectory(bench)
##	endif()
#endif()

#INSTALL(FILES build/libglm.so DESTINATION lib)
