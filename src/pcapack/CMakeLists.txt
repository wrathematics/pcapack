cmake_minimum_required(VERSION 2.6)
project(pcapack)
enable_language(C Fortran)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_Fortran_MODULE_DIRECTORY "include/")

### Flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfatal-errors -O3 -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wfatal-errors -O3 -fPIC")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


### External dependencies
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

#find_package(Check)
#find_package(Cobra)

find_package(OpenMP)


if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  
  # OpenMP flags for Fortran --- thanks internet man who I stole this from!
  # http://sethmmorton.github.io/How%20To/2013/04/28/using-cmake-to-find-openmp-flags-for-fortran-compilers/
  find_package(OpenMP_Fortran)
  set_target_properties(${EXECUTABLE_NAME} PROPERTIES
                      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                      LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
  
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



### Source files
set(RAND
	src/rand/old/rand.c
	src/rand/old/rand_mt.c
	src/rand/old/rand_normal.c
)

set(SUMSTATS
	src/sumstats/covariance_utils.f90
	src/sumstats/covariance.f90
	src/sumstats/means.f90
	src/sumstats/sweeps.f90
	src/sumstats/variances.f90
)

set(UTILS
	src/utils/cbool.f90
	src/utils/distances.f90
	src/utils/float_equals.f90
	src/utils/lapack.f90
	src/utils/rank.f90
	src/utils/rank_.c
	src/utils/sgns.f90
	src/utils/stolen_blas.f
	src/utils/transposition.f90
)

SET(PCA
	src/svd.c
#	src/pca_mod.f90
#	src/svd_mod.f90
#	src/fastmap_mod.f90
)


set(PCAPACK_FILES ${RAND} ${SUMSTATS} ${UTILS} ${PCA})


### Make and install
add_library(pcapack SHARED ${PCAPACK_FILES})
add_library(pcapack_static STATIC ${PCAPACK_FILES})


target_link_libraries (pcapack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_target_properties(pcapack PROPERTIES PREFIX "")
if(NOT WIN32)
  set_target_properties(pcapack PROPERTIES SUFFIX ".so")
  install(TARGETS pcapack LIBRARY DESTINATION ${CMAKE_SOURCE_DIR})
else(NOT WIN32)
  install(TARGETS pcapack RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
endif(NOT WIN32)


#if (CHECK_FOUND OR COBRA_FOUND)
#	enable_testing()

#	if(CHECK_FOUND)
#		add_subdirectory(test)
#	endif()

##	if(COBRA_FOUND)
##		add_subdirectory(bench)
##	endif()
#endif()

#INSTALL(FILES build/libglm.so DESTINATION lib)
