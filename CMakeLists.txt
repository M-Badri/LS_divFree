
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
project (divergence_free_intrpolation LANGUAGES Fortran)

set(PROJECT_DIRECTORY "/home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2")
#add_subdirectory(third_party)

list(APPEND _sources src/m_interpolations.f90
		     		 src/m_least_square.f90
                     src/m_math.f90
                     src/mod_test.f90
                     src/mod_write.f90
                     src/m_tests.f90 )


add_executable (main src/main.f90 ${_sources})
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_DIRECTORY}/bin)

target_link_libraries(main /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/third_party/liblapack.a /home/mehdi/Projects/codes/afivo_run/div_free/Simplified_2/third_party/libblas.a )
