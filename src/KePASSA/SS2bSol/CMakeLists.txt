cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)


#>> Set up Project and Subdirectories
project (Propagator)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)



#>> Source Code
add_executable (Propagator
				src/Main.f90

				src/DECSOL/decsol.f
				src/DECSOL/dc_decsol.f
				
				src/Integrator/RADAU_Core.f
				src/Integrator/RADAU_Driver.f90

				src/SS2b.f90
				)



#>> Compilation Options
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")

	#target_compile_options (Propagator PUBLIC "-g")                 # DEBUG flag
	#target_compile_options (Propagator PUBLIC "-double-size 64")    # Set 'DOUBLE PRECISION' to 8 bytes long.
	target_compile_options (Propagator PUBLIC "-r16")	            # Default 'REAL' is 8 bytes long. Alternatively:  -real_size 64
	target_compile_options (Propagator PUBLIC "-O2") 	            # Optimization level -O2

elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")

	#target_compile_options (Propagator PUBLIC "-g")                 # DEBUG flag
	#target_compile_options (Propagator PUBLIC "-fdefault-double-8") # 'DOUBLE PRECISION' not promoted to 16 bytes, but stays 8 bytes instead
	target_compile_options (Propagator PUBLIC "-fdefault-real-8")	# Default 'REAL' is 8 bytes long
	target_compile_options (Propagator PUBLIC "-O2") 				# Optimization level -O2

endif ()

