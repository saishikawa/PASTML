cmake_minimum_required(VERSION 3.8)
project(pastml)

set(CMAKE_C_STANDARD 99)
set(CMAKE_INCLUDE_PATH .)

set(SOURCE_FILES main.c likelihood.c tree.c marginal_approximation.c marginal_likelihood.c
        states.c runpastml.c parsimony.c likelihood.h marginal_likelihood.h tree.h
        marginal_approximation.h states.h pastml.h runpastml.h param_minimization.c
        param_minimization.h scaling.c scaling.h logger.c logger.h parsimony.h models.h models.c eigen.h eigen.c)
add_executable(pastml ${SOURCE_FILES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(pastml GSL::gsl -lm -lgsl -lgslcblas)
