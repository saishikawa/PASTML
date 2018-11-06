cmake_minimum_required(VERSION 3.6)
project(pastml)

set(CMAKE_C_STANDARD 11)
set(CMAKE_INCLUDE_PATH .)
set(CMAKE_C_FLAGS "--std=c11")

set(SOURCE_FILES main.c likelihood.c tree.c marginal_approximation.c marginal_likelihood.c
        states.c runpastml.c parsimony.c likelihood.h marginal_likelihood.h tree.h
        marginal_approximation.h states.h pastml.h runpastml.h param_minimization.c
        param_minimization.h scaling.c scaling.h logger.c logger.h parsimony.h metric.h metric.c)
add_executable(pastml ${SOURCE_FILES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(pastml GSL::gsl -lm -lgsl -lgslcblas)
