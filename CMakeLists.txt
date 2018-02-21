cmake_minimum_required(VERSION 3.8)
project(pastml)

set(CMAKE_C_STANDARD 99)
set(CMAKE_INCLUDE_PATH .)

set(SOURCE_FILES main.c lik.c make_tree.c marginal_approxi.c marginal_lik.c nrutil.c
        output_states.c output_tree.c runpastml.c lik.h marginal_lik.h make_tree.h
        marginal_approxi.h output_tree.h output_states.h pastml.h runpastml.h param_minimization.c param_minimization.h)
add_executable(pastml ${SOURCE_FILES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(pastml GSL::gsl -lm)
