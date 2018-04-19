cmake_minimum_required(VERSION 3.8)
project(pastml)

set(CMAKE_C_STANDARD 99)
set(CMAKE_INCLUDE_PATH .)

set(SOURCE_FILES main.c likelihood.c make_tree.c marginal_approximation.c marginal_likelihood.c
        output_states.c output_tree.c runpastml.c likelihood.h marginal_likelihood.h make_tree.h
        marginal_approximation.h output_tree.h output_states.h pastml.h runpastml.h param_minimization.c
        param_minimization.h scaling.c scaling.h logger.c logger.h output_simulation.h output_simulation.c
        joint_likelihood.h joint_likelihood.c models.c models.h eigen.c eigen.h)
add_executable(pastml ${SOURCE_FILES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(pastml GSL::gsl -lm)
