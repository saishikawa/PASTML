cmake_minimum_required(VERSION 3.8)
project(pastml)

set(CMAKE_C_STANDARD 99)
set(CMAKE_INCLUDE_PATH .)

set(SOURCE_FILES main.c fletcher.c fletcherJC.c golden.c lik.c make_tree.c marginal_approxi.c marginal_lik.c nrutil.c
        output_states.c output_tree.c runpastml.c lik.h marginal_lik.h make_tree.h golden.h fletcher.h fletcherJC.h
        marginal_approxi.h output_tree.h output_states.h asrml.h runpastml.h)
add_executable(pastml ${SOURCE_FILES})
target_link_libraries(pastml -lm)