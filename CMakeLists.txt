cmake_minimum_required(VERSION 3.28)
project(ProblemaDeLosNeutrones C)

set(CMAKE_C_STANDARD 11)

add_executable(ProblemaDeLosNeutrones main.c
        cJSON.c)
