cmake_minimum_required(VERSION 3.8)
project(scalerel)

add_subdirectory(./mylib)
add_executable(base base.cpp)
add_executable(main main.cpp)

target_link_libraries(main reliabilitylib)
target_link_libraries(base reliabilitylib ${CMAKE_SOURCE_DIR}/SAPPOROBDD/lib/BDD64.a)
target_compile_definitions(base PUBLIC B_64)