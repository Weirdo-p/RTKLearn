add_executable(navi test.cpp)

find_package(OpenMP)

target_link_libraries(navi common pnt optimal ambiguity OpenMP::OpenMP_CXX)

