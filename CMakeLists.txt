cmake_minimum_required(VERSION 3.16)
project(ComputerGraph)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# this two line don't work new.
set(INCLUDE_FILE
    ${INCLUDE_DIR}/OBJLoader.hpp
    ${INCLUDE_DIR}/DrawableObject.hpp
    ${INCLUDE_DIR}/Transform.hpp
)
set(SRC_FILES
    ${SRC_DIR}/DrawableObject.cpp
    ${SRC_DIR}/Transform.cpp
)
add_subdirectory(thirdparty/freeglut)

# you need pre-install OpenGL, GLU, please google
add_executable(${PROJECT_NAME} ${SRC_FILES} ${SRC_DIR}/main.cpp)
target_link_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${OPENGL_glu_LIBRARY})
target_link_libraries(${PROJECT_NAME} freeglut_static)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/thirdparty/freeglut/include/GL ${INCLUDE_DIR})
target_compile_definitions(${PROJECT_NAME} PRIVATE RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/TestFile")