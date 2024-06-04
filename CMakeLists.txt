cmake_minimum_required(VERSION 3.16)
project(ComputerGraph)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(THIRDPARTY ${CMAKE_SOURCE_DIR}/thirdparty)

# this two line don't work new.
set(INCLUDE_FILES
    ${THIRDPARTY}/stbImage/stb_image.h
)
set(SRC_FILES
    ${THIRDPARTY}/gltools/src/GLBatch.cpp
    ${THIRDPARTY}/gltools/src/glew.c
    ${THIRDPARTY}/gltools/src/GLShaderManager.cpp
    ${THIRDPARTY}/gltools/src/GLTools.cpp
    ${THIRDPARTY}/gltools/src/GLTriangleBatch.cpp
    ${THIRDPARTY}/gltools/src/math3d.cpp
)
add_subdirectory(thirdparty/freeglut)

# you need pre-install OpenGL, GLU, please google
add_executable(${PROJECT_NAME} ${SRC_DIR}/shadow.cpp ${INCLUDE_FILES} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${OPENGL_glu_LIBRARY})
target_link_libraries(${PROJECT_NAME} freeglut_static)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/thirdparty/freeglut/include/GL ${THIRDPARTY}/stbImage ${THIRDPARTY}/gltools/include)
target_compile_definitions(${PROJECT_NAME} PRIVATE RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/texture")