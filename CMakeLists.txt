cmake_minimum_required(VERSION 3.16)
project(ComputerGraph)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(THIRDPARTY ${CMAKE_SOURCE_DIR}/thirdparty)
add_subdirectory(thirdparty/freeglut)

# this two line don't work new.
set(INCLUDE_FILES
    ${THIRDPARTY}/stbImage/stb_image.h
    ${INCLUDE_DIR}/OBJLoader.hpp
    ${INCLUDE_DIR}/DrawableObject.hpp
    ${INCLUDE_DIR}/Transform.hpp
    ${INCLUDE_DIR}/CGMath.hpp
    ${INCLUDE_DIR}/Texture.hpp
    ${INCLUDE_DIR}/Dp.hpp
    ${INCLUDE_DIR}/ImageLoader.hpp
)
set(SRC_FILES
    ${SRC_DIR}/DrawableObject.cpp
    ${SRC_DIR}/Transform.cpp
    ${SRC_DIR}/Texture.cpp
)

add_executable(${PROJECT_NAME}
    ${SRC_DIR}/main.cpp
    ${INCLUDE_FILES}
    ${SRC_FILES}
)

target_link_libraries(${PROJECT_NAME} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
target_link_libraries(${PROJECT_NAME} freeglut_static)
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE
    ${CMAKE_SOURCE_DIR}/thirdparty/freeglut/include/GL
    ${THIRDPARTY}/stbImage
)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)

target_compile_definitions(${PROJECT_NAME} PRIVATE RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Resource")
