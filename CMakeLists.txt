#
# CMAKE file for TestRunner
#
set(CMAKE_OSX_ARCHITECTURES x86_64)
set(CMAKE_OSX_SYSROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk)
cmake_minimum_required(VERSION 2.8)


project(testrunner)

set(CMAKE_SOURCE_DIR ./src)
set(CMAKE_BINARY_DIR ./bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ./lib)


set(CMAKE_BUILD_TYPE Release)

include_directories("${PROJECT_SOURCE_DIR}")

#
# Example test libraries
#
list(APPEND meshopt src/main.cpp)
list(APPEND meshopt src/vec.c)
list(APPEND meshopt src/vecutil.c)

# set_property(TARGET trunner PROPERTY CXX_STANDARD 11)


# osx testing
add_executable(moo ${meshopt})
target_include_directories(moo PRIVATE ./src)
target_include_directories(moo PRIVATE .)
target_include_directories(moo PRIVATE /usr/local/include)


set_property(TARGET moo PROPERTY CXX_STANDARD 11)
# set_target_properties(trun PROPERTIES PUBLIC_HEADER src/testrunner/testinterface.h)



set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib /usr/lib/system)
find_library(COCOA_FRAMEWORK Cocoa)
find_library(IOKIT_FRAMEWORK IOKit)
find_library(CORE_FRAMEWORK CoreFoundation)
find_library(CORE_VIDEO_FRAMEWORK CoreVideo)
find_library(SYSTEM_FRAMEWORK System)
find_library(OPENGL_FRAMEWORK OpenGL)
find_library(GLFW NAMES libglfw3.a PATHS)


message(STATUS ${COCOA_FRAMEWORK})
message(STATUS ${CMAKE_OSX_ARCHITECTURES})
message(STATUS ${GLFW})

target_link_libraries(moo ${COCOA_FRAMEWORK} ${IOKIT_FRAMEWORK} ${CORE_VIDEO_FRAMEWORK} ${OPENGL_FRAMEWORK} ${CORE_FRAMEWORK} ${GLFW})

# install(TARGETS trun DESTINATION bin)
# install(FILES src/testrunner/testinterface.h DESTINATION include)