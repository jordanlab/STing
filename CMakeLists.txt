# # Minimum cmake version
# cmake_minimum_required(VERSION 3.0)
# # Name of project and that it is C++ only.
# project(cmake_sting CXX)

# set(CMAKE_CXX_STANDARD 11)

# # Custom location of libraries and headers (SeqAn )
# set(SEQAN_INCLUDE_DIR "include")
# set(SEQAN_CMAKE_DIR "share/cmake")
# set(CMAKE_MODULE_PATH "share/cmake/Modules")

# set(CMAKE_PREFIX_PATH ${SEQAN_CMAKE_DIR})

# set(SEQAN_INCLUDE_PATH ${SEQAN_INCLUDE_DIR})

# # Find includes in corresponding build directories (for ui_xxxxxx.h)
# set(CMAKE_INCLUDE_CURRENT_DIR ON)
# # Instruct CMake to run moc automatically when needed.
# # set(CMAKE_AUTOMOC ON)

# # ----------------------------------------------------------------------------
# # Dependencies
# # ----------------------------------------------------------------------------
# find_package(Threads)
# # Search for zlib as a dependency for SeqAn.
# find_package(ZLIB)

# # set(STING_CXX_FLAGS "-std=c++11 -g -O3 ")
# # set(STING_CXX_FLAGS "-W -Wall -O3 -pedantic -std=c++11")
# set(STING_CXX_FLAGS "-W -Wall -O3 -pedantic")
# set(STING_FIRST_CXX_FLAGS "-std=c++11")
# # Add CXX flags found by find_package (SeqAn).
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STING_FIRST_CXX_FLAGS}")
# # set(STING_CXX_DEFINES "-DSEQAN_HAS_EXECINFO=1 -DSEQAN_HAS_ZLIB=1 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE")

# # Load the SeqAn module and fail if not found.
# find_package(SeqAn REQUIRED)
# # Add include directories.
# include_directories(
#   ${SEQAN_INCLUDE_DIRS}
# )
# # Add definitions set by find_package (SeqAn).
# add_definitions(
#   ${SEQAN_DEFINITIONS}
# )
# # Add CXX flags found by find_package (SeqAn).
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${STING_CXX_FLAGS}")
# # ----------------------------------------------------------------------------
# # Build Setup
# # ----------------------------------------------------------------------------
# # set(INDEXER_SOURCES src/Indexer.h src/Indexer.cpp src/IndexerOptions.h src/IndexerOptions.cpp src/types.hpp src/custom_utils.hpp src/data_structures.hpp)
# # set(DETECTOR_SOURCES src/Detector.h src/Detector.cpp src/KmerMatcher.h src/KmerMatcher.cpp src/KmerMatcherOptions.h src/KmerMatcherOptions.cpp src/types.hpp src/custom_utils.hpp src/data_structures.hpp)

# # SeqAn typer application
# set(typer "typer")
# set(${typer}_HEADERS
#   src/Typer.h
#   src/KmerMatcher.h
#   src/KmerMatcherOptions.h 
#   src/SequenceInfo.h
#   src/Timer.h
# )
# set(${typer}_SOURCES
#   src/Typer.cpp
#   src/KmerMatcher.cpp
#   src/KmerMatcherOptions.cpp 
#   src/SequenceInfo.cpp
#   src/Timer.cpp
#   src/types.hpp 
#   src/custom_utils.hpp 
#   src/data_structures.hpp
# )
# # message("SEQAN_INCLUDE_DIRS: ${SEQAN_INCLUDE_DIRS}")
# # get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# # foreach(dir ${dirs})
# #   message(STATUS "dir='${dir}'")
# # endforeach()
# # qt5_wrap_cpp(${seqanApp}_SOURCES_MOC ${${seqanApp}_HEADERS})

# # Add executable and link against SeqAn dependencies.
# # add_executable(seqanApp ${seqanApp_SOURCES})
# add_executable(${typer} 
#   ${${typer}_HEADERS}
#   ${${typer}_SOURCES}
#   # ${${seqanApp}_SOURCES_MOC}
# )

# target_link_libraries(${typer}
#     ${SEQAN_LIBRARIES}
#     ${CMAKE_THREAD_LIBS_INIT}
# )

# message("CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS}) 
# # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STING_CXX_FLAGS}")
# # include(add_ide_support.cmake)
# # target_add_ide_support(${qtGuiApp})

# =================================================================

project(sting_cmake)
cmake_minimum_required (VERSION 3.0.0)
#aux_source_directory(. SRC_LIST)

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
#find_package (RT)

## Load the SeqAn module and fail if not found.
## find_package (SeqAn REQUIRED)

include_directories(include)

set(STING_CXX_FLAGS "-std=c++11 -W -Wall -pedantic -O3 -DNDEBUG")
# set(STING_LD_FLAGS "-D_GLIBCXX_USE_C99=1 -lpthread -lz -lrt")
set(STING_LD_FLAGS "-static-libstdc++ -static-libgcc -D_GLIBCXX_USE_C99=1 -lpthread -lz -lrt")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STING_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STING_LD_FLAGS}")

#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")

set(TYPER_SOURCES
  src/Typer.h 
  src/Typer.cpp
  src/KmerMatcher.h
  src/KmerMatcher.cpp
  src/Timer.h
  src/Timer.cpp
  src/SequenceInfo.h
  src/SequenceInfo.cpp
  src/KmerMatcherOptions.h
  src/KmerMatcherOptions.cpp
  src/types.hpp
  src/custom_utils.hpp
  src/data_structures.hpp)

set(INDEXER_SOURCES
  src/Indexer.h
  src/Indexer.cpp
  src/IndexerOptions.h
  src/IndexerOptions.cpp
  src/types.hpp
  src/custom_utils.hpp
  src/data_structures.hpp)

set(DETECTOR_SOURCES
  src/Detector.h
  src/Detector.cpp
  src/KmerMatcher.h
  src/KmerMatcher.cpp
  src/KmerMatcherOptions.h
  src/KmerMatcherOptions.cpp
  src/types.hpp
  src/custom_utils.hpp
  src/data_structures.hpp)

#add_executable(hello src/main.cpp)

add_executable(indexer ${INDEXER_SOURCES})
target_link_libraries (indexer ${STING_LD_FLAGS})

add_executable(typer ${TYPER_SOURCES})
target_link_libraries (typer ${STING_LD_FLAGS})

add_executable(detector ${DETECTOR_SOURCES})
target_link_libraries (detector ${STING_LD_FLAGS})