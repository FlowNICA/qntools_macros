cmake_minimum_required(VERSION 3.12)
project(QnToolsInterface CXX)
set(CMAKE_CXX_STANDARD 17)
set(PROJECT_VERSION 1.0)

find_package(ROOT REQUIRED)
find_package(QnTools REQUIRED COMPONENTS QnTools::Base QnTools::Correction)

include_directories(
        ${CMAKE_SOURCE_DIR}
        ${ROOT_INCLUDE_DIRS}
        ${QnTools_INCLUDE_DIR}/QnTools
)

message(STATUS ${ROOT_INCLUDE_DIRS})
message(STATUS ${QnTools_INCLUDE_DIR})

set(SOURCES
        src/variable_manager.cc
        src/correction_task.cc
        src/vector_config.cc
        src/t_file_ptr.cc
        )
set(HEADERS
        src/variable_manager.h
        src/correction_task.h
        src/vector_config.h
        src/t_file_ptr.h
        )

add_library(Correction SHARED ${SOURCES} G__Correction.cxx )
ROOT_GENERATE_DICTIONARY( G__Correction ${HEADERS} LINKDEF src/CorrectionLinkDef.h )
target_link_libraries(Correction PUBLIC QnTools::Correction QnTools::Base ROOT::RIO ROOT::Tree ROOT::ROOTDataFrame )

add_executable( correct src/correct.cc )
target_link_libraries( correct PUBLIC Correction ROOT::RIO )

set( QNTOOLS_INCLUDE_PATH ${QnTools_INCLUDE_DIR}/QnTools )
configure_file(src/correlate.cc.in ${CMAKE_SOURCE_DIR}/src/correlate.cc @ONLY)

add_executable( correlate src/correlate.cc )
set_target_properties( correlate PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries( correlate PRIVATE ROOT::RIO ROOT::Geom QnTools::DataFrame QnTools::Base )

