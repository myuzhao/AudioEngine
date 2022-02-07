cmake_minimum_required(VERSION 3.0)
project(AudioEngineTest)

include_directories(${PROJECT_SOURCE_DIR}/Include/)

file(GLOB_RECURSE SRC_FILES
	${PROJECT_SOURCE_DIR}/Src/*.cpp
)

file(GLOB_RECURSE MAIN_SRC_FILES
	${PROJECT_SOURCE_DIR}/Test/main/test_main.c
)

file(GLOB_RECURSE INCLUDE_FILES
	${PROJECT_SOURCE_DIR}/Include/*.h
)

add_executable(${PROJECT_NAME} ${MAIN_SRC_FILES} ${SRC_FILES} ${INCLUDE_FILES})

