cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(name    createTree)
#set(name    filterDS)
#set(name    filterSS)
#set(name    unfold_Double)
#set(name    unfold_Single)

set(mymain ${name}.cc )
set(myexe  ${name}    )
set(myproj ${name}_   )

project(${myproj})

# You need to tell CMake where to find the ROOT installation.
#list(APPEND CMAKE_PREFIX_PATH /opt/applications/root6/root-6.12.06-install/bin)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

include_directories(${PROJECT_SOURCE_DIR}/include)
        EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


#include(${ROOT_USE_FILE})
#include_directories(/opt/applications/root6/)

#---Create  a main program
add_executable(${myexe} ${mymain})
        EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})

target_link_libraries(${myexe} ${CMAKE_EXE_LINKER_FLAGS})

install(TARGETS ${myexe} DESTINATION .)
