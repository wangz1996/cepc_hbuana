include_directories(include)
include_directories(src)

message("pid_path_cmake:" ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(PIDTool pid.cxx src/PIDTool.cxx src/HTTool.cxx)
target_link_libraries(PIDTool ${ROOT_LIBRARIES})
