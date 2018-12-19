add_definitions(-DBOOST_ALL_NO_LIB)
find_package(PNG QUIET REQUIRED)
find_package(Boost QUIET REQUIRED COMPONENTS system filesystem)

add_subdirectory(deps/alglib)
add_subdirectory(deps/jsoncpp)

include_directories(src include deps/alglib/cpp/src deps/jsoncpp/include ${PNG_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

file(GLOB ARCSIM_SOURCES "src/*")
add_library(arcsim ${ARCSIM_SOURCES})
target_link_libraries(arcsim alglib json ${PNG_LIBRARIES} ${Boost_LIBRARIES})
