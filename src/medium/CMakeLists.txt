add_mts_plugin(homogeneous homogeneous.cpp)
add_mts_plugin(heterogeneous heterogeneous.cpp)

add_mts_plugin(nanovdbmedium nanovdb.cpp)
target_link_libraries(nanovdbmedium z tbb)

#add_mts_executable(open2nano open2nano.cpp)
#target_link_libraries(open2nano openvdb)