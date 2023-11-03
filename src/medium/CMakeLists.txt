add_mts_plugin(homogeneous homogeneous.cpp)
add_mts_plugin(heterogeneous heterogeneous.cpp)

add_mts_plugin(nanovdbmedium nanovdb.cpp majorantgrid.cpp)
target_link_libraries(nanovdbmedium z tbb)