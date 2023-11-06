#define NANOVDB_USE_OPENVDB
#include <nanovdb/NanoVDB.h>
#include <nanovdb/util/IO.h>
#include <nanovdb/util/OpenToNanoVDB.h>

//! No use, use nanovdb_converter
int main(int argc, char **argv) {
  openvdb::initialize();

  openvdb::io::File file(argv[1]);
  file.open();

  openvdb::GridBase::Ptr baseGrid;
  for (openvdb::io::File::NameIterator nameIter = file.beginName();
       nameIter != file.endName(); ++nameIter) {
    // TODO: deal with other channels if you have them
    if (nameIter.gridName() == "density") {
      baseGrid = file.readGrid(nameIter.gridName());
    }
  }
  file.close();

  auto handle = nanovdb::openToNanoVDB(baseGrid);

  // Write the NanoVDB grid to file
  nanovdb::io::writeGrid(argv[2], handle);
}