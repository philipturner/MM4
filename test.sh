# Activate the library linking mechanism in swift-openmm.
export OPENMM_LIBRARY_PATH="/C:/Users/phili/miniconda3/Library/lib"

# Run in release mode with incremental compilation.
swift test -Xswiftc -Ounchecked -Xswiftc -DRELEASE
