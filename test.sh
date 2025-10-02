# Activate the library linking mechanism in swift-openmm.
export OPENMM_LIBRARY_PATH="$(pwd)"

# Run in release mode with incremental compilation.
swift test -Xswiftc -Ounchecked -Xswiftc -DRELEASE
