#!/bin/bash

# Whether to build the tests in release mode.
build_release=false

# Check whether the arguments parse correctly.
if [[ $# -gt 1 ]]; then
  echo "Too many arguments."
  invalid_input=true
elif [[ $# != 0 ]]; then
  invalid_input=false
 
  if [[ $invalid_input == false ]]; then
    for param in "$@"; do
      if [[ $param == "--release" ]]; then
        if [[ $build_release == true ]]; then
          echo "Duplicate argument '--release'."
          invalid_input=true
        else
          build_release=true
        fi
      else
        echo "Unrecognized argument '${param}'."
        invalid_input=true
      fi
    done
  fi
else
  invalid_input=false
fi

# Return early if the arguments are incorrect.
if [[ $invalid_input == true ]]; then
  echo "Usage: test.sh [--release]"
  exit -1
fi

export OPENMM_LIBRARY_PATH="$(pwd)"

# Compile the executable.
if [[ $build_release == true ]]; then
  swift test -Xswiftc -Ounchecked -Xswiftc -DRELEASE
  export XCTEST_FILE="$(pwd)/.build/debug/MM4PackageTests.xctest"
else
  swift test
  export XCTEST_FILE="$(pwd)/.build/debug/MM4PackageTests.xctest"
fi

# Workaround for the tests not running on macOS.
if [[ "$OSTYPE" == "darwin"* ]]; then
  export XCTEST_EXEC="$XCTEST_FILE/Contents/MacOS/MM4PackageTests"
  
  # Link libc++ and every plausible OpenMM version for the next 5 years.
  install_name_tool -change @rpath/libOpenMM.dylib "$OPENMM_LIBRARY_PATH/libOpenMM.dylib" $XCTEST_EXEC
  
  # Actually run the Swift package tests.
  swift test --skip-build
fi
