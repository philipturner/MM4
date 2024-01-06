#!/bin/bash

# Whether the next argument is the user. This is primarily to make testing more
# ergonomic on macOS.
# Example: "philipturner" in "/Users/philipturner/miniforge3/lib"
next_argument_user=false

# Whether the user has been specified.
user_specified=false

# The name of the user.
user_name=""

# Whether the next argument is the path.
next_argument_path=false

# Whether the path has been specified.
path_specified=false

# The name of the path.
path_name=""

# Whether to build the tests in release mode.
build_release=false

# Check whether the arguments parse correctly.
if [[ $# -gt 5 ]]; then
  echo "Too many arguments."
  invalid_input=true
elif [[ $# != 0 ]]; then
  invalid_input=false
 
  if [[ $invalid_input == false ]]; then
    for param in "$@"; do
      if [[ $param == "--user" ]]; then
        if [[ $next_argument_user == true ]]; then
          echo "Duplicate argument '--user'."
          invalid_input=true
        elif [[ $user_specified == true ]]; then
          echo "Duplicate argument '--user'."
          invalid_input=true
        elif [[ $next_argument_path == true ]]; then
          echo "Invalid use of '--path'."
          invalid_input=true
        else
          next_argument_user=true
        fi
      elif [[ $param == "--path" ]]; then
        if [[ $next_argument_path == true ]]; then
          echo "Duplicate argument '--path'."
          invalid_input=true
        elif [[ $path_specified == true ]]; then
          echo "Duplicate argument '--path'."
          invalid_input=true
        elif [[ $next_argument_user == true ]]; then
          echo "Invalid use of '--user'."
          invalid_input=true
        else
          next_argument_path=true
        fi
      elif [[ $param == "--release" ]]; then
        if [[ $build_release == true ]]; then
          echo "Duplicate argument '--release'."
          invalid_input=true
        elif [[ $next_argument_user == true ]]; then
          echo "Invalid use of '--user'."
          invalid_input=true
        elif [[ $next_argument_path == true ]]; then
          echo "Invalid use of '--path'."
          invalid_input=true
        else
          build_release=true
        fi
      elif [[ $next_argument_user == true ]]; then
        next_argument_user=false
        user_specified=true
        user_name="${param}"
      elif [[ $next_argument_path == true ]]; then
        next_argument_path=false
        path_specified=true
        path_name="${param}"
      else
        echo "Unrecognized argument '${param}'."
        invalid_input=true
      fi
    done
  fi
else
  echo "No arguments found."
  invalid_input=true
fi

# If the arguments parse correctly, check whether the path is correct.
if [[ $invalid_input == false ]]; then
  if [[ $user_specified == true ]]; then
    if [[ $path_specified == true ]]; then
      echo "Choose either user or path."
      invalid_input=true
    else
      # Automatically detect the OpenMM library path. This is only tested on
      # arm64 macOS.
      export OPENMM_LIBRARY_PATH="/Users/${user_name}/miniforge3/lib"
    fi
  elif [[ $path_specified == true ]]; then
    export OPENMM_LIBRARY_PATH="${path_name}"
  else
    echo "No user or path specified."
    invalid_input=true
  fi
fi

# Return early if the arguments are incorrect.
if [[ $invalid_input == true ]]; then
  echo "Usage: test.sh [--user USER] [--path OPENMM_LIBRARY_PATH] [--release]"
  exit -1
fi

# Compile the executable. The tests shouldn't actually run, and no error should
# appear, because the command is broken.
if [[ $build_release == true ]]; then
#    swift test -c release
#    export XCTEST_FILE="$(pwd)/.build/release/MM4PackageTests.xctest"
  
  swift test -Xswiftc -Ounchecked
  export XCTEST_FILE="$(pwd)/.build/debug/MM4PackageTests.xctest"
else
  swift test
  export XCTEST_FILE="$(pwd)/.build/debug/MM4PackageTests.xctest"
fi

# Unsure what to do for other platforms. Only the exec path for macOS is known.
export XCTEST_EXEC="$XCTEST_FILE/Contents/MacOS/MM4PackageTests"

# Link libc++ and every plausible OpenMM version for the next 5 years.
install_name_tool -change @rpath/libOpenMM.8.0.dylib "$OPENMM_LIBRARY_PATH/libOpenMM.8.0.dylib" $XCTEST_EXEC
install_name_tool -change @rpath/libOpenMM.8.1.dylib "$OPENMM_LIBRARY_PATH/libOpenMM.8.1.dylib" $XCTEST_EXEC
install_name_tool -change @rpath/libOpenMM.8.2.dylib "$OPENMM_LIBRARY_PATH/libOpenMM.8.2.dylib" $XCTEST_EXEC
install_name_tool -change @rpath/libOpenMM.8.3.dylib "$OPENMM_LIBRARY_PATH/libOpenMM.8.3.dylib" $XCTEST_EXEC
install_name_tool -change @rpath/libOpenMM.8.4.dylib "$OPENMM_LIBRARY_PATH/libOpenMM.8.4.dylib" $XCTEST_EXEC
install_name_tool -change @rpath/libOpenMM.8.5.dylib "$OPENMM_LIBRARY_PATH/libOpenMM.8.5.dylib" $XCTEST_EXEC
install_name_tool -change @rpath/libc++.1.dylib "$OPENMM_LIBRARY_PATH/libc++.1.dylib" $XCTEST_EXEC

# Actually run the Swift package tests.
xcrun xctest $XCTEST_FILE
