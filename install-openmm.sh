mkdir .build

## Download each binary dependency from the Internet.

cd .build

curl -L -o "openmm-macos.zip" "https://github.com/philipturner/molecular-renderer-simulator-binaries/releases/download/v1.0.2/openmm-macos.zip"
unzip -o openmm-macos.zip

cd ../ # balance 'cd .build'

## Establish the vendors directory.

mkdir .build/vendors
cp ".build/openmm-macos/vendors/apple.icd" ".build/vendors/apple.icd"

## Copy binaries into the package folder.

cp ".build/openmm-macos/libc++.1.dylib" libc++.1.dylib
cp ".build/openmm-macos/libocl_icd_wrapper_apple.dylib" libocl_icd_wrapper_apple.dylib
cp ".build/openmm-macos/libOpenCL.1.dylib" libOpenCL.1.dylib
cp ".build/openmm-macos/libOpenMM.dylib" libOpenMM.dylib
cp ".build/openmm-macos/libOpenMMOpenCL.dylib" libOpenMMOpenCL.dylib

## Fix the self-link of each library.

# The "-id" command is required for "libOpenMM.dylib". The Workspace executable
# records this name when it links OpenMM. A bare file name is resolved relative
# to the current directory, so 'run.sh' must start from the package folder.

install_name_tool -id "libc++.1.dylib" libc++.1.dylib
install_name_tool -id "libocl_icd_wrapper_apple.dylib" libocl_icd_wrapper_apple.dylib
install_name_tool -id "libOpenCL.1.dylib" libOpenCL.1.dylib
install_name_tool -id "libOpenMM.dylib" libOpenMM.dylib
install_name_tool -id "libOpenMMOpenCL.dylib" libOpenMMOpenCL.dylib

## Fix the search path of the OpenMM plugin.

install_name_tool -rpath "@loader_path/../" "@loader_path/" libOpenMMOpenCL.dylib

## Repair the code signature.

codesign -fs - libc++.1.dylib
codesign -fs - libocl_icd_wrapper_apple.dylib
codesign -fs - libOpenCL.1.dylib
codesign -fs - libOpenMM.dylib
codesign -fs - libOpenMMOpenCL.dylib

echo "These code-signs should report success:"
codesign --verify --verbose libc++.1.dylib
codesign --verify --verbose libocl_icd_wrapper_apple.dylib
codesign --verify --verbose libOpenCL.1.dylib
codesign --verify --verbose libOpenMM.dylib
codesign --verify --verbose libOpenMMOpenCL.dylib
