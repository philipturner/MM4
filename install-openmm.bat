mkdir .build
cd .build

:: Download each binary dependency from the Internet.
powershell -c "Invoke-WebRequest -Uri 'https://github.com/philipturner/molecular-renderer-simulator-binaries/releases/download/v1.0.2/openmm-windows.zip' -OutFile 'openmm-windows.zip'"

:: Decompress each ZIP file into a hierarchy of folders.
powershell -c "Import-Module Microsoft.Powershell.Archive"
powershell -c "Expand-Archive -Force 'openmm-windows.zip' 'openmm-windows'"

copy /Y openmm-windows\openmm-windows\OpenMM.dll "../OpenMM.dll"
copy /Y openmm-windows\openmm-windows\OpenMM.lib "../OpenMM.lib"
copy /Y openmm-windows\openmm-windows\OpenMMOpenCL.dll "../OpenMMOpenCL.dll"
