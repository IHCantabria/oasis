
rem DELETE PREVIOUS BUILD ---------------------------------
rmdir /Q /S build
rem BUILD COMPILATION RUN ---------------------------------
cmake -B./build/ -H./
rem CHOOSE COMPILATION OPTION -----------------------------
rem DEBUG COMPILATION ---------------------------------
rem cmake --build ./build
rem copy ".\build\Debug\OASIS.exe" ".\bin\oasis.exe"
rem RELEASE COMPILATION -------------------------------
cmake --build ./build --config Release
copy ".\build\Release\OASIS.exe" ".\bin\oasis.exe"