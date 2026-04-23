@echo off
REM Build OASIS GUI as standalone Windows executable
REM Requires: pip install pyinstaller

pyinstaller --onedir --windowed --name OASIS_GUI ^
  --add-data "resources;resources" ^
  --hidden-import pyvistaqt ^
  --hidden-import pyvista ^
  --hidden-import vtkmodules ^
  --hidden-import vtkmodules.all ^
  --collect-all pyvista ^
  --collect-all pyvistaqt ^
  --collect-all vtkmodules ^
  --hidden-import ruamel.yaml ^
  --hidden-import ruamel.yaml.comments ^
  --hidden-import matplotlib.backends.backend_qt5agg ^
  main.py

echo.
echo Build complete. Find executable in dist\OASIS_GUI\
