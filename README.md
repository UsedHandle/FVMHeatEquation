# FVMHeatEquation
![Image of output](https://github.com/UsedHandle/FVMHeatEquation/blob/main/output.PNG?raw=true)
The initial boundary condition has the center at a higher temperature than the surroundings. Using the FVM method on the equation and dividing the area into an orthogonal grid, the program simulates how heat diffuses in a substrate 
## build (CMake which fetches glfw3 and glm)
For mac or linux:
```
cmake .
make
```
On windows, ```make``` can be replaced with compilation through Visual Studio after setting FVMHeatEquation as the startup project
