**PROJECT: Apreture robust optical flow for event driven cameras**

Dependencies: Eigen3, gcc, g++
Install dependencies:
sudo apt-get install libeigen3-dev gcc g++ cmake

To install project (recommended):
mkdir build
cd build
cmake ..
make

To run the project:
cd build
./FARMSOut

**Currently the parameters need to be changed from within the main.cpp or vFlow.cpp source, this will be update in future updates**



Input event files: txt files with N rows for each event storing events as: x y p t.

Output files: with suffix _FARMSOut_ with N rows for each event storing flow events as: x y p t globalR globalTheta localR localTheta


Send Queries to : himanshu.akolkar@gmail.com (or use the bitbucket ticketing)