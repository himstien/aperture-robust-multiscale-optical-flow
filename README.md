**PROJECT: Apreture robust optical flow for event driven sensors**

This code performs aperture corrected flow computation over events obtained from DVS or ATIS sensors.

The events must be written in a txt file as:  N x 4 array where N = number of events 

Each row must be in form:

x y timestamp polarity


********************************************************************************************************************************************************************************************************************************
** INSTALLATION **

Dependencies: Eigen3, gcc, g++

**Install dependencies:**

sudo apt-get install libeigen3-dev gcc g++ cmake


**To install project (recommended):**

mkdir _build

cd _build

cmake ..

make


**To run the project:**

cd _build

Type ./FARMS_Flow --help for parameters list .. 


e.g:
./FARMS_Flow --width 320 --height 320 --filename ~/path/to/file/filename *note: filename is without extension

********************************************************************************************************************************
Parameter list:

  	--help                Displays this message

	--filename arg        add events file name without extension (.txt)

	--height arg          set sensor height

	--width arg           set sensor width

	--filtersize arg      set size of neighbor for plane fitting

	--inlierCheck arg     set minimum number of inliers to validate plane

	--numEvents arg       set max number of events to process


Input event files: txt files with N rows for each event storing events as: x y t p.

Output files: with suffix _FARMSOut_ with N rows for each event storing flow events as: x y t p globalR globalTheta localR localTheta

************************************************************************************************************************************
Please cite as:
H. Akolkar, S. H. Ieng and R. Benosman, "Real-time high speed motion prediction using fast aperture-robust event-driven visual flow," in IEEE Transactions on Pattern Analysis and Machine Intelligence, doi: 10.1109/TPAMI.2020.3010468.

************************************************************************************************************************************
Send Queries to : himanshu.akolkar@gmail.com (or use the bitbucket ticketing)
