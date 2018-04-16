/*
 * Copyright (C) 2015 iCub Facility - Istituto Italiano di Tecnologia
 * Author: arren.glover@iit.it
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * A copy of the license can be found at
 * http://www.robotcub.org/icub/license/gpl.txt
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

/*
 * @file main.cpp
 * @brief main code for the computation of the optical flow
 */

#include "../include/vFlow.h"
//#include <yarp/os/all.h>
//#include <yarp/sig/all.h>


int main(int argc, char * argv[])
{

    /* instantiate the module */
    int height = 320;
    int width = 320;
    int filterSize = 3;
    int minEvtsOnPlane = 5;

    int NUMEVENTS = 1000000;

    std::string fileNameInput = "/home/himanshu/POST_DOC/DATA/atisData/bar_square/multiPattern1_fixed_";

    vFlowManager vFlowM(height, width, filterSize,
                                     minEvtsOnPlane,  fileNameInput);

    std::cout << "[debug Main] : size of lastFlowTime is [sx sy]: [" << vFlowM.returnFlowTime().dim_a() << " " << vFlowM.returnFlowTime().dim_b() << "]" << std::endl;
    vFlowM.run(NUMEVENTS);

    return 0;
}
