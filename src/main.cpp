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
#include <boost/program_options.hpp>
#include <exception>

namespace options = boost::program_options;

int main(int argc, char * argv[])
{

  /* instantiate default values for parameters */
  int height = 320;
  int width = 320;
  int filterSize = 3;
  int minEvtsOnPlane = 5;


  unsigned long int NUMEVENTS = pow(2, 63);

  std::string fileNameInput = "/home/himanshu/POST_DOC/DATA/atisData/bar_square/multiPattern1_fixed_";

  try
  {
    options::options_description desc("Allowed options");
    desc.add_options()
          ("help", "Displays this message")
          ("filename", options::value<std::string>(), "add events file name without extension (.txt)")
          ("height", options::value<int>(), "set sensor height")
          ("width", options::value<int>(), "set sensor width")
          ("filtersize", options::value<int>(), "set size of neighbor for plane fitting")
          ("inlierCheck", options::value<int>(), "set minimum number of inliers to validate plane")
          ("numEvents", options::value<int>(), "set max number of events to process");


    options::variables_map vm;
    options::store(options::parse_command_line(argc, argv, desc), vm);
    options::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 0;
    }

//set filename
    if (vm.count("filename"))
    {
        std::cout << "filename set to "
             << vm["filename"].as<std::string>() << ".\n";
        fileNameInput = vm["filename"].as<std::string>();
    }
    else
    {
//
    }

//set height
    if (vm.count("height"))
    {
        std::cout << "height set to "
             << vm["height"].as<int>() << ".\n";
        height = vm["height"].as<int>();

    }
    else
    {
    }

//set width
    if (vm.count("width"))
    {
        std::cout << "width set to "
             << vm["width"].as<int>() << ".\n";
        width = vm["width"].as<int>();

    }
    else
    {
    }

//set filtersize
    if (vm.count("filtersize"))
    {
        std::cout << "filtersize set to "
             << vm["filtersize"].as<int>() << ".\n";
        filterSize = vm["filtersize"].as<int>();

    }
    else
    {
    }

//set inlierCheck
    if (vm.count("inlierCheck"))
    {
        std::cout << "inlierCheck set to "
             << vm["inlierCheck"].as<int>() << ".\n";
        minEvtsOnPlane = vm["inlierCheck"].as<int>();

    }
    else
    {
    }

//set numEvents
    if (vm.count("numEvents"))
    {
        std::cout << "numEvents set to "
             << vm["numEvents"].as<int>() << ".\n";
        NUMEVENTS = vm["numEvents"].as<int>();

    }
    else
    {
    }

  }
  catch(std::exception& e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...)
  {
      std::cerr << "Exception of unknown type!\n";
  }



    vFlowManager vFlowM(height, width, filterSize,
                                     minEvtsOnPlane,  fileNameInput);

    std::cout << "[debug Main] : size of lastFlowTime is [sx sy]: [" << vFlowM.returnFlowTime().dim_a() << " " << vFlowM.returnFlowTime().dim_b() << "]" << std::endl;
    vFlowM.run(NUMEVENTS);

    return 0;
}
