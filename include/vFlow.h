/*
 * Author: himanshu.akolkar@gmail.com
 *
*/


#ifndef __VFLOW__
#define __VFLOW__

#include <string>
#include <math.h>
#include <Eigen/Core>
#include <vector>
#include <array>

#include "Event.h"
#include "FlowEvent.h"
#include "EventMatrix.h"

///
///
class vFlowManager
{
private:

  bool DEBUGMODE  = false;
	const double MAXSTAMP = pow(2, 32);
	const double TSTOSEC = 1e-6;

	int windowJump = 5;
	int maxWindow = 20;

	FlowEvent trueFlow;
  Event lastEvent;

  std::string fileNameInput;
  int height = 320; //! sensor height
  int width = 240;  //! sensor width
  int fRad = 3;   //! filter radius
  unsigned int planeSize; //! edge length of fitted plane
  int halfCount; //! plane area divided by 2
  int minEvtsOnPlane = 5; //! minimum number of events for plane validity

  EventMatrix<Event> surfaceOnL;
  EventMatrix<Event> surfaceOnR;
  EventMatrix<Event> surfaceOfL;
  EventMatrix<Event> surfaceOfR;

  EventMatrix<Event> cSurf;


	//extern double flowSurfaceLengthOn[][];     //!
	EventMatrix<double> flowSurfaceLengthOn; //(width, height);

	//extern double flowSurfaceLengthOf[width][height];// [width][height];     //!
	EventMatrix<double> flowSurfaceLengthOf; //(width, height);

	//double flowSurfaceThetaOn [width][height];     //!
	EventMatrix<double> flowSurfaceThetaOn; //(width, height);

	//extern double flowSurfaceThetaOf[][]; //[width][height]; //[height][width]; //!
	EventMatrix<double> flowSurfaceThetaOf; //(width, height);

	//double flowSurfaceVx [width][height]; 	//!
	EventMatrix<double> flowSurfaceVx; //(width, height);

	//double flowSurfaceVy [width][height];	//!
	EventMatrix<double> flowSurfaceVy; //(width, height);

	//double lastFlowTime [width][height];     //!
	EventMatrix<double> lastEventTime;// = new EventMatrix<double> (width, height, 0);
	//std::vector< double > lastFlowTime; //(width, height);

  Eigen::MatrixXd At;
  Eigen::MatrixXd AtA;
  Eigen::MatrixXd A2;
  std::vector< double > abc;


  int eventsComputed;
  int eventsPotential;
  int bottleCount;

  FlowEvent computeLocalFlow();
  int computeGrads(Eigen::MatrixXd A, Eigen::MatrixXd Y,
                    double cx, double cy, double cz,
                    double &dtdy, double &dtdx);
  int computeGrads(std::vector<Event> subsurf, Event &cen,
                    double &dtdy, double &dtdx);

	FlowEvent computeTrueFlow(int x, int y, unsigned int time, int pol);
	FlowEvent computeTrueFlow(int x, int y, int pol);

public:

    vFlowManager(int height, int width, int filterSize, int minEvtsOnPlane);
    vFlowManager(int height, int width, int filterSize, int minEvtsOnPlane, std::string fileName);

    bool run(unsigned long int); // run over given number of events only
    void close(); // close the program and clear memory. *not defined yet*

	EventMatrix<double> returnFlowTime(){return lastEventTime; }; //[DEBUG] test function to check if arrays are initialized correctly.

};

/******************************************************************************/
#endif //__VFLOW__
