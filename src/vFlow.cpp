/*
 * Author: himanshu.akolkar@gmail.com
*/

#include "../include/vFlow.h"

//#include "../include/FlowEvent.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <Eigen/LU>

#define DEBUG
/******************************************************************************/

/******************************************************************************/
//vFlowManager
/******************************************************************************/
vFlowManager::vFlowManager(int height, int width, int filterSize,
                                     int minEvtsOnPlane, std::string fileNameInput)
{

	std::cout << "[debug] : Begin creating vFlowManager object" << std::endl;
    this->height=height;
    this->width=width;

	this->fileNameInput = fileNameInput;
    //ensure sobel size is at least 3 and an odd number
    if(filterSize < 5) filterSize = 3; // 3
    if(!(filterSize % 2)) filterSize--;
    this->fRad = filterSize / 2; // 1.5
    this->halfCount = pow(filterSize, 2.0) / 2.0 + 0.5;
    this->planeSize = pow(filterSize, 2.0); // 9

    this->minEvtsOnPlane = minEvtsOnPlane;

    //for speed we predefine the mememory for some matricies
    At = Eigen::MatrixXd(3, filterSize * filterSize);
    AtA = Eigen::MatrixXd(3,3);
    abc = std::vector<double>(3);
    A2 = Eigen::MatrixXd(3, 3);

	//flowSurfaceThetaOf = new double[width][height];
	flowSurfaceThetaOf = EventMatrix<double> (width, height);

	//double flowSurfaceThetaOn[width][height] = 0;
	flowSurfaceThetaOn = EventMatrix<double> (width, height);

	//double flowSurfaceLengthOf[width][height] = 0;
	flowSurfaceLengthOf = EventMatrix<double> (width, height);

	//double flowSurfaceLengthOn[width][height] = 0;
	flowSurfaceLengthOn = EventMatrix<double> (width, height);

	//double flowSurfaceVx[width][height] = 0;;
	flowSurfaceVx = EventMatrix<double> (width, height);

	//flowSurfaceVy = new double[width] [height];
	flowSurfaceVy = EventMatrix<double> (width, height, 0);


	//lastEventTime = Eigen::MatrixXd(width, height);
	lastEventTime = EventMatrix<double>  (width, height, 0.0);
	lastEventTime.resize(width, height, 0.0);


	windowJump = 10; //width/40; // 10
	maxWindow = 100; //width/2; //100

  eventsComputed = 0;
  eventsPotential = 0;
  bottleCount = 0;

	Event *tmp = new Event(0, 0, 0, 0);
    //create our surface in synchronous mode
    surfaceOnL = EventMatrix<Event> (width, height, *tmp);
    //surfaceOnL = new Eigen::MatrixXd(width, height);

    surfaceOfL = EventMatrix<Event>  (width, height, *tmp);
    //surfaceOfL = new Eigen::MatrixXd(width, height);

    surfaceOnR = EventMatrix<Event>  (width, height, *tmp);
    //surfaceOnR = new Eigen::MatrixXd(width, height);

    surfaceOfR = EventMatrix<Event> (width, height, *tmp);

    cSurf = EventMatrix<Event> (width, height, *tmp);
    //surfaceOfR = new Eigen::MatrixXd(width, height);

    //this->run();

	DEBUGMODE = true;
	
	if(DEBUGMODE)
	{
 		std::cout << "[debug in Run] : size of lastEventTime is [sx sy]: [" << lastEventTime.dim_a() << " " << lastEventTime.dim_b() << "]" << std::endl;
    	std::cout << "[debug] : End creating vFlowManager object" << std::endl;
	}



}

/******************************************************************************/
bool vFlowManager::run(unsigned long int NUMEVENTS)
{

	if(DEBUGMODE)
	{
		std::cout << "[debug] : Begin run" << std::endl;

		std::cout << "[debug in Run] : Setting up initial values " << std::endl;
	}

	std::string outFileName = "";
	outFileName = fileNameInput + "_FARMSOut_bench.txt";

	std::ofstream eventsFileOut;
	eventsFileOut.open (outFileName.c_str());

	double i = 0;
	int x; int y; unsigned int time_; int pol;

	fileNameInput = fileNameInput + ".txt";
	std::string line;

	std::cout << fileNameInput << std::endl;

	std::ifstream eventsFile (fileNameInput.c_str());

  	eventsFile.seekg(0, std::ios::beg);
	std::streampos fsize = eventsFile.tellg();
	eventsFile.seekg( 0, std::ios::end );
	fsize = eventsFile.tellg() - fsize;


  	NUMEVENTS = (unsigned long int) (std::min(double(NUMEVENTS), double(fsize/18)));

	eventsFile.seekg(0, std::ios::beg);

	if(DEBUGMODE)
	{
		std::cout << "[debug in Run] : Done setting initial values " << std::endl;

		std::cout << "[debug in Run] : Max window value " << maxWindow << std::endl;


	}

	Event *currentEvent = new Event(0, 0, 0, 0);
	cSurf = EventMatrix<Event> (width, height, *currentEvent);

	  if (eventsFile.is_open())
	  {
		getline (eventsFile,line);
		std::stringstream stream(line);

		stream >> x;
		stream >> y;
		stream >> time_;
		stream >> pol;

		currentEvent->setX(x);
		currentEvent->setY(y);
		currentEvent->setStamp(time_);
		currentEvent->setPolarity(pol);

		if(DEBUGMODE)
		{
			std::cout << "[debug in Run] : First event " << std::endl;
		}
		unsigned int t0 = time_;

		std::cout << "First time = " << t0 << std::endl;

		if(DEBUGMODE)
		{
			std::cout << "[debug in Run] : First event - initial event [x y p t]: [" << x << " " << y << " " << pol << " " << time_ << "]" << std::endl;
		}
		lastEventTime[x][y] = time_;

		if(DEBUGMODE)
		{
			std::cout << "[debug in Run] : While events loop - Set up lastEventTime[x][y] " << std::endl;
		}

	    while (getline (eventsFile,line) && eventsComputed <= NUMEVENTS)
	    {

			std::stringstream stream(line);

			stream >> x;
			stream >> y;
			stream >> time_;
			time_ = time_ - t0;

			stream >> pol;
			if(pol < 0)
				pol = 0;
			//pol = pol+1;

			currentEvent->setX(x);
			currentEvent->setY(y);
			currentEvent->setStamp(time_);
			currentEvent->setPolarity(pol);

			if(DEBUGMODE)
			{
				std::cout << "[debug in Run] : Set current event " << std::endl;
			}

			lastEvent = *currentEvent;

			//cSurf[x][y] = *currentEvent;
			//cSurf.setMostRecent(*currentEvent);

			if(pol == 1)
			{
				surfaceOnL[x][y] = *currentEvent;
        surfaceOfL[x][y] = *currentEvent;
				//cSurf.clear();
				//cSurf.resize(width, height, 0);
				cSurf = surfaceOnL;
			}
			else
			{
				surfaceOfL[x][y] = *currentEvent;
        surfaceOnL[x][y] = *currentEvent;
				//cSurf.resize(width, height, currentEvent);
		    cSurf = surfaceOfL;
			}
			i++;


			//cSurf.setMostRecent(currentEvent);
			//std::cout << "Computing flow" << std::endl;

			//if(pol == 0)
			FlowEvent ofe = computeLocalFlow();

	    if(!std::isnan(abs(ofe.getVx())) && !std::isnan(abs(ofe.getVy())) && ofe.getVx() != 0 && ofe.getVy() != 0) // check that optical flow event is valid
			{

				//std::cout << pol << std::endl;
				//std::cout << "Got flow" << " ";
				//std::cout << "vFlow event occurd " << ofe[0] << std::endl;
				//std::cout << "Percentage complete: " << int((i*100*18.37)/(fsize)) << "%" << " "  << '\r';
				// basic file operations

				double length = sqrt( (ofe.getVx()*ofe.getVx() + ofe.getVy()*ofe.getVy()));
				double theta = atan2(ofe.getVy(), ofe.getVx());

				//std::cout << x << " " << y << " " << pol << " " << length << " " << theta << std::endl;


 			  if(pol == 1) // ON EVENTS
				{
					//cSurf = surfaceOnL;
					flowSurfaceLengthOn[x][y] = length;
					flowSurfaceThetaOn[x][y] = theta;

					flowSurfaceVx[x][y] = ofe.getVx();
					flowSurfaceVy[x][y] = ofe.getVy();

					flowSurfaceLengthOf[x][y] = length;
					flowSurfaceThetaOf[x][y] = theta;

					//double theta_0_2pi = theta - floor(theta/(2*M_PI))*theta;
					//flowSurfaceThetaOn(x, y) = theta_0_2pi;

				}
				else // OFF EVENTS
				{
					//cSurf = surfaceOfL;
					flowSurfaceLengthOn[x][y] = length;
					flowSurfaceThetaOn[x][y] = theta;

					flowSurfaceLengthOf[x][y] = length;
					flowSurfaceThetaOf[x][y] = theta;

					flowSurfaceVx[x][y] = ofe.getVx();
					flowSurfaceVy[x][y] = ofe.getVy();
				}

				FlowEvent trueFlow = computeTrueFlow(x, y, time_, pol);

        double trueLength = sqrt(trueFlow.getVy()*trueFlow.getVy() + trueFlow.getVx()*trueFlow.getVx());
        double trueAngle = atan2(trueFlow.getVy(), trueFlow.getVx());

				//std::cout << "{DEBUG}: Compute true flow (" <<  trueFlow.getVx() << " " << trueFlow.getVy() << ")" << std::endl;

				if(pol == 1)
				{
					eventsFileOut << x << " " << y << " " << time_ << " " << 1 << " " << length << " " << trueAngle << " " << ofe.getVx() << " " << ofe.getVy() << " " << length << " " << theta << std::endl;
					//std::cout << "{DEBUG}: Writing to file " << x << " " << y << " " << time_ << " " << 1 << " " << length << " " << trueFlow.getVy() << " " << ofe.getVx() << " " << ofe.getVy() << " " << length << " " << theta << std::endl;
				}
				else
        {
					eventsFileOut << x << " " << y << " " << time_ << " " << -1 << " " << length << " " << trueAngle << " " << ofe.getVx() << " " << ofe.getVy() << " " << length << " " << theta << std::endl;
					//std::cout << "{DEBUG}: Writing to file " << x << " " << y << " " << time_ << " " << 1 << " " << length << " " << trueFlow.getVy() << " " << ofe.getVx() << " " << ofe.getVy() << " " << length << " " << theta << std::endl;
				}

				 /*
				if(pol == 1)
				{
					flowSurfaceLengthOn[x][y] = trueLength;
					flowSurfaceThetaOn[x][y] = trueAngle;
				}
				else
				{

					flowSurfaceLengthOf[x][y] = trueLength;
					flowSurfaceThetaOf[x][y] = trueAngle;
				}
				// */
			}
			else{

				if(pol ==1)
					eventsFileOut << x << " " << y << " " << time_ << " " << 1 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
				else
					eventsFileOut << x << " " << y << " " << time_ << " " << -1 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;


				 /*
				if(pol == 1)
				{
					flowSurfaceLengthOn[x][y] = 0;
					flowSurfaceThetaOn[x][y] = 0;
				}
				else
				{
					flowSurfaceLengthOf[x][y] = 0;
					flowSurfaceThetaOf[x][y] = 0;
				}
					// */
			}

//			std::cout << aep->getPolarity() << std::endl;

			lastEventTime[x][y] = time_;			
			eventsComputed++;			
			this->numEvents = eventsComputed;
			//std::cout << "Percentage complete: " << double((100*eventsComputed)/(NUMEVENTS)) << "%" << " "  << '\r';
	    }
	    eventsFile.close();
	  }
	 else std::cout << "Unable to open file" << std::endl;

	 eventsFileOut.close();
	 std::cout << std::endl << "Done!" << std::endl;
	 //exit(0);


    //std::string inPortName = "/" + moduleName + "/vBottl:i";
    //bool check1 = BufferedPort<emorph::vBottle>::open(inPortName);

    //if(strictness) outPort.setStrict();
    //std::string outPortName = "/" + moduleName + "/vBottle:o";
    //bool check2 = outPort.open(outPortName);

//    Eigen::MatrixXd q = new Eigen::MatrixXd()
//    vf = new Eigen::VectorXd();
 //       vf->setVx(dtdx);
 //      vf->setVy(dtdy);
 //       vf->setDeath();

    return true; //check1 && check2;


}


/******************************************************************************/
void vFlowManager::close()
{
//    delete surfaceOnL;
//    delete surfaceOfL;
//    delete surfaceOnR;
//    delete surfaceOfR;

}

/******************************************************************************/
#define COS135on2 0.38268
FlowEvent vFlowManager::computeLocalFlow()
{

	if(DEBUGMODE)
	{
		std::cout << "[debug in (compute)] : Entering compute()" << std::endl;
	}

    FlowEvent *vf;
    vf =  new FlowEvent(0, 0, 0, 0, 0, 0);
    double dtdy = 0, dtdx = 0;

    //get the most recent event
	if(DEBUGMODE)
	{
		std::cout << "[debug in (compute)] : made proxy flow event vf" << std::endl;
	}
    //Event vr = cSurf.getMostRecent(); //->getAs<Event>();

	Event vr = lastEvent;

    //find the side of this event that has the collection of temporally nearby
    //events
    double bestscore = MAXSTAMP+1;
    int besti = 0, bestj = 0;

    for(int i = vr.getX()-fRad; i <= vr.getX()+fRad; i+=fRad)
    {
        for(int j = vr.getY()-fRad; j <= vr.getY()+fRad; j+=fRad)
        {
            //get the surface around the recent event
            double sobeltsdiff = 0;
            std::vector<Event> subsurf;

            for(int cx_ = i-fRad; cx_ <= i+fRad; cx_++)
            {
              for(int cy_ = j-fRad; cy_ <= j+fRad; cy_++)
              {
                subsurf.push_back(cSurf[cx_][cy_]);
              }
            } // = cSurf->getSurf(i, j, fRad);

//			std::cout << "[DEBUG in compute]: cSurf of center [cx cy]: " << cSurf[besti][bestj].getX() << " " << std::endl;

//std::cout << "{DEBUG} " << subsurf.size() << " " << planeSize << std::endl;
            if(subsurf.size() < planeSize) continue;

            for(unsigned int k = 0; k < subsurf.size(); k++)
            {
//std::cout << "best event " << subsurf[k]->getPolarity() << " ";
              sobeltsdiff += vr.getStamp() - subsurf[k].getStamp();

//std::cout << " " << vr.getStamp() - subsurf[k].getStamp() << std::endl;
              if(subsurf[k].getStamp() > vr.getStamp())
              {
              //subsurf[k]->setStamp(subsurf[k][2] -
              //                     emorph::vtsHelper::maxStamp());
                sobeltsdiff += MAXSTAMP;
              }
            }
            //std::cout << std::endl;
            sobeltsdiff /= subsurf.size();
            if(sobeltsdiff < bestscore)
            {
              bestscore = sobeltsdiff;
              besti = i; bestj = j;
            }
          }
        }

//	std::cout << "[DEBUG in compute]: bestscore: " << bestscore << std::endl;
        if(bestscore > MAXSTAMP)
        {
          return *vf;
        }


//std::cout << "[DEBUG in compute]: cSurf of center [cx cy]: " << besti << " " << bestj << std::endl;

        std::vector< Event > subsurf;
        for(int cx_ = besti-fRad; cx_ <= besti+fRad; cx_++)
        {
          for(int cy_ = bestj-fRad; cy_ <= bestj+fRad; cy_++)
          {
            subsurf.push_back(cSurf[cx_][cy_]);
          }
        }
	// = cSurf(besti, bestj, fRad);

    //and compute the flow
        if(computeGrads(subsurf, vr, dtdy, dtdx) >= minEvtsOnPlane) //if inliers more than threshold set flow
        {
          //std::cout << "compute grad " << vr[0] << " " << vr[1] << " " << vr->getPolarity() << std::endl;
          vf = new FlowEvent(vr.getX(), vr.getY(), vr.getStamp(), vr.getPolarity(), 0, 0);
          vf->setVx(dtdx);
          vf->setVy(dtdy);
        //vf.setDeath();
        }

        if(DEBUGMODE)
        {
          std::cout << "[debug] : Exiting compute()" << std::endl;
        }
        return *vf;
}

/******************************************************************************/
FlowEvent vFlowManager::computeTrueFlow(int x, int y, unsigned int timeEvent_, int pol)
{

	if(DEBUGMODE)
	{
		std::cout << "[debug] : Entering computeTrueFlow(timed)" << std::endl;
	}

	//std::cout << "in compute true flow time_ based" << std::endl;
	double KILL_OLD_FLOW_TIME = 5000;

	FlowEvent *outFlow =  new FlowEvent(x, y, timeEvent_, pol, 0, 0);
	//outFlow = Eigen::VectorXd(2);

	std::vector<double> spatialPool(maxWindow, 0);

	std::vector<double> spatialTheta(maxWindow, 0);

	std::vector<double> spatialVectorX(maxWindow, 0);

	std::vector<double> spatialVectorY(maxWindow, 0);
	//int windowJump = 1;
	//int maxWindow = 8; // declared as global now


	//std::cout << "On event true copmpute" << std::endl;
	double maxLength = 0;
	double bestTheta = 0;

	double numWindows = 0;

  if(pol == 1)
  {
  	for (int spatialSize = 0; spatialSize <= maxWindow; spatialSize+=windowJump)
  	{

  		double thetaSpatial = 0;
  		double lengthSpatial = 0;
  		double spatialX = 0;
  		double spatialY = 0;

  		double numNeighbors = 0;


  		for(int i = std::max(0, x-spatialSize); i <= std::min(x+spatialSize, width-1) ; i++) //search in spatial range
  		{
  			for(int j = std::max(0, y-spatialSize); j <= std::min(y+spatialSize, width-1); j++)
  			{
  				if(flowSurfaceLengthOn[i][j] > 0 && (abs(timeEvent_-lastEventTime[i][j]) < KILL_OLD_FLOW_TIME))
          //std::cout << "[Debug (computeTrueFlow) ] timeDiff: " << abs(timeEvent_-lastEventTime[i][j]) << std::endl;
  				{
  					lengthSpatial = lengthSpatial + flowSurfaceLengthOn[i][j];

  					spatialX = spatialX + flowSurfaceLengthOn[i][j]*cos(flowSurfaceThetaOn[i][j]); //flowSurfaceVx[i][j] ; //flowSurfaceLengthOn(i, j)*cos(flowSurfaceThetaOn(i, j));
  					spatialY = spatialY + flowSurfaceLengthOn[i][j]*sin(flowSurfaceThetaOn[i][j]); //flowSurfaceVy[i][j] ; //flowSurfaceLengthOn(i, j)*sin(flowSurfaceThetaOn(i, j));

  					thetaSpatial = thetaSpatial + flowSurfaceThetaOn[i][j];

  					numNeighbors++;

  					if(flowSurfaceLengthOn[i][j] > maxLength)
  					{
  						maxLength = flowSurfaceLengthOn[i][j];
  						bestTheta = flowSurfaceThetaOn[i][j];
  					}
  				}
  			}
  		}

  		if(numNeighbors > 0)
  		{
  			spatialPool.at(numWindows) = lengthSpatial/numNeighbors;
  			spatialTheta.at(numWindows) = thetaSpatial/numNeighbors; //bestTheta; //
  			spatialVectorX.at(numWindows) = spatialX/numNeighbors;
  			spatialVectorY.at(numWindows) = spatialY/numNeighbors;
  		}
  		else
  		{
  			spatialPool.at(numWindows) = 0;
  			spatialTheta.at(numWindows) = 0;
  			spatialVectorX.at(numWindows) = 0;
  			spatialVectorY.at(numWindows) = 0;
  		}

  		//std::cout << spatialPool(numWindows) << " " ;
  		numWindows++;
  	}

  	double maxVal = 0;
  	int maxValIndex = 0;



  	for (int spt = 0; spt < numWindows; spt++)
  	{
  		if( spatialPool.at(spt) > maxVal)
  		{
  			maxVal = spatialPool.at(spt);
  			maxValIndex = spt;
  		}
  	}

  	//set flow to max rather than max(mean)
  	//outFlow(0) = maxLength;
  	//outFlow(1) = bestTheta;

  	//use max of mean
  	// /*
  	if(maxVal > 0)
  	{
//  		outFlow->setVx(flowSurfaceLengthOn[x][y]*cos());
//  		outFlow->setVy(atan2(spatialVectorY.at(maxValIndex), spatialVectorX.at(maxValIndex)));

      outFlow->setVx(spatialVectorX.at(maxValIndex));
      outFlow->setVy(spatialVectorY.at(maxValIndex));
      //outFlow->setVx(maxLength);
  		//outFlow->setVy(bestTheta);

  		if(DEBUGMODE)
  		{
  			std::cout <<  "[debug] in computeTrueFlow(timed): Set angle " << " " << atan2(spatialVectorY.at(maxValIndex), spatialVectorX.at(maxValIndex)) << std::endl;
  		}

  	}
  	else
  	{
  		outFlow->setVx(flowSurfaceLengthOn[x][y]*cos(flowSurfaceThetaOn[x][y]));
  		outFlow->setVy(flowSurfaceLengthOn[x][y]*sin(flowSurfaceThetaOn[x][y]));

  		if(DEBUGMODE)
  		{
  			std::cout <<  "[debug] in computeTrueFlow(timed): Set angle " << " " << flowSurfaceLengthOn[x][y] << " " << flowSurfaceThetaOn[x][y] << std::endl;
  		}
  	}
  }

  else //off events

  {
    for (int spatialSize = 0; spatialSize <= maxWindow; spatialSize+=windowJump)
  	{

  		double thetaSpatial = 0;
  		double lengthSpatial = 0;
  		double spatialX = 0;
  		double spatialY = 0;

  		double numNeighbors = 0;


  		for(int i = std::max(0, x-spatialSize); i <= std::min(x+spatialSize, width-1) ; i++) //search in spatial range
  		{
  			for(int j = std::max(0, y-spatialSize); j <= std::min(y+spatialSize, width-1); j++)
  			{
  				if(flowSurfaceLengthOf[i][j] > 0 && (abs(timeEvent_-lastEventTime[i][j]) < KILL_OLD_FLOW_TIME))
          //std::cout << "[Debug (computeTrueFlow) ] timeDiff: " << abs(timeEvent_-lastEventTime[i][j]) << std::endl;
  				{
  					lengthSpatial = lengthSpatial + flowSurfaceLengthOf[i][j];

  					spatialX = spatialX + flowSurfaceLengthOf[i][j]*cos(flowSurfaceThetaOf[i][j]); //flowSurfaceVx[i][j] ; //flowSurfaceLengthOn(i, j)*cos(flowSurfaceThetaOn(i, j));
  					spatialY = spatialY + flowSurfaceLengthOf[i][j]*sin(flowSurfaceThetaOf[i][j]); //flowSurfaceVy[i][j] ; //flowSurfaceLengthOn(i, j)*sin(flowSurfaceThetaOn(i, j));

  					thetaSpatial = thetaSpatial + flowSurfaceThetaOf[i][j];
  					numNeighbors++;

  					if(flowSurfaceLengthOf[i][j] > maxLength)
  					{
  						maxLength = flowSurfaceLengthOf[i][j];
  						bestTheta = flowSurfaceThetaOf[i][j];
  					}
  				}
  			}
  		}

  		if(numNeighbors > 0)
  		{
  			spatialPool.at(numWindows) = lengthSpatial/numNeighbors;
  			spatialTheta.at(numWindows) = thetaSpatial/numNeighbors; //bestTheta; //
  			spatialVectorX.at(numWindows) = spatialX/numNeighbors;
  			spatialVectorY.at(numWindows) = spatialY/numNeighbors;
  		}
  		else
  		{
  			spatialPool.at(numWindows) = 0;
  			spatialTheta.at(numWindows) = 0;
  			spatialVectorX.at(numWindows) = 0;
  			spatialVectorY.at(numWindows) = 0;
  		}

  		//std::cout << spatialPool(numWindows) << " " ;
  		numWindows++;
  	}

  	double maxVal = 0;
  	int maxValIndex = 0;

  	for (int spt = 0; spt < numWindows; spt++)
  	{
  		if( spatialPool.at(spt) > maxVal)
  		{
  			maxVal = spatialPool.at(spt);
  			maxValIndex = spt;
  		}
  	}

  	//set flow to max rather than max(mean)
  	//outFlow(0) = maxLength;
  	//outFlow(1) = bestTheta;

  	//use max of mean
  	// /*
  	if(maxVal > 0) // if there was flow activity in a spatial scale
  	{
  		outFlow->setVx(spatialVectorX.at(maxValIndex));
  		outFlow->setVy(spatialVectorY.at(maxValIndex));

      //outFlow->setVx(maxLength);
  		//outFlow->setVy(bestTheta);

  		if(DEBUGMODE)
  		{
  			std::cout <<  "[debug] in computeTrueFlow(timed): Set angle " << " " << atan2(spatialVectorY.at(maxValIndex), spatialVectorX.at(maxValIndex)) << std::endl;
  		}

  	}
  	else  // if now assign the same as original flow
  	{
  		outFlow->setVx(flowSurfaceLengthOf[x][y]*cos(flowSurfaceThetaOf[x][y]));
  		outFlow->setVy(flowSurfaceLengthOf[x][y]*sin(flowSurfaceThetaOf[x][y]));

  		if(DEBUGMODE)
  		{
  			std::cout <<  "[debug] in computeTrueFlow(timed): Set angle " << " " << flowSurfaceLengthOn[x][y] << " " << flowSurfaceThetaOn[x][y] << std::endl;
  		}
  	}

  }
	// */

	//std::cout << maxVal << " " << maxValIndex << " " << std::endl;

	if(DEBUGMODE)
	{
		std::cout << "[debug] : Exiting computeTrueFlow(timed)" << std::endl;
	}
	return *outFlow;

}


/******************************************************************************/
int vFlowManager::computeGrads(std::vector<Event> subsurf,
                                     Event &cen,
                                     double &dtdy, double &dtdx)
{

//std::cout << "[Debug (computeGrads_1)]: size of Event vectors: " << subsurf.size() << std::endl;

    Eigen::MatrixXd A(subsurf.size(), 3);
//std::cout << "[Debug (computeGrads_1)]: size of A matrix : " << A.rows() << " " << A.cols() << std::endl;
    Eigen::VectorXd Y(subsurf.size());
    for(unsigned int vi = 0; vi < subsurf.size(); vi++) {
        Event v = subsurf.at(vi);
        A(vi, 0) = v.getX();
        A(vi, 1) = v.getY();
        A(vi, 2) = 1;
        if(v.getStamp() > cen.getStamp()) {
            Y(vi) = (v.getStamp() - this->MAXSTAMP) * this->TSTOSEC;
        } else {
            Y(vi) = v.getStamp() * this->TSTOSEC;
        }
    }

    return computeGrads(A, Y, cen.getX(), cen.getY(), cen.getStamp() *
                        this->TSTOSEC, dtdy, dtdx);
}

/******************************************************************************/
int vFlowManager::computeGrads(Eigen::MatrixXd A, Eigen::MatrixXd Y,
                                     double cx, double cy, double cz,
                                     double &dtdy, double &dtdx)
{
  double meanX=0, meanY=0, meanT=0;
  double sx2=0, sy2=0, sxy=0, stx=0, sty=0;

	Eigen::VectorXd iV(2);

	double tol = 0;

// get mean of point cloud
 /*    for (int i = 0; i < A.rows(); i++)
    {
		meanX = meanX + A(i, 0);
		meanY = meanY + A(i,1);
		meanT = meanT + Y(i);

		sx2 = sx2 + A(i,0)*A(i,0);
		sy2 = sy2 + A(i,1)*A(i,1);
		sxy = sxy + A(i,0)*A(i,1);

		stx = stx + A(i,0)*Y(i);
		sty = sty + A(i,1)*Y(i);
	}
	meanX = meanX/(A.rows()+1);
	meanY = meanY/(A.rows()+1);
	meanT = meanT/(A.rows()+1);
	*/
// center events to [0 0 0]
 /*
	for (int i = 0; i < A.rows(); i++)
    {
		A(i, 0) = A(i, 0) - cx; //meanX
		A(i, 1) = A(i, 1) - cy; // meanY
		Y(i) = Y(i) - cz; // meanT
	}
//	*/
// calculate covariance
 /*
    for (int i = 0; i < A.rows(); i++)
    {
		sx2 = sx2 + A(i,0)*A(i,0);
		sy2 = sy2 + A(i,1)*A(i,1);
		sxy = sxy + A(i,0)*A(i,1);

		stx = stx + A(i,0)*Y(i);
		sty = sty + A(i,1)*Y(i);
	}


	double detmnt = sx2*sy2-sxy*sxy;

	if(detmnt != tol)
	{
		iV(0) = (sy2*stx-sxy*sty)/detmnt;
		iV(1) = (sx2*sty-sxy*stx)/detmnt;
		abc(0) = -1*iV(0);
		abc(1) = -1*iV(1);
	}
// */

// /*

//std::cout << "[Debug (computeGrads_2)]: num rows A before Transpose : " << A.rows() << std::endl;

  At=A.transpose();

//std::cout << "[Debug (computeGrads_2)]: num rows A after Transpose : " << A.rows() << std::endl;

  AtA=At*A;



  double* dataATA=AtA.data();
  double DET = AtA.determinant();
  //*dataATA*( *(dataATA+8)**(dataATA+4)-*(dataATA+7)**(dataATA+5))-
    //        *(dataATA+3)*(*(dataATA+8)**(dataATA+1)-*(dataATA+7)**(dataATA+2))+
      //      *(dataATA+6)*(*(dataATA+5)**(dataATA+1)-*(dataATA+4)**(dataATA+2));

//  std::cout << "[Debug: (computeGrads_2)]: " << AtA.determinant() << " " << DET << std::endl;

  if(DET < 1) return 0;

//  A2 = AtA;
  double *dataA=A2.data();
  DET=1.0/DET;
    *(dataA+0)=DET*(*(dataATA+8)**(dataATA+4)-*(dataATA+7)**(dataATA+5));
    *(dataA+1)=DET*(*(dataATA+7)**(dataATA+2)-*(dataATA+8)**(dataATA+1));
    *(dataA+2)=DET*(*(dataATA+5)**(dataATA+1)-*(dataATA+4)**(dataATA+2));
    *(dataA+3)=DET*(*(dataATA+6)**(dataATA+5)-*(dataATA+8)**(dataATA+3));
    *(dataA+4)=DET*(*(dataATA+8)**(dataATA+0)-*(dataATA+6)**(dataATA+2));
    *(dataA+5)=DET*(*(dataATA+3)**(dataATA+2)-*(dataATA+5)**(dataATA+0));
    *(dataA+6)=DET*(*(dataATA+7)**(dataATA+3)-*(dataATA+6)**(dataATA+4));
    *(dataA+7)=DET*(*(dataATA+6)**(dataATA+1)-*(dataATA+7)**(dataATA+0));
    *(dataA+8)=DET*(*(dataATA+4)**(dataATA+0)-*(dataATA+3)**(dataATA+1));

	Eigen::MatrixXd temp = A2*At*Y;
	abc.at(0)=temp(0);
  abc.at(1)=temp(1);
  abc.at(2)=temp(2);
// */

//std::cout << "[Debug (computeGrads_2)]: size of A matrix : " << A.rows() << " " << A.cols() << std::endl;
//std::cout << "[Debug (computegradient)]: Size of A: " << A.size() << std::endl;

//std::cout << "[Debug (computeGrads_2)]: num rows A after Transpose : " << A.rows() << std::endl;

  double dtdp = sqrt(pow(abc.at(0), 2.0) + pow(abc.at(1), 2.0));
  int inliers = 0;
//std::cout << "[Debug (computeGrads_2)]: num rows A : " << A.rows() << std::endl;
  for(int i = 0; i < A.rows(); i++)
  {
        //so I think that abc(0) and abc(1) are already scaled to the magnitude
        //of the slope of the plane. E.g. when only using abc(0) and abc(1) and
        //fitting a 3-point plane we always get 0 error. Therefore the differ-
        //ence in time is perfect with only abc(0,1) and the speed should also
        //be.

        double planedt = (abc.at(0) * (A(i, 0) - cx) + abc.at(1) * (A(i, 1) - cy));
        double actualdt =  Y(i) - cz;

        //double planedt = (abc.at(0) * (A(i, 0) ) + abc.at(1) * (A(i, 1) ));
        //double actualdt =  Y(i) ;

        if(abs(planedt - actualdt) < dtdp/2 && Y(i) > 0) inliers++;

//        std::cout << "[Debug (computeGrads_2)]: num inliers: " << inliers << " " << cx << " " << cy << " " << cz << " " << Y(i) << " " << planedt << " " << actualdt << " " << fabs(planedt - actualdt) << " " << dtdp/2 << " " << std::endl;
  }



  double speed = 1.0 / dtdp;

  double angle = atan2(abc.at(0), abc.at(1));
  dtdx = speed * cos(angle);
  dtdy = speed * sin(angle);

  return inliers;

}
