/* 
 * Author: himanshu.akolkar@gmail.com
 * 
*/


#ifndef __FLOWEVENT__
#define __FLOWEVENT__

#include <string>
#include <math.h>
///
///
class FlowEvent
{
private:    

int x;
int y;
int pol;
double t;
double Vx;
double Vy;
int scale;

public:

    FlowEvent(int x, int y, double t, int p, double Vx, double Vy);
    FlowEvent();

	FlowEvent operator=(const FlowEvent&);

    ///
    /// \brief open the port and start callback
    /// \param moduleName defines port names (default vFlow)
    /// \param strictness read and write strcitly such that events are not lost
    /// \return true on successful open
    ///    
    bool setX(int);
    bool setY(int);
    bool setStamp(double);
    bool setPolarity(int);
    bool setVx(double);
    bool setVy(double);
    bool setScale(int);
    
    int getX();
    int getY();
    double getVx();
    double getVy();
    double getStamp();
    int getPolarity();
    int getScale();

};

/******************************************************************************/
#endif //__FLOWEVENT__

