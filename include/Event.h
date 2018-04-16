/* 
 * Author: himanshu.akolkar@gmail.com
 * 
*/


#ifndef __EVENT__
#define __EVENT__

#include <string>
#include <math.h>
///
///
class Event
{
private:    

int x;
int y;
int pol;
double t;


public:

    Event(int x, int y, double t, int p);
	Event();
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
    
    
    int getX();
    int getY();
    double getStamp();
    int getPolarity();

};

/******************************************************************************/
#endif //__EVENT__

