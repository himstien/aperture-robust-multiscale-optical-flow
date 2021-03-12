/* 
 * Author: himanshu.akolkar@gmail.com
 * 
*/

#include "../include/FlowEvent.h"
///
///
FlowEvent::FlowEvent(int x, int y, double t, int p, double Vx, double Vy){

this->x = x;
this->y = y;
this->pol = pol;
this->t = t;
this->Vx = Vx;
this->Vy = Vy;
this->scale = 0;

};


FlowEvent::FlowEvent(){

this->x = 0;
this->y = 0;
this->pol = 0;
this->t = 0;
this->Vx = 0;
this->Vy = 0;
this->scale = 0;

};

bool FlowEvent::setX(int value){

 this->x = value;
	
}

bool FlowEvent::setY(int value){

 this->y = value;
	
}

bool FlowEvent::setStamp(double value){

 this->t = value;
	
}

bool FlowEvent::setPolarity(int value){

 this->pol = value;
 return true;
	
}

bool FlowEvent::setVx(double value){

 this->Vx = value;
 return true;
 
}

bool FlowEvent::setVy(double value){

 this->Vy = value;
 return true;	
}

bool FlowEvent::setScale(int value){

 this->scale = value;
 return true;	
}
      
FlowEvent FlowEvent::operator=(const FlowEvent& copy) {

FlowEvent tmp;
tmp.setX(copy.x);
tmp.setY(copy.y);
tmp.setStamp(copy.t);
tmp.setPolarity(copy.pol);
tmp.setVx(copy.Vx);
tmp.setVy(copy.Vy);
tmp.setScale(copy.scale);
return tmp;

}
    
    
int FlowEvent::getX(){ return x;};
int FlowEvent::getY(){ return y;};
double FlowEvent::getVx(){ return Vx;};
double FlowEvent::getVy(){ return Vy;};
double FlowEvent::getStamp(){ return t;};
int FlowEvent::getPolarity(){ return pol;};
int FlowEvent::getScale(){ return scale;};

