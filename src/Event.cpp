/* 
 * Author: himanshu.akolkar@gmail.com
 * 
*/

#include "../include/Event.h"
///
///
Event::Event(int x, int y, double t, int pol){

this->x = x;
this->y = y;
this->t = t;
this->pol = pol;
};

Event::Event(){

this->x = 0;
this->y = 0;
this->t = 0;
this->pol = 0;
};

bool Event::setX(int value){
	this->x = value;
	return true;
};

bool Event::setY(int value){
	this->y = value;
	return true;
};

bool Event::setStamp(double value){

	this->t = value;
	return true;
	
	};

bool Event::setPolarity(int value){

	this->pol = value;
	return true;
};
    
    
int Event::getX(){ return x; };
int Event::getY(){ return y; };
double Event::getStamp(){ return t; };
int Event::getPolarity(){ return pol; };


