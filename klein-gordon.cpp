#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

//simulation parameters
const int N_x = 100; //how many points there are on the 1D grid
const double domain_size = 10.0; //size of the physical domain
const double dx = domain_size/N_x; //distance between grid points on the physical domain (spatial step size)
const double dt = 0.9 * dx;

//since we're working in natural units, c = h_bar = 1, so the only numerical value we need to worry about is mass
const double m = 1.0; //mass

//gaussian parameters
const double A = 1.0; //amplitude of wave packet
const double x_0 = domain_size/2.0; //initial position of the wave packet, initialized to center of domain
const double sigma = 0.5;  //gaussian width 
const double k = 3.0; //wave number


void initial_conditions()
{

}
