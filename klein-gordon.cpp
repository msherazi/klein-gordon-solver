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


void initial_conditions(std::vector<double>& phi, std::vector<double>& prev_phi)
{
    for(int i = 0; i < N_x; i++)
    {
        double x_i = i * dx; //physical coordinate which corresponds to the ith grid point
        phi[i] = A * std::exp(-1 * std::pow((x_i - x_0), 2)/(2.0 * sigma * sigma)) * std::cos(k * x_i); //initialization of each field point based on gaussian wave packet
        prev_phi[i] = phi[i]; //each future time step requires knowledge of the previous and present time steps 
    }
}

void time_evolution(std::vector<double>& phi, std::vector<double>& prev_phi, std::vector<double>& prev_new)
{

}