#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

//simulation parameters
const int N_x = 100; //how many points there are on the 1D grid
const int N_steps = 500; //amount of time steps in the simulation
const int steps_per_data_write = 50; //amount of steps transpiring per each time the data is saved

const double domain_size = 10.0; //size of the physical domain
const double dx = domain_size/N_x; //distance between grid points on the physical domain (spatial step size)
const double dt = 0.9 * dx; //time step 

//since we're working in natural units, c = h_bar = 1, so the only numerical value in the Klein-Gordon equation we need to worry about is mass
const double m = 1.0; //mass

//Gaussian parameters
const double A = 1.0; //amplitude of wave packet
const double x_0 = domain_size/2.0; //initial position of the wave packet, initialized to center of domain
const double sigma = 0.5;  //Gaussian width 
const double k = 3.0; //wave number


void set_initial_conditions(std::vector<double>& phi, std::vector<double>& prev_phi)
{
    for(int i = 0; i < N_x; i++)
    {
        double x_i = i * dx; //physical coordinate which corresponds to the ith grid point
        phi[i] = A * std::exp(-1 * std::pow((x_i - x_0), 2)/(2.0 * sigma * sigma)) * std::cos(k * x_i); //initialization of each field point based on gaussian wave packet
        prev_phi[i] = phi[i]; //each future time step requires knowledge of the previous and present time steps 
    }
}

void time_step_evolution(std::vector<double>& phi, std::vector<double>& prev_phi, std::vector<double>& new_phi)
{   
    //computes the next time step for each field point excluding the boundaries to enforce periodic boundary conditions
    for(int i = 1; i < N_x - 1; i++)
    {
        double laplacian = (phi[i - 1] - (2.0 * phi[i]) + phi[i + 1])/(dx * dx); //computes the approximate laplacian given by the finite difference method
        new_phi[i] = 2.0 * phi[i] - prev_phi[i] + (dt * dt) * (laplacian - (m * m) * phi[i]); //next field point calculation
    }

    //since we want to enforce periodic boundary conditions it is necessary to make sure the edges of the field wrap around
    new_phi[0] = new_phi[N_x - 2];
    new_phi[N_x - 1] = new_phi[1]; 

    //after one time step the previous field becomes the "current" one and the "current" one becomes the one that was calculated in this function
    prev_phi = phi;
    phi = new_phi;

}

void write_data(const std::vector<double>& phi, int step) 
{
    std::ofstream file("klein_gordon_output_" + std::to_string(step) + ".dat");
    for (int i = 0; i < N_x; i++) 
    {
        file << i * dx << " " << phi[i] << "\n";
    }
    file.close();
}

int main() 
{
    //field arrays 
    std::vector<double> prev_phi(N_x, 0.0);
    std::vector<double> phi(N_x, 0.0);
    std::vector<double> new_phi(N_x, 0.0);

    //sets initial conditions 
    set_initial_conditions(phi, prev_phi);

    //time evolution loop
    for(int step = 0; step < N_steps; step++)
    {
        time_step_evolution(phi, prev_phi, new_phi);
        std::cout << "completed step" << step << "\n";

        //writes data for each step that is steps_per_data_write steps away from the previous data write
        if(step % steps_per_data_write == 0)
        {
            write_data(phi, step);
            std::cout << "wrote data for step" << step << "\n";
        }
    }
}

