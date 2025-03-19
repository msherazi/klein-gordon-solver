#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

//simulation parameters
const int N_x = 500; //how many points there are on the 1D grid
const int N_steps = 500; //amount of time steps in the simulation
const int steps_per_data_write = 50; //amount of steps transpiring per each time the data is saved

const double domain_size = 10.0; //size of the physical domain
const double dx = domain_size / N_x; //distance between grid points on the physical domain (spatial step size)
const double dt = 0.9 * dx; //time step

//since we're working in natural units, c = h_bar = 1,
//so the only numerical value in the Klein-Gordon equation we need to worry about is mass
const double m = 5.0; // mass

//Gaussian parameters
const double A = 1.0; // amplitude of wave packet
const double x_0 = domain_size / 2.0; // initial position of the wave packet, initialized to center of domain
const double sigma = 0.5;  // Gaussian width
const double k = 3.0; // wave number


//error checking function
static void checkCudaError(cudaError_t result, const char* msg)
{
    if (result != cudaSuccess) {
        std::cerr << "CUDA error (" << msg << "): "
                  << cudaGetErrorString(result) << std::endl;
        exit(EXIT_FAILURE);
    }
}

void set_initial_conditions(std::vector<double>& phi, std::vector<double>& prev_phi)
{
    for (int i = 0; i < N_x; i++)
    {
        double x_i = i * dx; //physical coordinate which corresponds to the i-th grid point
        phi[i] = A * std::exp(-1.0 * std::pow((x_i - x_0), 2) / (2.0 * sigma * sigma)) * std::cos(k * x_i); // initialization of each field point based on Gaussian wave packet
        prev_phi[i] = phi[i]; // each future time step requires knowledge of the previous and present time steps
    }
}

//CUDA kernel that calculates each time step
__global__ void time_step_kernel(const double* __restrict__ phi, const double* __restrict__ prev_phi, double* __restrict__ new_phi, double dx, double dt, double m, int N_x)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= 1 && i < (N_x - 1))
    {
        double laplacian = (phi[i - 1] - 2.0 * phi[i] + phi[i + 1]) / (dx * dx); //computes the approximate laplacian given by the finite difference method

        new_phi[i] = 2.0 * phi[i] - prev_phi[i] + (dt * dt) * (laplacian - (m * m) * phi[i]); //next field point calculation
    }
}


//CUDA kernel that applies periodic boundary conditions
__global__ void apply_periodic_bc_kernel(double* __restrict__ new_phi, int N_x)
{
    //since we want to enforce periodic boundary conditions it is necessary to make sure the edges of the field wrap around
    new_phi[0] = new_phi[N_x - 2];
    new_phi[N_x - 1] = new_phi[1];
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
    std::vector<double> h_phi(N_x, 0.0);
    std::vector<double> h_prev_phi(N_x, 0.0);
    std::vector<double> h_new_phi(N_x, 0.0);

    //sets initial conditions
    set_initial_conditions(h_phi, h_prev_phi);

    //pointer declaration
    double* d_phi       = nullptr;
    double* d_prev_phi  = nullptr;
    double* d_new_phi   = nullptr;


    //gpu memory allocation
    checkCudaError(cudaMalloc((void**)&d_phi, N_x * sizeof(double)), "alloc d_phi");
    checkCudaError(cudaMalloc((void**)&d_prev_phi, N_x * sizeof(double)), "alloc d_prev_phi");
    checkCudaError(cudaMalloc((void**)&d_new_phi, N_x * sizeof(double)), "alloc d_new_phi");

    //copies data from the CPU to GPU
    checkCudaError(cudaMemcpy(d_phi, h_phi.data(), N_x * sizeof(double), cudaMemcpyHostToDevice), "cpy phi H->D");
    checkCudaError(cudaMemcpy(d_prev_phi, h_prev_phi.data(), N_x * sizeof(double), cudaMemcpyHostToDevice), "cpy prev_phi H->D");

    
    //kernel parameters
    const int blockSize = 128;
    const int gridSize  = (N_x + blockSize - 1)/blockSize;

    //time evolution loop
    for (int step = 0; step < N_steps; step++)
    {
        //launches time_step_kernel
        time_step_kernel<<<gridSize, blockSize>>>(d_phi, d_prev_phi, d_new_phi, dx, dt, m, N_x);
        checkCudaError(cudaGetLastError(), "time_step_kernel");

        //launches periodic boundary conditions kernel
        //<<<1, 1>>> as arguments because only two values are being updated 
        apply_periodic_bc_kernel<<<1, 1>>>(d_new_phi, N_x);
        checkCudaError(cudaGetLastError(), "apply_periodic_bc_kernel");

        //makes sure that everything is working
        checkCudaError(cudaDeviceSynchronize(), "kernel sync");

        
        //variable swapping
        double* temp = d_prev_phi;
        d_prev_phi = d_phi;
        d_phi      = d_new_phi;
        d_new_phi  = temp;
    

        std::cout << "completed step " << step << "\n";
     
        if (step % steps_per_data_write == 0 || step == N_steps - 1)
        {
            checkCudaError(cudaMemcpy(h_phi.data(), d_phi, N_x * sizeof(double), cudaMemcpyDeviceToHost), "cpy phi D->H");

            write_data(h_phi, step);
            std::cout << "wrote data for step " << step << "\n";
        }
    }

    //cleans up device memory
    cudaFree(d_phi);
    cudaFree(d_prev_phi);
    cudaFree(d_new_phi);

    return 0;
}
