/*
This program performs a simple Particle-Particle N-body simulation
Initial conditions for position and velocity are initially randomized over a small range
Position and velocity data is written to a csv file on each time step for later analysis
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define NUM_PARTICLES 50

//Constant values, dm_mass and epsilon currently chosen to give decent results
const double dm_mass = 1e6; //kg
const double G = 6.6743015e-11; //N*m^2/kg^2
const double epsilon = 1e-1; //m
const double dt = 0.1;
const double tf = 10000.0;
const double sigmaPos = 5.0;
const double sigmaVel = 1e-7;

//Each particle is represented in a struct that holds position and velocity values
struct particle {
    double position[3];
    double velocity[3];
};

double randfrom(double min, double max);
void initialize(struct particle *parts);
double fv(struct particle *parts, double k, double multiplier, int dim, int part_num);
double sum_term_func(struct particle *parts, int part_num, int dim);

int main() {
    srand(time(0)); //Seed the random number generator

    //First call to function gives consistent value
    double rand_init = randfrom(-1.0, 1.0);

    //Runge Kutta variables
    double k1_p;
    double k2_p;
    double k3_p;
    double k4_p;

    double k1_v;
    double k2_v;
    double k3_v;
    double k4_v;

    //Analysis file
    FILE *fp = NULL;
    fp = fopen("Data\\pos_data.csv", "w");
    fprintf(fp, "time,particle,px,py,pz,vx,vy,vz\n");

    struct particle particles[NUM_PARTICLES]; //Main struct containing particle data
    struct particle parts_copy[NUM_PARTICLES]; //Copy of main struct

    //Normally initializing the particles
    initialize(particles);

    //Main time loop
    double t = 0;
    while (t < tf) {
        if ((int)t % 100 == 0) {
            printf("Time: %f\n", t);
        }
        //Dumping data into file on each iteration
        for (int i=0; i < NUM_PARTICLES; i++) {
            fprintf(fp, "%f,%d,%f,%f,%f,%f,%f,%f\n", t, i, particles[i].position[0],
                    particles[i].position[1], particles[i].position[2], 
                    particles[i].velocity[0], particles[i].velocity[1],
                    particles[i].velocity[2]);
        }

        //Updating copy on each iteration
        for (int i = 0; i < NUM_PARTICLES; i++) {
            parts_copy[i] = particles[i];
        }

        //4th order Runge Kutta on each particle-dimension pair
        for (int i = 0; i < NUM_PARTICLES; i++) {
            for (int j = 0; j < 3; j++) {
                k1_p = dt*parts_copy[i].velocity[j]; //m
                k2_p = dt*(parts_copy[i].velocity[j] + 0.5*k1_p);
                k3_p = dt*(parts_copy[i].velocity[j] + 0.5*k2_p);
                k4_p = dt*(parts_copy[i].velocity[j] + k3_p);
                
                //Updating positions in main struct
                particles[i].position[j] += (k1_p + 2*k2_p + 2*k3_p + k4_p)/6;

                k1_v = dt*fv(parts_copy, 0, 0, j, i); //m/s
                k2_v = dt*fv(parts_copy, k1_v, 0.5, j, i);
                k3_v = dt*fv(parts_copy, k2_v, 0.5, j, i);
                k4_v = dt*fv(parts_copy, k3_v, 1, j, i);

                //Updating velocities in main struct
                particles[i].velocity[j] += (k1_v + 2*k2_v + 2*k3_v + k4_v)/6;
            }
        }

        t += dt;
    }

    //Final dump for last iteration
    for (int i=0; i < NUM_PARTICLES; i++) {
        fprintf(fp, "%f,%d,%f,%f,%f,%f,%f,%f\n", t, i, particles[i].position[0],
            particles[i].position[1], particles[i].position[2], 
            particles[i].velocity[0], particles[i].velocity[1],
            particles[i].velocity[2]);
    }

    fclose(fp);
    return 0;
}

//Generates a random double in a range
double randfrom(double min, double max) {
    double range = max - min;
    double div = RAND_MAX/range;
    return min + rand()/div;
}

void initialize(struct particle *parts) {
    for (int i=0; i < NUM_PARTICLES; i++) {
        for (int j=0; j < 3; j++) {
            double rand1 = randfrom(0.0, 1.0);
            double rand2 = randfrom(0.0, 1.0);
            double randNorm1 = sqrt(-2*log(rand1))*cos(2*M_PI*rand2);
            double randNorm2 = sqrt(-2*log(rand1))*sin(2*M_PI*rand2);
            parts[i].position[j] = sigmaPos*randNorm1;
            parts[i].velocity[j] = sigmaVel*randNorm2;
        }
    }
}

//Updates inputs according to the Runge Kutta step and performs acceleration calculation
double fv(struct particle *parts, double k, double multiplier, int dim, int part_num) {
    struct particle parts_use[NUM_PARTICLES]; //Copy of the copy

    //Updating each particle's position for the Runge Kutta step
    for (int i = 0; i < NUM_PARTICLES; i++) {
        parts_use[i] = parts[i];
        parts_use[i].position[dim] = parts[i].position[dim] + k*multiplier;
    }

    //Acceleration calculation
    return -G*dm_mass*sum_term_func(parts_use, part_num, dim); //m/s^2
}

//Superposition of contributions from each particle
double sum_term_func(struct particle *parts, int part_num, int dim) {
    double sum_term = 0;
    double numerator;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        if (i == part_num) {
            continue;
        }
        
        //Calculating distances
        double dx = parts[part_num].position[0] - parts[i].position[0];
        double dy = parts[part_num].position[1] - parts[i].position[1];
        double dz = parts[part_num].position[2] - parts[i].position[2];
        double r = sqrt(dx*dx + dy*dy + dz*dz + epsilon*epsilon);

        //Setting appropriate numerator for the given dimension
        switch(dim) {
            case 0:
                numerator = dx;
                break;
            case 1:
                numerator = dy;
                break;
            case 2:
                numerator = dz;
                break;
            default:
                //Particle does not contribute if there is an issue
                numerator = 0;
        }

        //Adding individual term of the superposition
        sum_term += numerator/(r*r*r);
    }

    return sum_term;
}