// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>


double t = 0;
double tFinal = 0;
double tPlot = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;
//todo delete
bool collided = false;


/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double **x;

/**
 * Equivalent to x storing the velocities.
 */
double **v;

/**
 * One mass entry per molecule/particle.
 */
double *mass;

/**
 * Global time step size used.
 */
double timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double maxV;

/**
 * Minimum distance between two elements.
 */
double minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char **argv) {
    NumberOfBodies = (argc - 4) / 7;

    x = new
    double*[NumberOfBodies];
    v = new
    double*[NumberOfBodies];
    mass = new
    double [NumberOfBodies];

    int readArgument = 1;

    tPlotDelta = std::stof(argv[readArgument]);
    readArgument++;
    tFinal = std::stof(argv[readArgument]);
    readArgument++;
    timeStepSize = std::stof(argv[readArgument]);
    readArgument++;

    for (int i = 0; i < NumberOfBodies; i++) {
        x[i] = new
        double[3];
        v[i] = new
        double[3];

        x[i][0] = std::stof(argv[readArgument]);
        readArgument++;
        x[i][1] = std::stof(argv[readArgument]);
        readArgument++;
        x[i][2] = std::stof(argv[readArgument]);
        readArgument++;

        v[i][0] = std::stof(argv[readArgument]);
        readArgument++;
        v[i][1] = std::stof(argv[readArgument]);
        readArgument++;
        v[i][2] = std::stof(argv[readArgument]);
        readArgument++;

        mass[i] = std::stof(argv[readArgument]);
        readArgument++;

        if (mass[i] <= 0.0) {
            std::cerr << "invalid mass for body " << i << std::endl;
            exit(-2);
        }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

    if (tPlotDelta <= 0.0) {
        std::cout << "plotting switched off" << std::endl;
        tPlot = tFinal + 1.0;
    } else {
        std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
        tPlot = 0.0;
    }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
    videoFile.open("result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"
              << std::endl
              << "<Collection>";
}


/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
    static int counter = -1;
    counter++;
    std::stringstream filename;
    filename << "result-" << counter << ".vtp";
    std::ofstream
    out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int i = 0; i < NumberOfBodies; i++) {
        out << x[i][0]
            << " "
            << x[i][1]
            << " "
            << x[i][2]
            << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>"
              << std::endl;
}


/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {


//    std::cout << "Position 1: " << x[0][0] << " " << x[0][1] << " " << x[0][2] << std::endl;
//    std::cout << "Position 2: " << x[1][0] << " " << x[1][1] << " " << x[1][2] << std::endl;
//    std::cout << "Position 3: " << x[2][0] << " " << x[2][1] << " " << x[2][2] << std::endl;

    if (NumberOfBodies == 1) {
        std::cout << "Position: " << x[0][0] << " " << x[0][1] << " " << x[0][2] << std::endl;
        tPlot = tFinal;
        tFinal = t;
        if (tPlotDelta != 0) {
            // Plot the final timestep if plotting is active
            t += tPlot;
        }
        t += timeStepSize;
        return;
    }

    // 10 total buckets
    int totalBuckets = 10;
    // Pointer that points to each bucket
    int **buckets = new
    int*[totalBuckets];
    // Each bucket can hold up to NumberOfBodies values in it
    for (int i = 0; i < totalBuckets; i++) {
        buckets[i] = new
        int[NumberOfBodies];
    }
    int *bucketLocation = new
    int[totalBuckets]();
    if (t == 0) {
        for (int particle = 0; particle < NumberOfBodies; particle++) {
            maxV = std::max(maxV, std::sqrt(v[particle][0] * v[particle][0] + v[particle][1] * v[particle][1] +
                                            v[particle][2] * v[particle][2]));
        }
    }

    // The velocity difference between each bucket
    double vBucket = maxV / totalBuckets;

    for (int i = 0; i < NumberOfBodies; i++) {
        // Get the velocity of each particle and sort it into the correct bucket
        double velocity = std::sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
        int bucketChosen = std::floor(velocity / vBucket);
        if (bucketChosen < 0) {
            bucketChosen = 0;
        }

        if (bucketChosen == totalBuckets) {
            buckets[bucketChosen - 1][bucketLocation[bucketChosen - 1]++] = i;
        } else {
            buckets[bucketChosen][bucketLocation[bucketChosen]++] = i;
        }
    }



    if (collided == 1){

        std::cout << "Checking start" << std::endl;
        exit(0);

    }
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();

    // force0 = force along x direction
    // force1 = force along y direction
    // force2 = force along z direction
    double *force0 = new
    double[NumberOfBodies]();
    double *force1 = new
    double[NumberOfBodies]();
    double *force2 = new
    double[NumberOfBodies]();


//todo check if this works
    double *force0New = new
    double[NumberOfBodies]();
    double *force1New = new
    double[NumberOfBodies]();
    double *force2New = new
    double[NumberOfBodies]();


    for (int bucket = 0; bucket < totalBuckets; bucket++) {
        if (bucketLocation[bucket] == 0) {
            continue;
        }
        int totalSteps = std::pow(2, bucket);
        double newTimeStepSize = timeStepSize / totalSteps;
        for (int step = 0; step < totalSteps; step++) {
            for (int i = 0; i < bucketLocation[bucket]; i++) {
                int particle = buckets[bucket][i];

                //todo go through particles in the order they are in in the buckets
                force0New[particle] = force0[particle];
                force1New[particle] = force1[particle];
                force2New[particle] = force2[particle];
//                std::cout << "Forces on " << particle << " " << force0New[particle] << " " << force1New[particle] << " "
//                          << force2New[particle] << std::endl;

                for (int bucketOfJ = bucket; bucketOfJ < totalBuckets; bucketOfJ++) {
                    int start = 0;
                    if ((bucket == bucketOfJ) && (i != bucketLocation[bucket] - 1)) {
                        start = i + 1;
                    }
                    if ((bucket != bucketOfJ) || (start != 0)) {
                        for (int j = start; j < bucketLocation[bucketOfJ]; j++) {
                            int jParticle = buckets[bucketOfJ][j];

                            if (collided == 1){
                                std::cout << "Checking particle " << particle << " with " << jParticle << std::endl;

                            }
//                            std::cout << "Comparing i=" << particle << " with j=" << jParticle << std::endl;
                            double dist0 = x[jParticle][0] - x[particle][0], dist1 =
                                    x[jParticle][1] - x[particle][1], dist2 =
                                    x[jParticle][2] - x[particle][2];
                            double squareDistance = dist0 * dist0 + dist1 * dist1 + dist2 * dist2;

                            while ((squareDistance <= (0.01 * 0.01)) && (bucketLocation[bucketOfJ] > 0)) {
                                std::cout << "Particle i=" << particle << " in bucket " << bucket << " collided with j=" << jParticle << " in bucket " << bucketOfJ << std::endl;

                                const double NewMass = mass[particle] + mass[jParticle];

                                v[particle][0] = (v[particle][0] * mass[particle] / NewMass) +
                                                 (v[jParticle][0] * mass[jParticle] / NewMass);
                                v[particle][1] = (v[particle][1] * mass[particle] / NewMass) +
                                                 (v[jParticle][1] * mass[jParticle] / NewMass);
                                v[particle][2] = (v[particle][2] * mass[particle] / NewMass) +
                                                 (v[jParticle][2] * mass[jParticle] / NewMass);

                                mass[particle] = NewMass;
                                std::cout << "New mass=" << mass[particle] << std::endl;
                                std::cout << "Bucket " << bucketOfJ << " has this many particles: " << bucketLocation[bucketOfJ]<< std::endl;
                                std::cout << "J is particle " << j << std::endl;


                                if (j != bucketLocation[bucketOfJ] - 1) {
                                    std::cout << "Reached here" << std::endl;

                                    x[jParticle] = x[buckets[bucketOfJ][bucketLocation[bucketOfJ]]];
                                    v[jParticle] = v[buckets[bucketOfJ][bucketLocation[bucketOfJ]]];
                                    mass[jParticle] = mass[buckets[bucketOfJ][bucketLocation[bucketOfJ]]];

                                    force0[jParticle] = force0[buckets[bucketOfJ][bucketLocation[bucketOfJ]]];
                                    force1[jParticle] = force1[buckets[bucketOfJ][bucketLocation[bucketOfJ]]];
                                    force2[jParticle] = force2[buckets[bucketOfJ][bucketLocation[bucketOfJ]]];

                                    jParticle = buckets[bucketOfJ][bucketLocation[bucketOfJ]];

                                    dist0 = x[jParticle][0] - x[particle][0], dist1 =
                                            x[jParticle][1] - x[particle][1], dist2 =
                                            x[jParticle][2] - x[particle][2];
                                    squareDistance = dist0 * dist0 + dist1 * dist1 + dist2 * dist2;
                                }


                                bucketLocation[bucketOfJ]--;
                                NumberOfBodies--;
                                std::cout << "Particle " << jParticle << " and " << buckets[bucketOfJ][bucketLocation[bucketOfJ]] << std::endl;
                                exit(0);
                                collided = true;
//                                std::cout << "Reached here with " << collided << std::endl;


                            }

                            double distance = std::sqrt(squareDistance);

                            if ((bucketLocation[bucketOfJ] != 0) && (j != bucketLocation[bucketOfJ])) {
                                const double forces = mass[jParticle] * mass[particle] / distance / distance / distance;
                                const double f0 = dist0 * forces;
                                const double f1 = dist1 * forces;
                                const double f2 = dist2 * forces;

                                // x,y,z forces acting on particle i from j
                                force0New[particle] += f0;
                                force1New[particle] += f1;
                                force2New[particle] += f2;

                                // In the final step add the forces on particle j from i
                                if (step == totalSteps - 1) {
                                    force0[jParticle] -= f0;
                                    force1[jParticle] -= f1;
                                    force2[jParticle] -= f2;
                                }

                                minDx = std::min(minDx, distance);
                            }
                        }
                    }
                }

                if (collided == 1){

                    std::cout << "Reached here with particle " << particle << std::endl;
                    std::cout << "Step " << step << " of " << totalSteps << std::endl;

                }


//                std::cout << "Forces now on " << particle << " " << force0New[particle] << " " << force1New[particle] << " "
//                          << force2New[particle] << std::endl;
                x[particle][0] = x[particle][0] + newTimeStepSize * v[particle][0];
                x[particle][1] = x[particle][1] + newTimeStepSize * v[particle][1];
                x[particle][2] = x[particle][2] + newTimeStepSize * v[particle][2];

                v[particle][0] = v[particle][0] + newTimeStepSize * force0New[particle] / mass[particle];
                v[particle][1] = v[particle][1] + newTimeStepSize * force1New[particle] / mass[particle];
                v[particle][2] = v[particle][2] + newTimeStepSize * force2New[particle] / mass[particle];

                maxV = std::max(maxV, std::sqrt(v[particle][0] * v[particle][0] + v[particle][1] * v[particle][1] +
                                                v[particle][2] * v[particle][2]));
                if (collided == 1){

                    std::cout << "Step completed with this many particles left " << NumberOfBodies << std::endl;

                }

            }
        }
    }


    t += timeStepSize;


    delete[]
    buckets;
    delete[]
    force0;
    delete[]
    force1;
    delete[]
    force2;
    delete[]
    force0New;
    delete[]
    force1New;
    delete[]
    force2New;
    delete[]
    bucketLocation;

    if (collided == 1){

        std::cout << "Checking end" << std::endl;


    }
}

// Old Step3 Position: 2.00419158400086 0.994833655104409 0
// New Step3 Position: 2.00104175486776 0.996304306225511 0


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char **argv) {
    if (argc == 1) {
        std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
                  << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting"
                  << std::endl
                  << "  final-time      simulated time (greater 0)" << std::endl
                  << "  dt              time step size (greater 0)" << std::endl
                  << std::endl
                  << "Examples:" << std::endl
                  << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1"
                  << std::endl
                  << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one"
                  << std::endl
                  << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture"
                  << std::endl
                  << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup"
                  << std::endl
                  << std::endl
                  << "In this naive code, only the first body moves" << std::endl;

        return -1;
    } else if ((argc - 4) % 7 != 0) {
        std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)"
                  << std::endl;
        std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
        std::cerr << "run without arguments for usage instruction" << std::endl;
        return -2;
    }

    std::cout << std::setprecision(15);

    setUp(argc, argv);

    openParaviewVideoFile();

    int snapshotCounter = 0;
    if (t > tPlot) {
        printParaviewSnapshot();
        std::cout << "plotted initial setup" << std::endl;
        tPlot = tPlotDelta;
    }

    int timeStepCounter = 0;
    while (t <= tFinal) {
        updateBody();
        timeStepCounter++;
        if (t >= tPlot) {
            printParaviewSnapshot();
            std::cout << "plot next snapshot"
                      << ",\t time step=" << timeStepCounter
                      << ",\t t=" << t
                      << ",\t dt=" << timeStepSize
                      << ",\t v_max=" << maxV
                      << ",\t dx_min=" << minDx
                      << std::endl;

            tPlot += tPlotDelta;
        }
    }

    closeParaviewVideoFile();

    return 0;
}
