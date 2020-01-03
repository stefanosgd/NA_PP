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
        if (velocity == maxV) {
            buckets[totalBuckets - 1][bucketLocation[totalBuckets - 1]] = i;
            bucketLocation[totalBuckets - 1] += 1;
        } else {
            int bucketChosen = std::floor(velocity / vBucket);
            buckets[bucketChosen][bucketLocation[bucketChosen]] = i;
            bucketLocation[bucketChosen] += 1;
        }
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

    for (int bucket = 0; bucket < totalBuckets; bucket++) {
        int totalSteps = std::pow(2, bucket);
        double newTimeStepSize = timeStepSize / totalSteps;
        for (int step = 0; step < totalSteps; step++) {
            for (int i = 0; i < bucketLocation[bucket]; i++) {
                int particle = buckets[bucket][i];
                force0[particle] = 0.0;
                force1[particle] = 0.0;
                force2[particle] = 0.0;
                for (int j = 0; j < NumberOfBodies; j++) {
                    if (particle != j) {
                        double dist0 = x[j][0] - x[particle][0], dist1 = x[j][1] - x[particle][1], dist2 =
                                x[j][2] - x[particle][2];
                        double squareDistance = dist0 * dist0 + dist1 * dist1 + dist2 * dist2;
                        double distance = std::sqrt(squareDistance);

                        while ((squareDistance <= (0.01 * 0.01)) && (NumberOfBodies > 1)) {

                            const double NewMass = mass[particle] + mass[j];

                            v[particle][0] =
                                    (v[particle][0] * mass[particle] / NewMass) + (v[j][0] * mass[j] / NewMass);
                            v[particle][1] =
                                    (v[particle][1] * mass[particle] / NewMass) + (v[j][1] * mass[j] / NewMass);
                            v[particle][2] =
                                    (v[particle][2] * mass[particle] / NewMass) + (v[j][2] * mass[j] / NewMass);
                            mass[particle] = NewMass;

                            if (j != NumberOfBodies - 1) {
                                x[j] = x[NumberOfBodies - 1];
                                v[j] = v[NumberOfBodies - 1];
                                mass[j] = mass[NumberOfBodies - 1];

                                dist0 = x[j][0] - x[particle][0], dist1 = x[j][1] - x[particle][1], dist2 =
                                        x[j][2] - x[particle][2];
                                squareDistance = dist0 * dist0 + dist1 * dist1 + dist2 * dist2;
                                distance = std::sqrt(squareDistance);
                            }
                            for (int k = 0; k < totalBuckets; k++) {
                                for (int kLoc = 0; kLoc < bucketLocation[k]; kLoc++) {
                                    if (buckets[k][kLoc] == j) {
                                        buckets[k][kLoc] = buckets[k][bucketLocation[k]];
                                        bucketLocation[k]--;
                                    }
                                }
                            }
                            NumberOfBodies -= 1;
                        }

                        if ((NumberOfBodies != 1) && (j != NumberOfBodies)) {
                            const double forces = mass[j] * mass[particle] / distance / distance / distance;

                            // x,y,z forces acting on particle i from j
                            force0[particle] += dist0 * forces;
                            force1[particle] += dist1 * forces;
                            force2[particle] += dist2 * forces;

                            minDx = std::min(minDx, distance);
                        }
                    }
                }
                x[particle][0] = x[particle][0] + newTimeStepSize * v[particle][0];
                x[particle][1] = x[particle][1] + newTimeStepSize * v[particle][1];
                x[particle][2] = x[particle][2] + newTimeStepSize * v[particle][2];

                v[particle][0] = v[particle][0] + newTimeStepSize * force0[particle] / mass[particle];
                v[particle][1] = v[particle][1] + newTimeStepSize * force1[particle] / mass[particle];
                v[particle][2] = v[particle][2] + newTimeStepSize * force2[particle] / mass[particle];

                maxV = std::max(maxV, std::sqrt(v[particle][0] * v[particle][0] + v[particle][1] * v[particle][1] +
                                                v[particle][2] * v[particle][2]));
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
    bucketLocation;

}


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
