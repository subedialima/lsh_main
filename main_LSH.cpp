#include <iostream>
#include <vector>
#include <string>
#include "parse_geodata.h"
#include "query.h"
#include <geos_c.h>
#include <cfloat>
#include <mpi.h>
#include <unordered_set>

using namespace std;

// Function to calculate the union of local MBRs for a process
GEOSGeometry *calculateAggregateLocalMBR(GEOSContextHandle_t ctx, const vector<GEOSGeometry *> &localMBRs)
{
    // Combine individual MBRs to get the local MBR for the process
    GEOSGeometry *processLocalMBR = nullptr;

    for (const auto &mbr : localMBRs)
    {
        if (!processLocalMBR)
        {
            processLocalMBR = GEOSGeom_clone_r(ctx, mbr);
        }
        else
        {
            GEOSGeometry *unionResult = GEOSUnion_r(ctx, processLocalMBR, mbr);
            GEOSGeom_destroy_r(ctx, processLocalMBR);
            processLocalMBR = unionResult;
        }
    }

    return processLocalMBR;
}

// Return the MBR for a single geometry
GEOSGeometry *minimumBoundingRectangle(GEOSContextHandle_t ctx, GEOSGeometry *geometry)
{
    GEOSGeometry *curMBR = GEOSEnvelope_r(ctx, geometry);
    return curMBR;
}

void printCoordinates(GEOSContextHandle_t ctx, const GEOSGeometry *geometry)
{
    int geometryType = GEOSGeomTypeId_r(ctx, geometry);

    if (geometryType != GEOS_POLYGON)
    {
        cerr << "Invalid geometry type. Expected Polygon." << endl;
        return;
    }

    // Get the exterior ring of the polygon
    const GEOSGeometry *exteriorRing = GEOSGetExteriorRing_r(ctx, geometry);
    if (exteriorRing == nullptr)
    {
        cerr << "Error getting exterior ring." << endl;
        return;
    }

    const GEOSCoordSequence *coordSeq = GEOSGeom_getCoordSeq_r(ctx, exteriorRing);
    if (coordSeq == nullptr)
    {
        cerr << "Error getting coordinate sequence." << endl;
        return;
    }

    unsigned int numPoints;
    GEOSCoordSeq_getSize_r(ctx, coordSeq, &numPoints);

    cout << "Number of Points: " << numPoints << endl;

    for (unsigned int i = 0; i < numPoints; ++i)
    {
        double x, y;
        GEOSCoordSeq_getX_r(ctx, coordSeq, i, &x);
        GEOSCoordSeq_getY_r(ctx, coordSeq, i, &y);

        cout << "Point " << i + 1 << ": (" << x << ", " << y << ")" << endl;
    }

    cout << endl;
}



void processGeometryParallel(const string &filename, int rank, int size, GEOSContextHandle_t ctx)
{
    // Read polygons from the file
    vector<GEOSGeometry *> *geos = read_wkt(filename, ctx);

    if (!geos)
    {
        cerr << "Error reading file: " << filename << endl;
        return;
    }

    // Broadcast the number of geometries to all processes
    int numGeometries = geos->size();
    MPI_Bcast(&numGeometries, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate the range of geometries to process for each process
    int chunkSize = numGeometries / size;
    int startIdx = rank * chunkSize;
    int endIdx = (rank == size - 1) ? numGeometries : (rank + 1) * chunkSize;

    cout << "Process " << rank << " will process geometries from index " << startIdx << " to " << endIdx - 1 << endl;

    // Declare variables outside the loop
    GEOSGeometry *aggregateLocalMBR = nullptr;
    vector<GEOSGeometry *> localMBRs;

    // Calculate MBRs for each polygon in the specified range
    for (int i = startIdx; i < endIdx; ++i)
    {
        GEOSGeometry *currentPolygon = geos->at(i);
        if (!currentPolygon)
        {
            cerr << "Error: Null geometry for polygon " << i << " in process " << rank << endl;
            continue;
        }

        GEOSGeometry *centeredPolygon = centerGeometry(ctx, currentPolygon);
        if (!centeredPolygon)
        {
            cerr << "Error: Null geometry after centering for polygon " << i << " in process " << rank << endl;
            continue;
        }

        GEOSGeometry *mbr = minimumBoundingRectangle(ctx, centeredPolygon);
        if (!mbr)
        {
            cerr << "Error: Null MBR for polygon " << i << " in process " << rank << endl;
            continue;
        }

        localMBRs.push_back(mbr); // Store the local MBR in the vector
        
        // Print coordinates for the MBR of the first polygon processed by each process
        if (i == startIdx)
        {
            cout << "Coordinates for MBR of the first polygon processed by process " << rank << ":" << endl;
            printCoordinates(ctx, mbr);
        }


        // Free memory for the MBR
        GEOSGeom_destroy_r(ctx, centeredPolygon); // Free centeredPolygon instead of mbr
    }

    // Calculate aggregate local MBR for each process
    aggregateLocalMBR = calculateAggregateLocalMBR(ctx, localMBRs);
    
    // Print coordinates for the aggregateLocalMBR of each process
    cout << "Aggregate MBR for process " << rank << ":" << endl;
    printCoordinates(ctx, aggregateLocalMBR);
    
      // MPI Barrier: Ensure all processes have completed reading before proceeding
    MPI_Barrier(MPI_COMM_WORLD);
    


    // Gather local aggregate MBRs from all processes to the root
    vector<GEOSGeometry *> rootAggregateMBRs(size, nullptr);
    MPI_Gather(&aggregateLocalMBR, 1, MPI_LONG_LONG, rootAggregateMBRs.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
   
   // Gather local aggregate MBRs from all processes to the root
//vector<GEOSGeometry *> rootAggregateMBRs;
//rootAggregateMBRs.resize(size, nullptr);

//MPI_Gather(&aggregateLocalMBR, 1, MPI_LONG_LONG, rootAggregateMBRs.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    
   // vector<vector<GEOSGeometry *>> rootAggregateMBRs(size);
   // MPI_Gather(localMBRs.data(), localMBRs.size(), MPI_LONG_LONG, rootAggregateMBRs[0].data(), localMBRs.size(), MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    
    // Gather local aggregate MBRs from all processes to the root
// Gather local aggregate MBRs from all processes to the root


// Output to check the gathered pointers (on each process)
for (int i = 0; i < size; ++i)
{
    if (rank == i)
    {
        cout << "Process " << rank << " gathered pointers: ";
        for (const auto &ptr : rootAggregateMBRs)
        {
            cout << ptr << " ";
        }
        cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


    //MPI_Gather(&aggregateLocalMBR, 1, MPI_BYTE, rootAggregateMBRs.data(), 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    
   

     if (rank == 0)
    {
    cout << "error 1 processes:" << endl;
        // Filter out null pointers
        vector<GEOSGeometry *> filteredRootAggregateMBRs;
        for (const auto &rootMBR : rootAggregateMBRs)
        {  cout << "error 2 processes:" << endl;
            if (rootMBR)
                filteredRootAggregateMBRs.push_back(rootMBR);
                 cout << "error 3 processes:" << endl;
        }
        
        cout << "error 4 processes:" << endl;
        GEOSGeometry *finalAggregateMBR = calculateAggregateLocalMBR(ctx, filteredRootAggregateMBRs); //error is here needs to be resolved
         cout << "error 5 processes:" << endl; 
        cout << "Final Aggregate MBR for all processes:" << endl;
        printCoordinates(ctx, finalAggregateMBR);
cout << "error 6 processes:" << endl;
        // Free memory for the final aggregate MBR
        GEOSGeom_destroy_r(ctx, finalAggregateMBR);
    }

    // Freeing memory
    delete geos;
}

int main(int argc, char **argv)
{
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3)
    {
        if (rank == 0)
        {
            cerr << "Incorrect number of arguments!\n";
            cerr << "Usage: mpirun -n <num_processes> ./parallel_main <dataFile> <inputFile>\n";
        }
        MPI_Finalize();
        return 1;
    }

    string datafile = string(argv[1]);

    // Create new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    // Process geometry in parallel
    processGeometryParallel(datafile, rank, size, ctx);

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before finalizing

    GEOS_finish_r(ctx);
    MPI_Finalize();

    return 0;
}