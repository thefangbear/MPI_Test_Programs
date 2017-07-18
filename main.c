/*
 * Copyright (C) 2017- Rui-Jie Fang <rfang@temple.edu>.
 * Cite-as: Rui-Jie Fang. MPI Test Programs. GitHub.com.
 * Https://github.com/thefangbear/MPI_Test_Programs. 2017.
 */
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <memory.h>

/*
 * Timing is only available on non-Mac machines
 * (OSX does not have C11 timing functions we need)
 */
#ifndef __MACH__

static long get_nanos(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (long) ts.tv_sec * 1000000000L + ts.tv_nsec;
}

#endif

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    /* determine which program to run */

    if (argc < 2) {
        printf("Err: malformed arguments\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "Addition") == 0) {

    } else if (strcmp(argv[1], "Matrix") == 0) {
        /* Init */
        double tStart, tEnd, tDiff;
        tStart = MPI_Wtime();
        int processCount;
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
        MPI_Request requests[processCount - 1];
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        char myProcessorName[MPI_MAX_PROCESSOR_NAME];
        int processorNameLen;
        MPI_Get_processor_name(myProcessorName, &processorNameLen);
        /* Do work */
        int SIZE = 100;
        int sizes[processCount];
        int rows[processCount];
        int displ[processCount];
        int displ_rows[processCount];
        int *M = NULL, *N, *M_Recv, *O_Parts, *O = NULL;
        {
            /* Partition */
            int avg = SIZE / processCount;
            int i;
            sizes[0] = avg * SIZE;
            rows[0] = avg;
            displ[0] = 0;
            displ_rows[0] = 0;
            if (processCount > 1)
                for (i = 1; i < processCount; i++) {
                    if (i != processCount - 1) {
                        sizes[i] = avg * SIZE;
                        rows[i] = avg;
                        displ[i] = displ[i - 1] + avg * SIZE;
                    } else {
                        int remainingRows = SIZE - i * avg;
                        sizes[i] = remainingRows * SIZE;
                        rows[i] = avg;
                        displ[i] = displ[i - 1] + avg * SIZE;
                    }
                }
        }
        /* Allocate: we store arrays in row-major order for the ease of scattering */
        N = malloc(sizeof(int) * SIZE * SIZE);
        if (myID == 0) {
            M = malloc(sizeof(int) * SIZE * SIZE);
            O = malloc(sizeof(int) * SIZE * SIZE);
            int i;
            for (i = 0; i < SIZE * SIZE; i++) {
                M[i] = rand();
                N[i] = rand();
            }
        }
        M_Recv = malloc(sizeof(int) * sizes[myID]);
        O_Parts = malloc(sizeof(int) * (sizes[myID] + 1));
        printf("%d: Size: %d, Displacement: %d.\n", myID, sizes[myID], displ[myID]);
        /* Distribute */
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatterv((myID == 0) ? (&(M[0])) : NULL, sizes, displ, MPI_INT, &(M_Recv[0]), sizes[myID], MPI_INT, 0,
                     MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&(N[0]), SIZE * SIZE, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        /* Do work! */
        {
            int i, k, j, l = 0;
            for (i = 0; i < sizes[myID] / SIZE; i++) {/* m */
                for (k = 0; k < SIZE; k++) {/* n */
                    for (j = 0; j < SIZE; j++) { /* p , m*n n*p */
                        O_Parts[l] += M_Recv[i * SIZE + k] * N[k * SIZE + j];
                        /*              O_Parts[i * SIZE + j] += M_Recv[i * SIZE + k] * N[k * SIZE + j]; */
                    }
                    l++;
                }
            }
        }
        MPI_Gatherv(&(O_Parts[0]), sizes[myID], MPI_INT, &(O[0]), sizes, displ, MPI_INT, 0, MPI_COMM_WORLD);
        /* Finalize */
        tEnd = MPI_Wtime();
        tDiff = tEnd - tStart;
        if (myID == 0)
            printf("Master: Costs %f sec.\n", tDiff);
        free(M);
        free(N);
        free(M_Recv);
        MPI_Finalize();
    } else if (strcmp(argv[1], "Quicksort") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "Mergesort") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "TSP") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "DFS") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "Search") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "PI") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "Test") == 0) {
        int processCount;
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
        MPI_Request requests[processCount - 1];
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        char myProcessorName[MPI_MAX_PROCESSOR_NAME];
        int processorNameLen;
        MPI_Get_processor_name(myProcessorName, &processorNameLen);
        printf("Hi! My ID is %d and my processor name is %s. I'm out of %d processes in total.\n", myID,
               myProcessorName, processCount);
        MPI_Finalize();
    } else if (strcmp(argv[1], "TestConcept") == 0) {
        printf("TestConcept\n");

        MPI_Finalize();
    } else {
        printf("No such program as \" %s \".\n", argv[1]);
        MPI_Finalize();
    }
    return 0;
}