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

        long tStart, tEnd, tDiff;
#ifndef __MACH__

        tStart = get_nanos();
#endif
        int processCount;
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
        MPI_Request requests[processCount - 1];
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        char myProcessorName[MPI_MAX_PROCESSOR_NAME];
        int processorNameLen;
        MPI_Get_processor_name(myProcessorName, &processorNameLen);
        unsigned long addUntil;
        unsigned long results[processCount - 1];
        MPI_Barrier(MPI_COMM_WORLD);
        if (myID == 0) {
            addUntil = 100000000l;
            printf("MPI_Master: %d, %s. Add until %lu\n", myID, myProcessorName, addUntil);
            MPI_Bcast((void *) &addUntil, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
            int i;
            for (i = 1; i < processCount; i++) {
                printf("Receiving from %d\n", i);
                MPI_Irecv(&results[i], 1, MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD, &requests[i]);
            }
            for (i = 1; i < processCount; i++) {
                printf("Waiting %d\n", i);
                MPI_Wait(&requests[i], NULL);
            }
#ifndef __MACH__
            tEnd = get_nanos();
            tDiff = tEnd - tStart;
            printf("Master: Total time: %ld ns\n", tDiff);
#endif
            /* verify */
            for (i = 1; i < processCount; i++) {
                if (results[i] != addUntil) {
                    printf("VERIFY: ERR: Node ID %d outputs: %lu\n", i, results[i]);
                } else {
                    printf("VERIFY: Node ID %d is just fine.\n", i);
                }
            }
        } else {
            MPI_Bcast((void *) &addUntil, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
            unsigned long i;
            for (i = 0; i < addUntil; i++);
            MPI_Send(&addUntil, 1, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
            printf("MPI_Slave: FIN %d, %s. Add until %lu\n", myID, myProcessorName, addUntil);
        }
        MPI_Finalize();
    } else if (strcmp(argv[1], "Matrix") == 0) {
        /* program_1 Matrix_2 <size>_3 -G_4 <-1|m>_5 */
        if (argc != 5 || strcmp(argv[3], "-G") != 0) {
            printf("Matrix: Error: Number of arguments. ./program Matrix <size> -G [-1|<m>] (-1=default)");
            MPI_Finalize();
            return EXIT_FAILURE;
        }
        long tStart, tEnd, tDiff;

        int size = atoi(argv[3]), G = atoi(argv[4]);
#ifndef __MACH__
        tStart = get_nanos();
#endif
        int processCount;
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
        MPI_Request requests[processCount - 1];
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        char myProcessorName[MPI_MAX_PROCESSOR_NAME];
        int processorNameLen;
        MPI_Get_processor_name(myProcessorName, &processorNameLen);

        int eachSize = size / processCount;
        int **A = NULL, **B = NULL, **C = NULL, *A_row = NULL, **B_copy;
        if (myID == 0) {

            A = malloc(sizeof(int *) * size);
            B = malloc(sizeof(int *) * size);
            C = malloc(sizeof(int *) * size);
            {
                int i;
                for (i = 0; i < size; i++) {
                    A[i] = malloc(sizeof(int) * size);
                    B[i] = malloc(sizeof(int) * size);
                    int j;
                    for (j = 0; j < size; j++) {
                        A[i][j] = rand();
                        B[i][j] = rand();
                    }
                    C[i] = malloc(sizeof(int) * size); /* don't know C yet!  */
                }
            }
        } else {
            A_row = malloc(sizeof(int) * size);
            B = malloc(sizeof(int *) * size);
            int i;
            for (i = 0; i < size; i++) {
                B[i] = malloc(sizeof(int) * size);
            }
        }

        /* Partition */
        int eachSizes[processCount];
        int displacement[processCount];
        if (G == -1) {
            /* default, use each size */
            int i;
            for (i = 0; i < processCount; i++) {
                if (i == processCount - 1) {
                    eachSizes[i] = size - eachSize * i;
                    displacement[i] = displacement[i - 1] + (size - eachSize * i);
                } else {
                    eachSizes[i] = eachSize;
                    if (i != 0) {
                        displacement[i] = displacement[i - 1] + eachSize;
                    } else
                        displacement[i] = 0;
                }
            }
        } else {
            int fairProcs = size / G;
            int i;
            if (fairProcs >= processCount) {
                /* if G is really small and distributes to a lot of processes, let n-1 nodes be fair
                 * and make the last node handle the rest. */
                for (i = 0; i < processCount; i++) {
                    if (i == processCount - 1) {
                        eachSizes[i] = size - G * i;
                        displacement[i] = displacement[i - 1] + (size - G * i);
                    } else {
                        eachSizes[i] = G;
                        if (i != 0)
                            displacement[i] = displacement[i - 1] + G;
                        else
                            displacement[i] = 0;
                    }
                }
            } else {
                /* if G is really big and only distributes to a small amount of processes, allocate to fairProcs-1
                 * and distribute the tasks fairly to the rest of the nodes. */
                for (i = 0; i < fairProcs - 2; i++) {
                    eachSizes[i] = G;
                    if (i != 0)
                        displacement[i] = displacement[i - 1] + G;
                    else
                        displacement[i] = 0;
                }
                int restSize = (processCount - fairProcs + 2) / (size - ((fairProcs - 2) * G));
                for (++i; i < processCount; i++) {
                    if (i == processCount - 1) {
                        int currentSize =
                                size - (restSize * (processCount - fairProcs + 1)) - (G * (fairProcs - 2));
                        eachSizes[i] = currentSize;
                        displacement[i] = displacement[i - 1] + currentSize;
                    } else {
                        eachSizes[i] = restSize;
                        displacement[i] = displacement[i - 1] + restSize;
                    }
                }

            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (myID == 0)
            MPI_Scatterv(&(A[0][0]), eachSizes, displacement, MPI_INT, A_row, eachSizes[myID], MPI_INT, 0,
                         MPI_COMM_WORLD);
        else
            MPI_Scatterv(NULL, NULL, NULL, MPI_INT, A_row, eachSizes[myID], MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(&(B[0][0]), size * size, MPI_INT, 0, MPI_COMM_WORLD);

        /* Perform work */

        /* Finalize */
        MPI_Finalize();
    } else if (strcmp(argv[1], "Strassen") == 0) {

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
        /* Init */
        int processCount;
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
        MPI_Request requests[processCount - 1];
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        char myProcessorName[MPI_MAX_PROCESSOR_NAME];
        int processorNameLen;
        MPI_Get_processor_name(myProcessorName, &processorNameLen);
        /* Do work */
        printf("Partition\n");
        int SIZE = 100;
        int sizes[processCount];
        int rows[processCount];
        int displ[processCount];
        int displ_rows[processCount];
        int *M = NULL, *N, *M_Recv, *O_Parts, *O;
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
        printf("Allocate\n");
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
        O_Parts = malloc(sizeof(int) * sizes[myID]);
        printf("%d: Size: %d, Displacement: %d.\n", myID, sizes[myID], displ[myID]);
        printf("Send\n");
        /* Distribute */
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatterv((myID == 0) ? (&(M[0])) : NULL, sizes, displ, MPI_INT, &(M_Recv[0]), sizes[myID], MPI_INT, 0,
                     MPI_COMM_WORLD);
        MPI_Bcast(&(N[0]), SIZE * SIZE, MPI_INT, 0, MPI_COMM_WORLD);
        /* Do work! */
        {
            int i, k, j, l=0;
            for (i = 0; i < sizes[myID] / SIZE; i++) /* m */
                for (k = 0; k < SIZE; k++) /* n */
                    for (j = 0; j < SIZE; j++) /* p , m*n n*p */
                        O_Parts[i * SIZE + j] += M_Recv[i * SIZE + k] * N[k * SIZE + j];
        }

        MPI_Gatherv(&(O_Parts[0]), sizes[myID], MPI_INT, &(O[0]), sizes, displ, MPI_INT, 0, MPI_COMM_WORLD);
        /* Finalize */
        printf("Hi! My ID is %d and my processor name is %s. I'm out of %d processes in total.\n", myID,
               myProcessorName, processCount);
        free(M);
        free(N);
        free(M_Recv);
        MPI_Finalize();
    } else {
        printf("No such program as \" %s \".\n", argv[1]);
        MPI_Finalize();
    }
    return 0;
}