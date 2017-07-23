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

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    /* determine which program to run */

    if (argc < 2) {
        printf("Err: malformed arguments\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "Addition") == 0) {
        printf("Addition\n");
        /* Init */
        if (argc < 3) {
            printf("Err: Incorrect # of args.\n");
            MPI_Finalize();
            return 0;
        }
        unsigned long total;

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
        if (myID == 0) {
            total = strtoul(argv[2], NULL, 10);
            printf("Total: %lu", total);
        }
        unsigned long results[processCount];
        MPI_Barrier(MPI_COMM_WORLD);
        /* Partition */
        MPI_Bcast(&total, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        /* Do work */
        unsigned long result;
        for (result = 0; result < total; result++);
        MPI_Gather(&result, 1, MPI_UNSIGNED_LONG, &results[0], processCount, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        if (myID == 0) {
            int i;
            unsigned long adder = 0;
            for (i = 0; i < processCount; i++)
                adder += results[i];
        }
        tEnd = MPI_Wtime();
        tDiff = tEnd - tStart;
        if (myID == 0)
            printf("Master: Addition finished in %f sec.\n", tDiff);
        MPI_Finalize();
    } else if (strcmp(argv[1], "Matrix") == 0) {
        /* Init */
        double ttStart, ttEnd, ttDiff, tInitStart, tInitEnd, tCompStart, tCompEnd, tAllocStart, tAllocEnd, tAllocDiff, tInitDiff, tCompDiff;
        ttStart = MPI_Wtime();
        int processCount;
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
        MPI_Request requests[processCount - 1];
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        char myProcessorName[MPI_MAX_PROCESSOR_NAME];
        int processorNameLen;
        MPI_Get_processor_name(myProcessorName, &processorNameLen);
        ttEnd = MPI_Wtime();
        /* Do work */
        int SIZE = 100;
        int sizes[processCount];
        int rows[processCount];
        int displ[processCount];
        int displ_rows[processCount];
        int G = -1;
        if (argc == 4 && strcmp(argv[2], "-G") == 0) {
            G = (int) strtol(argv[3], NULL, 10);
            printf("-G: %d\n", G);
        } else if (argc == 4 && strcmp(argv[2], "--size") == 0) {
            SIZE = (int) strtol(argv[3], NULL, 10);
            printf("--size: %d\n", SIZE);
        } else if (argc == 6
                   && strcmp(argv[2], "-G") == 0
                   && strcmp(argv[4], "--size") == 0) {
            G = (int) strtol(argv[3], NULL, 10);
            SIZE = (int) strtol(argv[5], NULL, 10);
            printf("-G: %d\n", G);
            printf("--size: %d\n", SIZE);
        } else if (argc == 6
                   && strcmp(argv[2], "--size") == 0
                   && strcmp(argv[4], "-G") == 0) {
            G = (int) strtol(argv[5], NULL, 10);
            SIZE = (int) strtol(argv[3], NULL, 10);
            printf("-G: %d\n", G);
            printf("--size: %d\n", SIZE);
        }
        tInitStart = MPI_Wtime();
        int *M = NULL, *N, *M_Recv, *O_Parts, *O = NULL;
        {
            /* Not completely fair partition */
            if (G == -1) {
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
            } else {
                if (G > SIZE) {
                    printf("Err: G is too large.\n");
                    MPI_Finalize();
                    return EXIT_FAILURE;
                }
                int numHosts = SIZE / G;
                int i;
                sizes[0] = G;
                displ[0] = 0;
                if (numHosts > 1)
                    for (i = 0; i < numHosts; i++) {
                        if (i != numHosts - 1) {
                            sizes[i] = G;
                            displ[i] = displ[i - 1] + G;
                        } else {
                            sizes[i] = SIZE - G * i;
                            displ[i] = displ[i - 1] + G;
                        }
                    }
                for (++i; i < processCount; i++) {
                    sizes[i] = 0;
                    displ[i] = displ[i - 1];
                }
            }
        }
        tInitEnd = MPI_Wtime();
        tAllocStart = MPI_Wtime();
        /* Allocate: we store arrays in row-major order for the ease of scattering */
        N = malloc(sizeof(int) * SIZE * SIZE);
        if (myID == 0) {
            M = malloc(sizeof(int) * SIZE * SIZE);
            O = malloc(sizeof(int) * SIZE * SIZE);
            int i, j;
            for (i = 0; i < SIZE; i++) {
                for (j = 0; j < SIZE; j++) {
                    M[i * SIZE + j] = i * j;
                    N[i * SIZE + j] = i * j;
                }
            }
        } else {
            int i;
            int j;
            for (i = 0; i < SIZE; i++) {
                for (j = 0; j < SIZE; j++)
                    N[i * SIZE + j] = i * j;
            }
        }
        M_Recv = malloc(sizeof(int) * sizes[myID]);
        O_Parts = malloc(sizeof(int) * (sizes[myID] + 1));
        printf("%d: Size: %d, Displacement: %d.\n", myID, sizes[myID], displ[myID]);
        tAllocEnd = MPI_Wtime();
        /* Distribute */
        MPI_Barrier(MPI_COMM_WORLD);
        tCompStart = MPI_Wtime();

        MPI_Scatterv((myID == 0) ? (&(M[0])) : NULL, sizes, displ, MPI_INT, &(M_Recv[0]), sizes[myID], MPI_INT, 0,
                     MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        //   MPI_Bcast(&(N[0]), SIZE * SIZE, MPI_INT, 0, MPI_COMM_WORLD);
        //   MPI_Barrier(MPI_COMM_WORLD);
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
        tCompEnd = MPI_Wtime();
        ttDiff = ttEnd - ttStart;
        tInitDiff = tInitEnd - tInitStart;
        tAllocDiff = tAllocEnd - tAllocStart;
        tCompDiff = tCompEnd - tCompStart;
        if (myID == 0) {
            printf("Master: T Cost: %f sec. Init cost: %f sec. Alloc Cost: %f sec. Comp cost: %f sec. Total cost: %f sec.\n",
                   ttDiff, tInitDiff, tAllocDiff, tCompDiff, ttDiff + tInitDiff + tAllocDiff + tCompDiff);
            free(M);
        }

        free(N);
        free(M_Recv);
        free(O);
        free(O_Parts);

        MPI_Finalize();
    } else if (strcmp(argv[1], "Quicksort") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "Mergesort") == 0) {

        MPI_Finalize();
    } else if (strcmp(argv[1], "Sort") == 0) {
        /* Init */
        printf("Sort\n");
        if (argc != 4 || strcmp(argv[3], "--size") != 0) {
            printf("Err: Num of arguments.\n");
            MPI_Finalize();
            return EXIT_FAILURE;
        }
        int size = (int) strtol(argv[3], NULL, 10);
        printf("Size is: %lu\n", size);
        int processCount;
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
        MPI_Request requests[processCount - 1];
        int myID;
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
        char myProcessorName[MPI_MAX_PROCESSOR_NAME];
        int processorNameLen;
        MPI_Get_processor_name(myProcessorName, &processorNameLen);
        int *A = NULL, *A_Recv = NULL;
        if (myID == 0)
            A = malloc(sizeof(int) * size);
        /* Partition */
        int eachSizes[processCount];
        int displ[processCount];
        {
            int i, avg = size / processCount;
            eachSizes[0] = avg;
            displ[0] = 0;
            if (processCount > 1)
                for (i = 1; i < processCount; i++) {
                    if (i != processCount - 1) {
                        eachSizes[i] = avg;
                        displ[i] = displ[i - 1] + avg;
                    } else {
                        int remaining = size - avg * i;
                        eachSizes[i] = remaining;
                        displ[i] = displ[i - 1] + avg;
                    }
                }
        }
        /* Allocate */
        if (myID == 0) {
            int i;
            for (i = 0; i < size; i++)
                A[i] = rand();
        }
        A_Recv = malloc(sizeof(int) * eachSizes[myID]);
        /* Do work! */
        MPI_Scatterv(&A, eachSizes, displ, MPI_INT, &A_Recv, eachSizes[myID], MPI_INT, 0, MPI_COMM_WORLD);
        {
            int i;
            for (i = 1; i < eachSizes[myID]; i++) {
                int j = i;
                while (A_Recv[j] > A_Recv[j - 1] && j > 0) {
                    int temp = A_Recv[j - 1];
                    A_Recv[j - 1] = A_Recv[j];
                    A_Recv[j] = temp;
                    j--;
                }
            }
        }
        MPI_Gatherv(&A_Recv, eachSizes[myID], MPI_INT, &A, eachSizes, displ, MPI_INT, 0, MPI_COMM_WORLD);
        // TODO merge
        {
            int i;
        }
        /* Cleanup */
        if (myID == 0)
            free(A);
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
