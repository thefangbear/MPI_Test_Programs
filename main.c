#include <mpi.h>
#include <stdio.h>
#include <time.h>

static long get_nanos(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (long) ts.tv_sec * 1000000000L + ts.tv_nsec;
}

int main(int argc, char **argv) {
    long tStart, tEnd, tDiff;
    tStart = get_nanos();
    MPI_Init(&argc, &argv);
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
        tEnd = get_nanos();
        tDiff = tEnd - tStart;
        printf("Master: Total time: %ld ns\n", tDiff);
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
    return 0;
}