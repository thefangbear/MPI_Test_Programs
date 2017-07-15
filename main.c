#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int processCount;
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);
    int myID;
    MPI_Comm_rank(MPI_COMM_WORLD, &myID);
    char myProcessorName[MPI_MAX_PROCESSOR_NAME];
    int processorNameLen;
    MPI_Get_processor_name(myProcessorName,&processorNameLen);
    if (myID == 0) {
        printf("MPI_Program: Hello world! I'm MPI master.\n");
        printf("My ID is %d and my processor name is %s.\n",myID,myProcessorName);
    } else {
        printf("MPI_Program: Hello world: I'm MPI slave.\n");
        printf("My ID is %d and my processor name is %s.\n",myID,myProcessorName);
    }
    MPI_Finalize();
    return 0;
}