#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv)
{
    int world_rank, world_size;
    int row_rank, row_size;
    int col_rank, col_size;
    int n_comm=2;
    int row_color, col_color;
    MPI_Comm row_comm;
    MPI_Comm col_comm;

    MPI_Init(&argc, &argv);   

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
  
    if (argc>1)
        n_comm=atoi(argv[1]);

    if (world_size%n_comm !=0 )
       {
       if ( world_rank == 0)
           {
           fprintf(stderr, " world_size%%n_comm!=0 \n ");
           fprintf(stderr, " Aborting... ");
           }
       MPI_Abort(MPI_COMM_WORLD, 1) ;
       }
    
    if(world_rank==0)
        printf(" Splitting %d process in %d communicators \n", world_size, n_comm);

    row_color = world_rank / n_comm;

    MPI_Comm_split(MPI_COMM_WORLD, row_color, world_rank, &row_comm); 

    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);

    printf(" World rank %d/%d  row rank %d/%d \n", world_rank, world_size, row_rank, row_size);

    col_color = world_rank % n_comm;

    MPI_Comm_split(MPI_COMM_WORLD, col_color, world_rank, &col_comm); 
    
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);
   
    printf(" World rank %d/%d  col rank %d/%d \n", world_rank, world_size, col_rank, col_size);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    
    MPI_Finalize();
    return 0;

}
