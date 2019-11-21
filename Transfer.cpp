#include <mpi.h>
#include <stdio.h>

int main (int argc, char** argv){
    int size; // кол-во процессов
    int rank; // уникальный идентификатор процессов
    MPI_Status status;

    int N;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        printf("rank: %d\nsize: %d\n\nEnter mesh size: ", rank, size);
        scanf("%d", &N);
    }

    
    
    MPI_Finalize();
    return 0;
}


void LV (double U[], double F[], double fun, double dt, double dx, int N){
//Potter p.97                             ����� �����-���������
double *U12 = new double [N-1];
double *F12 = new double [N-1];

for (int i = 0; i <=N-2; i++) {                     // ��������������� ���
	U12[i] = (U[i]+U[i+1])/2 - dt/2/dx*(F[i+1]-F[i]);
	F12[i] = fun*U12[i];
}
U[0] -= dt/dx*(F12[0]);
for (int i = 1; i <=N-2; i++) {                    // �������� ���
	U[i] -= dt/dx*(F12[i]-F12[i-1]);
	F[i] = U[i]*fun;
}
U[N-1] -= dt/dx*(-F12[N-2]);

delete [] U12;
delete [] F12;
U12 = NULL;
F12 = NULL;
}