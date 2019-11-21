#include <mpi.h>
#include <stdio.h>
int readFile (char*);
int main (int argc, char** argv){
    int size; // кол-во процессов
    int rank; // уникальный идентификатор процессов
    MPI_Status status;

    int N;
    double *U = new double [N-1];
    double *F = new double [N-1];

    char fileName[20] = "file1.txt"; 

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        printf("rank: %d\nsize: %d\n\nEnter mesh count: ", rank, size);
        scanf("%d", &N);
        readFile(fileName);
    }
    //начальные условия (ступенька)
    for (int i = 0; i < N/10; i++) {
		U[i] = 1;
	}
    for (int i = N/10; i < N; i++) {
		U[i] = 0;
	}
    
    

    delete [] U;
    delete [] F;
    U = NULL;
    F = NULL;
    
    MPI_Finalize();
    return 0;
}


void LW (double U[], double F[], double fun, double dt, double dx, int N){
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

int readFile (char* fileName) {
    FILE *file = fopen(fileName, "r");
    struct params {
		char name[20]; 
		double value; 
	};
    struct params option[10];
	int n=0;
    if (file == NULL) {
        printf("File not found!");
    }
    else {
        while (fscanf (file, "%s%lf", option[n].name, &(option[n].value)) != EOF) {
		printf("%s %.2f\n", option[n].name, option[n].value);
		n++;
	}
        fclose(file);
    }
    return 0;
}