#include <mpi.h>
#include <stdio.h>

#define DEBUG

int readFile (char*, double[]);
int consoleGraph(double[], int);
int writeFile (double[], int, char*);
void getFlow (double U[], double F[], int n, double velocity);
void LW (double U[], double F[], double fun, double dt, double dx, int N);
int boundExchange(double &left, double &right, int size, int rank, MPI_Request request, MPI_Status status);

int main (int argc, char** argv){
    char fileName[20] = "file1.txt", add[5] = "a", create[5] = "w"; 
    int N, n, dt = 0;  
    double* Uglobal;
    double velocity = 1;

    int size; // кол-во процессов
    int rank; // уникальный идентификатор процессов
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        double options[5];
        readFile(fileName, options); 
        N = (int) options[0];           // размер сетки
        while (N % size != 0) N++; 
        // Nlocal = (N / size) + 1;    //количество ячеек и размер массива !!
        // N=Nlocal*size;              //выравнивание
        N = 20;
        printf("Parameters are readed. Total mesh count is %d\n", N);
        Uglobal = new double[N];
        for (int i = 0; i < N/10; i++) {    // начальные условия (ступенька)
            Uglobal[i] = 1;
        }
        for (int i = N/10; i < N; i++) {
            Uglobal[i] = 0;
        }
        printf("Initial condition: \n");
        consoleGraph(Uglobal, N);             // график в консоль
        writeFile(Uglobal, N, create);        // вывод в новый файл
        n = N / size;
    }
    MPI_Bcast((void*)&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double* U = new double[n+2];
    double* F = new double[n+2];

    MPI_Scatter((void*)Uglobal, n, MPI_DOUBLE, (void*)&U[1], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    U[0] = U[1];
    U[n] = 2;
    U[n+1] = U[n];
    // if (rank == 0) MPI_Send((void*)&U[n+1], 1, MPI_DOUBLE, 1, 123, MPI_COMM_WORLD);  
    // if (rank == 1) MPI_Recv((void*)&U[0], 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, &status);  
    boundExchange(U[0], U[n+1], size, rank, request, status);
    getFlow(U, F, n, velocity);
    printf("left bound = %f\tright bound = %f\tcenter = %f\trank = %d\n", U[0], U[n+1], U[n/2], rank);
    
    // вычисление нового временного слоя
    while (dt < 10)
    {
        LW(U, F, velocity, 1, 1, n);
        dt++;
        // U[n] = U[n-1];
        boundExchange(U[0], U[n], size, rank, request, status);
        if (rank == 0) {
            U[0] = U[1];
            printf("U[n] = %f, U[n-1] = %f h = %d\n", U[n], U[n-1], dt);
            // printf("F[n] = %f, F[n-1] = %f h = %d\n", F[n], F[n-1], dt);
        } 
        else {
            U[1] = U[0];
            if (rank == 1)
            printf("rank %d, left boundary %f, first mesh = %f\n", rank, U[0], U[1]);
        }
        if (dt % 1 == 0){
            // сборка данных (gather)
            MPI_Gather((void*)U, n, MPI_DOUBLE, (void*)Uglobal, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if (rank == 0) {
                writeFile(Uglobal, N, add); 
                // printf("left bound = %f\tright bound = %f\tcenter = %f\trank = %d\t h = %d\n", U[0], U[n+1], U[n/2], rank, dt);
            }        
            // if (rank == 1) printf("left bound = %f\tright bound = %f\tcenter = %f\trank = %d\t h = %d\n", U[0], U[n+1], U[n/2], rank, dt);
        }
    }
    // swap (не использлуется в трёх-слойной схеме)

    // вывод в файл -> x; numerical solution; analytic solution
    // printf("rank = %d, pointer to F before del = \t%p\n", rank, F);
    delete[] U; //U = NULL;
    // printf("delete U - OK! rank = %d\n", rank);
    delete[] F; //F = NULL;
    // printf("delete F - OK! rank = %d\n", rank);
    MPI_Finalize();

    return 0;
}

int consoleGraph (double Y[], int x) {
    // printf("ConsoleGraph x = %d", x);
    for (int i = 0; i < x; i++)
    {
        if (Y[i]) {
            printf("|");
        }
        else
        {
            printf(".");
        }
        i++;
    }
    printf("\n");
    return 0;
}

void LW (double U[], double F[], double fun, double dt, double dx, int N){
    //Potter p.97                       
    double* U12 = new double [N];
    double* F12 = new double [N];

    for (int i = 0; i < N; i++) {                     // вспомогательный шаг
    	U12[i] = (U[i] + U[i+1])/2 - dt/2/dx*(F[i+1] - F[i]);
    	F12[i] = fun*U12[i];
    }
    // U[1] -= dt/dx*(F12[1]);
    for (int i = 0; i < N; i++) {                     // основной шаг
    	U[i] -= dt/dx*(F12[i]-F12[i-1]);
    	F[i] = U[i]*fun;
    }
    // U[N] -= dt/dx*(-F12[N-1]);

    delete [] U12; U12 = NULL;
    delete [] F12; F12 = NULL;
}

int readFile (char* fileName, double initArray[]) {
    FILE* file = fopen(fileName, "r");
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
		    initArray[n] = option[n].value;
		    n++;
	    }
        fclose(file);
    }
    return 0;
}

int boundExchange(double &left, double &right, int size, int rank, MPI_Request request, MPI_Status status) {
    for (int i = 0; i < size-1; i++)
    {
        if (rank == i) {
            MPI_Send((void*)&right, 1, MPI_DOUBLE, i+1, 123, MPI_COMM_WORLD);
            // MPI_Isend((void*)&right, 1, MPI_DOUBLE, i+1, 123, MPI_COMM_WORLD, &request);
            // if (rank == 0) printf("I send %f to %d\n", right, i+1);
        }
    }
    for (int i = 1; i < size; i++)
    {
        if (rank == i) {
            MPI_Recv((void*)&left, 1, MPI_DOUBLE, i-1, 123, MPI_COMM_WORLD, &status);  
            // printf("rank %d recieved left bound as %f\n", rank, left);
        }
    }
    return 0;
}

int writeFile (double A[], int N, char* option) {
    FILE* file = fopen("output.txt", option);
    for (int i = 0; i < N; i++)
    {
        fprintf(file, "%d\t%f\n", i, A[i]);
    }
    fprintf(file, "\n");
    fclose(file);
    return 0;
}

void getFlow (double U[], double F[], int n, double velocity) {
    for (int i = 0; i < n+1; i++)             
    {
        F[i] = U[i] * velocity;
    }
}
    