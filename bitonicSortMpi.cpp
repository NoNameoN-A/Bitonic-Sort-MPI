#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

using namespace std;

int potenzaDiDue(int);
bool eUnaPotenzaDiDue(int);
void inserisco(int*,int);
void riempi(int*,int);
void mergeBitonic(int *,int,int,int,int*);
void compareExchange(int *,int,int,int,int,int,int,int*);
int compareDescending(const void *, const void *);
int compareAscending(const void *, const void *);

//per parte sequenziale
void bitonicSort(int, int, int *);
void bitonicMerge(int,int,int,int *);
void downSwap(int,int,int *);
void upSwap(int,int,int*);

int main(int arg,char* arc[]){
  int dimensioneArray;
	int processore, processoriSize;
	double inizioTempo, fineTempo;

  MPI_Init(&arg,&arc);
	MPI_Comm_size(MPI_COMM_WORLD,&processoriSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&processore);

  if(processore == 0){
    cout << "Inserisci la dimensione, in ogni caso verrà approssimata alla potenza di 2 per eccesso." << endl;
    cin >> dimensioneArray;
    if(dimensioneArray < 2) dimensioneArray = 2;
    dimensioneArray = potenzaDiDue(dimensioneArray);
    cout << "Dimensione dell'Array: " << dimensioneArray << endl;
    inizioTempo = MPI_Wtime();
  }
  MPI_Bcast(&dimensioneArray,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  int* arrayIntero = new int[dimensioneArray];

  	if(processore==0)
	{
		if(processoriSize==1)
		{
			arrayIntero=new int [dimensioneArray];
			riempi(arrayIntero,dimensioneArray);
			bitonicSort(0, dimensioneArray-1,arrayIntero);
			for(int i = 0; i < dimensioneArray; i++)
        cout << "Array[" << i << "] -> " << arrayIntero[i] << endl;
      fineTempo = MPI_Wtime();
      cout << "Tempo impiegato: " << fineTempo-inizioTempo << endl;
		}
	}

  riempi(arrayIntero, dimensioneArray);
  MPI_Barrier(MPI_COMM_WORLD);
  int* arrayTemporaneo;

  mergeBitonic(arrayIntero,dimensioneArray,processoriSize,processore,arrayTemporaneo);
  int* arrayCompleto = NULL; //Inizializzazione
  if(processore == 0)
    arrayCompleto = new int [dimensioneArray * processoriSize];
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(arrayIntero, dimensioneArray, MPI_INT,
              arrayCompleto, dimensioneArray, MPI_INT,
              0, MPI_COMM_WORLD);
  if(processore == 0){
    fineTempo = MPI_Wtime();
    for(int i = 0; i < dimensioneArray; i++)
      cout << "Array[" << i << "] -> " << arrayCompleto[i] << endl;
    cout << "Tempo impiegato: " << fineTempo-inizioTempo << endl;
    delete[] arrayCompleto;
  }
  delete[] arrayIntero;

  MPI_Finalize();
  return 0;
}

int potenzaDiDue(int numero){
	if (eUnaPotenzaDiDue(numero)) return numero;
	if (numero == 0) return 1;
 	int iterazioni = 0;
	while (numero > 0){
		iterazioni++;
    numero /= 2;
  }
  numero = 1;
  for(int i = 0; i < iterazioni; i++) numero *= 2;
  return numero;
}

bool eUnaPotenzaDiDue(int numero){
  while (numero > 0){
    if (numero%2==1 && numero/2!=0) return false;
    numero /= 2;
  }
    return true;
}

void inserisco(int* arrayIntero, int dimensioneArray){
	cout << "Inserisci " << dimensioneArray << " numeri." <<endl;
	for(int i=0; i<dimensioneArray; i++)
		cin >> arrayIntero[i];
	cout << "array riempito correttamente." << endl;
}

void riempi(int* arrayIntero,int dimensioneArray){
	for(int i=0; i<dimensioneArray; i++)
		arrayIntero[i]=rand()%dimensioneArray;
}

void mergeBitonic(int* arrayIntero, int dimensioneArray,int processoriSize,int processore,int* arrayTemporaneo) {
  arrayTemporaneo = new int[dimensioneArray*2];
  int log = processoriSize;
  int pow2i = 2;
  int tag = 0;

  for(int i=1; log > 1 ; i++) {
    int pow2j = pow2i;
    for(int j=i; j >= 1; j--) {
      tag++;
      for(int node=0; node < processoriSize; node += pow2j) {
      	for(int k=0; k < pow2j/2; k++) {
      	  compareExchange(arrayIntero, dimensioneArray, node+k, node+k+pow2j/2, ((node+k) % (pow2i*2) >= pow2i),tag,processore,arrayTemporaneo);
      	}
      }
      pow2j /= 2;
    }
    pow2i *= 2;
    log /= 2;
  }
  delete[]arrayTemporaneo;
}

void compareExchange(int *arrayIntero, int dimensioneArray,
            		     int node1, int node2, int biggerFirst,
            		     int tag,int processore,int* arrayTemporaneo) {
  if ((node1 != processore) && (node2 != processore)) return;
  memcpy(arrayTemporaneo, arrayIntero, dimensioneArray*sizeof(int));	//copia l'array in 'arrayTemporaneo'
  MPI_Status status;
  //ottengo i numeri dall'altro nodo.
  //nodeFrom deve essere diverso dal processore quindi prende l'altro nodo, se node1 è uguale a processore allora nodeFrom = node2 e viceversa
  int nodeFrom = node1==processore ? node2 : node1;
  if (node1 == processore) {
    //Mando dal processo nodeFrom a tutti i processi l'array dalla posizione 0 fino alla metà
    MPI_Send(arrayIntero, dimensioneArray, MPI_INT, nodeFrom, tag, MPI_COMM_WORLD);
    MPI_Recv(&arrayTemporaneo[dimensioneArray], dimensioneArray, MPI_INT, nodeFrom, tag,
	     MPI_COMM_WORLD, &status);
  }
  else {
    //Altrimenti lo mando dalla metà fino alla fine
    MPI_Recv(&arrayTemporaneo[dimensioneArray], dimensioneArray, MPI_INT, nodeFrom, tag,
	     MPI_COMM_WORLD, &status);
    MPI_Send(arrayIntero, dimensioneArray, MPI_INT, nodeFrom, tag, MPI_COMM_WORLD);
  }

  //Ordino l'array
  if (biggerFirst)  qsort(arrayTemporaneo, dimensioneArray*2, sizeof(int), compareDescending);
  else  qsort(arrayTemporaneo, dimensioneArray*2, sizeof(int), compareAscending);

  //A seconda se è il processore è node1 o node2 mi salvo l'array temporaneo dall'inizio o dalla metà
  if (node1 == processore)
    memcpy(arrayIntero, arrayTemporaneo, dimensioneArray*sizeof(int));
  else
    memcpy(arrayIntero, &arrayTemporaneo[dimensioneArray], dimensioneArray*sizeof(int));
}

int compareDescending(const void *item1, const void *item2)
{
  int x = * ( (int *) item1), y = * ( (int *) item2);
  return y-x;
}

int compareAscending(const void *item1, const void *item2){
  int x = * ( (int *) item1), y = * ( (int *) item2);
  return x-y;
}

void quick(int* arrayIntero, int left, int right) {

  int i = left;
  int j = right;
  int temp;
  int pivot = arrayIntero[(left + right) / 2];

  while (i <= j) {
    while (arrayIntero[i] < pivot)
	i++;
    while (arrayIntero[j] > pivot)
	j--;

    if (i <= j) {
      temp = arrayIntero[i];
      arrayIntero[i] = arrayIntero[j];
      arrayIntero[j] = temp;
      i++;
      j--;
    }
  }

  if (left < j){ quick(arrayIntero, left, j); }
  if (i < right){ quick(arrayIntero, i, right); }
}

// qui per la parte sequenziale
void upSwap(int i, int j, int *A)
{
    if (A[j] < A[i])
    {
        int temp = A[j];
        A[j] = A[i];
        A[i] = temp;
    }
}

void downSwap(int i, int j, int *A)
{
    if (A[i] < A[j])
    {
        int temp = A[j];
        A[j] = A[i];
        A[i] = temp;
    }
}
// bitonicMerge(j, i+j-1, dir, A);
void bitonicMerge(int first, int last, int dir, int *A)
{
    //conto elementi di swap
    int cont = 0;
    int size = last - first + 1;


    for (int j=size/2; j>0; j/=2)
    {
        cont = 0;
        for (int i=first; i+j <= last; i++)
        {
            if (cont < j)
            {
                if (dir == 1) // 1
                    upSwap(i, i+j, A);
                else //
                    downSwap(i, i+j, A);

                cont++;
            }
            else
            {
                cont = 0;
                i += j-1;
            }
        }
    }
}

void bitonicSort(int first, int last, int *A)
{
    int size = last - first + 1;
    int dir;
    for (int i = 2; i <= size; i*= 2)
    {
        for (int j = 0; j < size; j = j + i)
        {
            if (((j / i) % 2) == 0)
                dir = 1; //up
            else
                dir = 0; //down

            bitonicMerge(j, i+j-1, dir, A);
        }
    }
}
