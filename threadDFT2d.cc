// Threaded two-dimensional Discrete FFT transform
// YOUR NAME HERE
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.

Complex* Data;
int Width;
int Height;
bool inv = false;
const int N = 1024;
const int numThreads = 16;
Complex* W = new Complex[N / 2];

bool* chainLock;
bool globalLock;
pthread_mutex_t Mutex;
int counter;


void Transpose(Complex* matrix, Complex* buffer)
{
  //Calculating the transpose.
  int pos = 0;
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      buffer[pos] = matrix[i + j * N];
      pos++;
    }
  } 

  for(int i = 0; i < N * N; i++)
  {
    matrix[i] = buffer[i];
  }    
}  

unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

void BitReverseOrder(Complex* matrix, Complex* buffer)
{
  for(int n = 0; n < N; n++) {
    int startPos = n * N;
    Complex* tmp = new Complex[N];
    for(int i = 0; i < N; i++) 
      tmp[i] = matrix[i + startPos];
    for(int i = startPos; i < startPos + N; i++)
      buffer[i] = tmp[ ReverseBits(i) ];
    delete[] tmp;
  }

  for(int i = 0; i < N * N; i++)
    matrix[i] = buffer[i];

}

void prepare_W()
{
  for(int i = 0; i < N / 2; i++)
    W[i] = Complex(cos(2*M_PI*i/N), -sin(2*M_PI*i/N));
}

void prepare_W_inv()
{
  for(int i = 0; i < N / 2; i++)
    W[i] = Complex(cos(2*M_PI*i/N), sin(2*M_PI*i/N));
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{
  counter = numThreads + 1;
  chainLock = new bool[numThreads + 1]; // +1 for main()
  pthread_mutex_init(&Mutex, NULL);
  for(int i = 0; i < counter; i++)
    chainLock[i] = true;
  globalLock = true;
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier(int id) // Again likely need parameters
{
  chainLock[id] = !chainLock[id];
  pthread_mutex_lock(&Mutex);
  int me = counter;
  counter--;
  pthread_mutex_unlock(&Mutex);
  if(me == 1) {
    globalLock = false;
  }
  else {
    while(globalLock != chainLock[id]) {}; // spin
  }
}
                    
void Transform1D(Complex* h, int N, Complex* W)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  Complex tmp;
  for(int np = 2; np <= N; np = np*2)
  {
    for(int i = 0; i < N; i = i + np)
    {
      for(int j = 0; j < np/2; j++)
      {
        int offset = np/2;
        tmp = h[i + j];
        h[i + j] = h[i + j] + W[j*N/np] * h[i + j + offset];
        h[i + j + offset] = tmp - W[j*N/np] * h[i + j + offset];
      }
    }
  }

  if(inv == true)
  {
    for(int i = 0; i < N; i++)
    {
      h[i].real = h[i].real / N;
      h[i].imag = h[i].imag / N;
    }
  }

}

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  unsigned long id = (unsigned long)v;

  int startRow = id * N / numThreads;
  for(int i = 0; i < N / numThreads; i++)
  {
    int startPos = N * (startRow + i);
    Transform1D(Data + startPos, N, W);
  }

  MyBarrier(id);

  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.

  // Create the helper object for reading the image
  InputImage image(inputFN);

  // Create the global pointer to the image array data  
  Data = image.GetImageData(); 
  Width = image.GetWidth();
  Height = image.GetHeight();
  Complex* buffer = new Complex[N * N];

  /****************************    2D FFT    ********************************/

  inv = false;

  /***** row fft *****/

  prepare_W();

  BitReverseOrder(Data, buffer);

  MyBarrier_Init();

  // Create 16 threads
  pthread_t threads[numThreads];
  for(int i = 0; i < numThreads; i++)
    pthread_create(&threads[i], 0, Transform2DTHread, (void*)i);

  // Wait for all threads complete
  MyBarrier(numThreads);

  // Write the transformed data
  cout<<"Generating Image File MyAfter1D.txt"<<endl;
  image.SaveImageData("MyAfter1D.txt", Data, Width, Height);

  // Transpose
  Transpose(Data, buffer);

  /***** column fft *****/

  BitReverseOrder(Data, buffer);

  MyBarrier_Init();

  for(int i = 0; i < numThreads; i++)
    pthread_create(&threads[i], 0, Transform2DTHread, (void*)i);

  MyBarrier(numThreads);

  Transpose(Data, buffer);

  cout<<"Generating Image File MyAfter2D.txt"<<endl;
  image.SaveImageData("MyAfter2D.txt", Data, Width, Height);


  /****************************   Inverse    ********************************/

  inv = true;

  prepare_W_inv();

  BitReverseOrder(Data, buffer);

  MyBarrier_Init();

  for(int i = 0; i < numThreads; i++)
    pthread_create(&threads[i], 0, Transform2DTHread, (void*)i);

  MyBarrier(numThreads);

  Transpose(Data, buffer);

  MyBarrier_Init();
  
  BitReverseOrder(Data, buffer);

  for(int i = 0; i < numThreads; i++)
    pthread_create(&threads[i], 0, Transform2DTHread, (void*)i);

  MyBarrier(numThreads);

  Transpose(Data, buffer);

  for(int i = 0; i < N * N; i++)
  {
    if(Data[i].real < 0.9  && Data[i].imag < 0.9)
    {
      Data[i].real = 0;
      Data[i].imag = 0;
    }
  }

  cout<<"Generating Image File MyAfterInverse.txt"<<endl;
  image.SaveImageData("MyAfterInverse.txt", Data, Width, Height);

  /************************** Garbage Collection ****************************/

  delete[] W;
  delete[] buffer;

}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  // MPI initialization here
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
