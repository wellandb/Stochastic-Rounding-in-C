#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

union DoubleToInt {
  double dVal;
  uint64_t iVal;
};

/*
  A function that rounds a binary64 value to a binary32 value
  stochastically. Implemented by treating FP number representations
  as integer values.
*/
float SR(double x) {

  union DoubleToInt temp;
  temp.dVal = x;
  uint32_t r = rand() & 0x1FFFFFFF;
  temp.iVal += r;
  temp.iVal = temp.iVal & 0xFFFFFFFFE0000000;

  return (float)temp.dVal;
}


/* --------------------------------- */
/*              PART 1               */
/* --------------------------------- */
  
// Implement SR_alternative according to the Eqn 1.
float SR_alternative(double x) {
  // TODO
  // Create Double to int
  union DoubleToInt temp;
  temp.dVal = x;
  // Work out closest binary32 number
  float closest = (float)x;
  // Initialize the up and down binary32 numbers from x
  float down, up;
  // Choose up and down
  if (closest > x) {
    down = nextafterf(closest, -INFINITY);
    up = closest;
  }
  else {
    down = closest;
    up = nextafterf(closest, INFINITY);
  }
  // Create random number from 0 to 99
  float rnd = rand() % 100;
  // P = x-RZ(x)/RA(x)-RZ(x), then x100 to get a percentage out of 100 to compare with random number
  double p = (x-down)/(up-down) * 100;

  // if P < p return RA else RZ
  if (rnd < p){
    temp.dVal = up;
  }else{
    temp.dVal = down;
  }

  return (float)temp.dVal;
}

// fastTwoSum, adds 2 numbers while getting the error
void fastTwoSum (float a, float b, float *s, float *t) {
  float temp;

  *s = a + b;
  temp = *s - a;
  *t = b - temp;
}


const long int K = 5000000;

int main() {
  
  // Create time and the seed random number to create random number
  time_t ti;
  srand((unsigned) time(&ti));

  // An arbitrary value for rounding.
  double sample = M_PI;
  double avg = 0;
  double avg_alternative = 0;

  // Calculate the neighbouring binary32 values.
  float closest = (float)sample;
  float down, up;
  if (closest > sample) {
    down = nextafterf(closest, -INFINITY);
    up = closest;
  }
  else {
    down = closest;
    up = nextafterf(closest, INFINITY);
  }

  // Round many times, and calculate the average values as well as count
  // the numbers of times rounding was up/down.
  int countup = 0;
  int countdown = 0;

  // Error inits
  double sr_err[K/1000], sralt_err[K/1000];
  int counter = 0;
  for (int i = 1; i <= K; i++) {
    // TODO
    float tempVal = SR(sample);
    // SR_alt
    float tempVal_alt = SR_alternative(sample);
    // Countup and countdown
    if (tempVal_alt == up){
      countdown += 1;
    }else{
      countup += 1;
    }
    // avg value calculations
    avg += tempVal;
    avg_alternative += tempVal_alt;

    // Error calcs
    if (i % 100000 == 0){
      sr_err[counter] = fabs(sample - avg/i);
      sralt_err[counter] = fabs(sample - avg_alternative/i);
      counter += 1;
    }
  }
  // Avg value final calculations
  avg = avg/K;
  avg_alternative = avg_alternative/K;

  // Print out some useful stats.
  printf("Value being rounded:           %.60f \n", sample);
  printf("SR average value:              %.60f \n", avg);
  printf("SR_alternative average value:  %.60f \n", avg_alternative);
  printf("Binary32 value before:         %.60f \n", down);
  printf("Binary32 value after:          %.60f \n", up);
  printf("Closest binary32:              %.60f \n", closest);

  
  // Check that SR_alternative function is correct by comparing the probabilities
  // of rounding up/down, and the expected probability. Print them out
  // below.
  // TODO
  // Count up and count down stats
  printf("Number of rounding up:         %d \n", countup);
  printf("Number of rounding down:       %d \n", countdown);

  // Probability of countup/down
  double p_down = (double) countdown * 100 / (K);
  double exp_p_down = 100*fabs(sample-down)/fabs(up-down);
  double p_up = (double) countup * 100 / (K) ;
  double exp_p_up = 100*fabs(up-sample)/fabs(up-down);

  printf("Probability of rounding up:    %f%% \n", p_up);
  printf("Probability of rounding down:  %f%% \n", p_down);
  printf("Expected Probability of rounding up:    %f%% \n", exp_p_up);
  printf("Expected Probability of rounding down:  %f%% \n", exp_p_down);

  
  /* --------------------------------- */
  /*              PART 2               */
  /* --------------------------------- */
  long int N = 500000000;
  float fharmonic = 0;
  float fharmonic_sr = 0;
  float fharmonic_comp = 0;
  double dharmonic = 0;

  // Error term in the compensated summation.
  float t = 0;


  // Other error terms
  double f_err[N/1000000]; 
  double fsr_err[N/1000000];
  double fcomp_err[N/1000000];
  counter = 0;
  int stag = 0;
  
  for (int i = 1; i <= N; i++) {
    // Recursive sum, binary32 RN
    fharmonic += (float)1/i;

    // Other summation methods, TODO.

    //Compensated summation, binary32
    fastTwoSum(fharmonic_comp, (float)1/i, &fharmonic_comp, &t);

    //Recursive summation, binary64
    dharmonic += (double)1/i;

    //Recursive summation with SR, binary32
    fharmonic_sr = SR(dharmonic);

    //Error calcs
    if(i % 1000000 == 0){
      f_err[counter] = fabs(dharmonic - fharmonic);
      fsr_err[counter] = fabs(dharmonic - fharmonic_sr);
      fcomp_err[counter] = fabs(dharmonic - fharmonic_comp);
      counter += 1;
    }
    // Stagnation calcs
    if(fharmonic == fharmonic + (float) 1/i & stag == 0){
      stag = i;
    }
  }

  printf("\n");
  printf("Values of the harmonic series after %ld iterations \n", N);
  printf("Recursive summation, binary32:          %.30f \n", fharmonic);
  printf("Recursive summation with SR, binary32:  %.30f \n", fharmonic_sr);
  printf("Compensated summation, binary32:        %.30f \n", fharmonic_comp);
  printf("Recursive summation, binary64:          %.30f \n", dharmonic);
  printf("Stagnation at i = %d", stag);

  /* --------------------------------- */
  /*              PART 3               */
  /* --------------------------------- */
  

  // TODO
  // Calculating the reciprocals of the fibonacci sequence
  long int X = 500000000;
  float fib_rep = 0;
  float fib_rep_sr = 0;
  float fib_rep_comp = 0;
  double fib_rep_d = 0;

  // error calc
  float error = 0;
  float err[X/1000000];
  // Other error terms
  double fib_err[X/1000000]; 
  double fibsr_err[X/1000000];
  counter = 0;

  // fibonacci sequence initalise
  int fib1 = 1;
  int fib2 = 1;

  for (int i = 1; i <= X; i++) {
    // reciprocals of the fibonacci sequence binary32
    fib_rep += (float) 1/fib1;

    // reciprocals of the fibonacci sequence binary64
    fib_rep_d += (double) 1/fib1;

    //Recursive summation with SR, binary32
    fib_rep_sr = SR(fib_rep_d);

    //Compensated summation, binary32
    fastTwoSum(fib_rep_comp, (float)1/fib1, &fib_rep_comp, &error);

    // fibonacci sequence
    int temp = fib2;
    fib2 += fib1;
    fib1 = temp;

     //Error calcs
    if(i % 1000000 == 0){
      fib_err[counter] = fabs(fib_rep_d - fib_rep);
      fibsr_err[counter] = fabs(fib_rep_d - fib_rep_sr);
      err[counter] = fabs(error);
      counter += 1;
    } 

  }
  


  printf("\n");
  printf("Values of the reciprocals of the fibonacci sequence after %ld iterations \n", X);
  printf("Recursive summation, binary32:          %.30f \n", fib_rep);
  printf("Recursive summation with SR, binary32:  %.30f \n", fib_rep_sr);
  printf("Compensated summation, binary32:        %.30f \n", fib_rep_comp);
  printf("Recursive summation, binary64:          %.30f \n", fib_rep_d);

  /* ----- GRAPH ----- */

  char * commandsForGnuplot[] = {"set title \"Stochastic Rounding: Absolute error vs Rounding iteration\"","set xlabel 'Rounding Iteration'","set ylabel 'Absolute error'","set style line 1 lt rgb '#FF0000' lw 3 pt 6","plot 'SR Error', 'SR alt Error'","exit"};
  double xvals1[K/100000];
  for(int i = 0; i < K/100000; i++){
    xvals1[i] = 100000*i;
  }
  FILE * temp = fopen("SR Error", "w");
  FILE * temp2 = fopen("SR alt Error", "w");
  /*Opens an interface that one can use to send commands as if they were typing into the
    *     gnuplot command line.  "The -persistent" keeps the plot open even after your
    *     C program terminates.
    */
  FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
  int i;
  for (i=0; i < K/100000; i++)
  {
    fprintf(temp, "%lf %.60lf \n", xvals1[i], sr_err[i]); //Write the data to a temporary file
    fprintf(temp2, "%lf %.60lf \n", xvals1[i], sralt_err[i]); //Write the data to a temporary file
  }
  for (i=0; i < 6; i++)
  {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
  }
 
  // GRAPH 2
  char * commandsForGnuplot1[] = {"set title \"Stochastic Rounding: Harmonic Series\"","set xlabel 'Rounding Iteration'","set ylabel 'Absolute error'","set style line 1 lt rgb '#FF0000' lw 3 pt 6","set logscale y 2", "plot 'harmonic error' ls 1, 'harmonic sr error' ls 2, 'harmonic compensated error' ls 3"};
  
  double xvals2[N/1000000];
  for(int i = 0; i < N/1000000; i++){
    xvals2[i] = 1000000*i;
  }


  FILE * temp3 = fopen("harmonic error", "w");
  FILE * temp4 = fopen("harmonic sr error", "w");
  FILE * temp5 = fopen("harmonic compensated error", "w");
 

  FILE * gnuplotPipe1 = popen ("gnuplot -persistent", "w");

  for(i=0; i < N/1000000; i++){
    fprintf(temp3, "%lf %.60lf \n", xvals2[i], f_err[i]); //Write the data to a temporary file
    fprintf(temp4, "%lf %.60lf \n", xvals2[i], fsr_err[i]); //Write the data to a temporary file
    fprintf(temp5, "%lf %.60lf \n", xvals2[i], fcomp_err[i]); //Write the data to a temporary file
  }
  for (i=0; i < 6; i++)
  {
    fprintf(gnuplotPipe1, "%s \n", commandsForGnuplot1[i]); //Send commands to gnuplot one by one.
  }

    // GRAPH 3
  char * commandsForGnuplot2[] = {"set title \"Stochastic Rounding: Reciprocols of the Fibonacci Series\"","set xlabel 'Rounding Iteration'","set ylabel 'Absolute error'","set style line 1 lt rgb '#FF0000' lw 3 pt 6","set logscale y 10","plot 'fibonacci error', 'fibonacci sr error', 'fibonacci compensated error'"};
  
  FILE * temp6 = fopen("fibonacci error", "w");
  FILE * temp7 = fopen("fibonacci sr error", "w");
  FILE * temp8 = fopen("fibonacci compensated error", "w");

  FILE * gnuplotPipe2 = popen ("gnuplot -persistent", "w");

  for(i=0; i < X/1000000; i++){
    fprintf(temp6, "%lf %.60lf \n", xvals2[i], fib_err[i]); //Write the data to a temporary file
    fprintf(temp7, "%lf %.60lf \n", xvals2[i], fibsr_err[i]); //Write the data to a temporary file
    fprintf(temp8, "%lf %.60lf \n", xvals2[i], err[i]); //Write the data to a temporary file
  }
  for (i=0; i < 6; i++)
  {
    fprintf(gnuplotPipe2, "%s \n", commandsForGnuplot2[i]); //Send commands to gnuplot one by one.
  }
  return 0;
}
