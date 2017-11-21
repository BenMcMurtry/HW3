//  Created by Ben McMurtry on 24/11/2016.

/*
 A simple program investigating random number generators (included currently is random(), the Approximate method, the Box-Muller method, and the Rejection method)
 Given a random number generator and its corresponding PDF, this program can:
 (a) Calculate the maximum, minimum, mean, and standard deviation of N samples provided by the RNG(), and the standard error of the mean.
 (b) Output data to a file to be plotted using Gnuplot and a script, in a form that allows it to produce a histogram of values and the function PDF(x) appropriately scaled to facilitate comparison on the same plot.
 (c) Time the random number generators to see how long it takes them to produce a number of random variates.
*/


#include <stdio.h>
#include <stdlib.h> /* srandom, random */
#include <string.h>
#include <limits.h> /* for PATH_MAX */
#include <math.h>
#include <errno.h>
#include <time.h> /* for clock(), clock_t */
#include <complex.h> /* for cexp */

/* Definitions for RNG function, and sample size. */
#define RANDOM_MAX 0x7FFFFFFF
#define LOWER 3.0
#define UPPER 7.0
#define SAMPLE_SIZE 10000

/* Definitions for storing sample output in a file and plotting it with gnuplot. */
//#define GNUPLOT_EXE       "/opt/local/bin/gnuplot" /* This is my gnuplot path on my mac */
#define GNUPLOT_EXE       "gnuplot"               /* Your path for gnuplot here! */
#define OUTPUT_FILE       "RNG_DATA.dat"
#define GNUPLOT_SCRIPT    "RNG_gnuplot.script"

/* This is chosen at run-time by the user */
typedef enum { NONE, RANDOM, APPROXIMATE, BOX_MULLER, REJECTION } Method;


/* Function to create a random number between min and max using random() from the POSIX library. */
double RNG(double min, double max) {
    return( (random() / (double)RANDOM_MAX) * (max - min) + min );
}


/* Function takes x values to produce a uniform probability density between LOWER and UPPER. */
double PDF(double x) {
    if ( ( LOWER < x ) && ( x < UPPER ) ) {
        return 1.0 / (UPPER - LOWER);
    } else if (x == LOWER || x == UPPER ) {
        return 0.5 / (UPPER - LOWER);
    }
    return 0.0; /* Note that zero is returned if x is NaN */
}


/* Function to create a random variate with gaussian(mean, variance) density using the Approximate method */
double approximate_method(double mean, double variance) {
    double sum = 0.0;
    for (int i = 0; i < 12; i++) {
        sum += sqrt(variance) * RNG(0.0, 1.0);
    }
    return (sum - (6.0 * sqrt(variance)) + mean);
}


/* Function to create 2 random variates with gaussian(mean, variance) density using the Box-Muller Method */
void box_muller_method(double *y1, double *y2, double mean, double variance) {
    double u1 = RNG(0.0, 1.0);
    double u2 = RNG(0.0, 1.0);
    double complex z1 = cexp(I * 2 * M_PI * u2);
    *y1 = (sqrt(- 2.0 * log(u1)) * creal(z1)) * sqrt(variance) + mean;
    *y2 = (sqrt(- 2.0 * log(u1)) * cimag(z1)) * sqrt(variance) + mean;
}


/* This is the function of the normal distribution */
double normal(double x) {
    double p = (1.0 / sqrt(2.0 * M_PI)) * exp(-pow(x, 2.0) / 2.0);
    return p;
}


/* This function is chosen to be larger than the normal distribution for all x, but its CDF is easier to invert */
double majorising_func(double x) {
    double a = 0.50 / (1.0 + pow(x, 2.0));
    return a;
}


/* This function performs the inverse of the CDF of the normalised majorising function on uniform U[0,1] variates, creating random variates with a density of the normalised majorising function. */
double inverseCDF() {
    double x = RNG(0.0, 1.0);
    double result = tan(M_PI * (x - 0.5));
    return result;
}


/* Function takes random variates Y with denisty of the normalised majorising function, and creates uniform U[0,1] variates. When the U[0,1] variate is less than the ratio of the normal function(Y) over the majorising function(Y) (Monte Carlo), it converts the Y values from normal to gaussian density and retuns them */
double rejection_method(double mean, double variance) {
    double y = inverseCDF();
    double u = RNG(0.0, 1.0);
    if (u <= normal(y) / majorising_func(y)) { // If the RV is <= normal(y) / majoriser(y), then add it to sample.
        return y * sqrt(variance) + mean;
    } else {
        return rejection_method(mean, variance);
    }
}


/* Function to produce the Probability Density Function of a Gaussian with desired mean and variance */
double PDF_Gaussian(double x, double mean, double variance) {
    double result = (1.0 / sqrt(2.0 * variance * M_PI)) * exp(-pow((x - mean), 2.0) / (2.0 * variance));
    return result;
}


/* This function is required for qsort, helping it decide in what order to place two numbers. */
int compare (const void * a, const void * b)
{
    if ( *(double*)a <  *(double*)b ) return -1;
    if ( *(double*)a == *(double*)b ) return 0;
    if ( *(double*)a >  *(double*)b ) return 1;
    return 0;
}


/* Function to find the mean of an array of doubles. */
double mean_finder(double array[]) {
    double mean = 0.0;
    for (int i = 0; i < SAMPLE_SIZE; i++) {
        mean += array[i];
    }
    mean = mean / SAMPLE_SIZE;
    return mean;
}


/* Function to find the standard deviation of an array of doubles. */
double standard_dev_finder(double array[], double mean) {
    double sigma = 0.0;
    double variance = 0.0;
    for (int i = 0; i < SAMPLE_SIZE; i++) {
        variance += pow((array[i] - mean), 2);
    }
    variance = variance / SAMPLE_SIZE;
    sigma = sqrt(variance);
    return sigma;
}


/* Finds mean, std dev and std error of a sample, and then will use GNUPLOT_SCRIPT to plot sample using GNUPLOT_EXE. */
void investigator(double sample[SAMPLE_SIZE]) {
    printf("The number of datapoints in the sample is %d\n", SAMPLE_SIZE);
    printf("The minimum of the sample is %lg\n", sample[0]);
    printf("The maximum of the sample is %lg\n", sample [SAMPLE_SIZE - 1]);
    
    double mean = mean_finder(sample);
    double std_dev = standard_dev_finder(sample, mean);
    double std_err = std_dev / sqrt(SAMPLE_SIZE);
    printf("The mean of the sample %lg\n", mean);
    printf("The standard deviation of the sample is %lg\n", std_dev);
    printf("The standard error of the mean of the sample is %lg\n", std_err);
    
    char command[PATH_MAX];
    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
    system( command );
}


/* Function to fill an array of doubles with a Random Variate generator, with chosen Method "mychoice". */
void array_maker(double sample[SAMPLE_SIZE], Method mychoice) {
    if (mychoice == RANDOM) {
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            sample[i] = RNG(UPPER, LOWER);
        }
    }
    if (mychoice == APPROXIMATE) {
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            sample[i] = approximate_method(5.0, 2.0);
        }
    }
    if (mychoice == BOX_MULLER) {
        for (int i = 0; i < SAMPLE_SIZE; i+=2) {
            box_muller_method(&sample[i], &sample[i+1], 5.0, 2.0);
        }
    }
    if (mychoice == REJECTION) {
        for (int i = 0; i < SAMPLE_SIZE; i++) {
            sample[i] = rejection_method(5.0, 2.0);
        }
    }
}



int main() {
    double sample[SAMPLE_SIZE];
    char m;
    printf("Please enter 'p' for RNG, 'a' for Approximate, 'b' for Box-Muller, or 'r' for Rejection method: ");
    m = getchar();
    while ((m != 'p') & (m!= 'a') * (m != 'b') & (m != 'r')) {
        printf("\nInvalid character for method choice.\n");
        printf("Please enter 'p' for RNG, 'a' for Approximate, 'b' for Box-Muller, or 'r' for Rejection method: ");
        while ((m = getchar()) != '\n') {}
        m = getchar();
    }
    Method mychoice = NONE;
    if (m == 'p') {mychoice = RANDOM;}
    if (m == 'a') {mychoice = APPROXIMATE;}
    if (m == 'b') {mychoice = BOX_MULLER;}
    if (m == 'r') {mychoice = REJECTION;}
    
    //srandom((unsigned int)time(NULL));
    clock_t begin = clock();
    array_maker(sample, mychoice);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / (double)CLOCKS_PER_SEC;
    printf("Time taken was %lg seconds\n", time_spent);
    
    qsort (sample, SAMPLE_SIZE, sizeof(double), compare);
    
    FILE *output = fopen(OUTPUT_FILE, "w");
    if (!output) {
        fprintf(stderr, "Error: Could not open file '%s'.\n", OUTPUT_FILE);
        exit(1);
    }
    for (int i = 0; i < SAMPLE_SIZE; i++) {
        if (mychoice == RANDOM) {
            fprintf(output, "%-12lg%-12lg\n", sample[i], PDF(sample[i]));
        } else {
            fprintf(output, "%-12lg%-12lg\n", sample[i], PDF_Gaussian(sample[i], 5.0, 2.0));
        }
    }
    fclose(output);
    
    investigator(sample);

    return 0;
}
