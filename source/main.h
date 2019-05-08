//
//  main.h
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 29/01/2019.
//  From an original R script of Eva KF Chan
//  Function to calculate the Rsb statistic for a given chromosome as described in:
//  Tang K, Thornton KR, Stoneking M (2007) A New Approach for Using Genome Scans to Detect Recent Positive Selection in the Human Genome . PLoS Biol 5(7): e171 doi:10.1371/journal.pbio.0050171
//

#ifndef main_h
#define main_h

#define TANG_SOFTW  "\nSoftware for calculating iES and Rsb statistics." \
"\nfollowing Tang, Thornton & Stoneking, PloS Biology 2007." \
"\nversion 20190508"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ran1.h"


int read_row(FILE *plink_file,char *chr_name, double *lox, int **geno, int geno_cols, long int row);

void impute_genotypes(int **geno, long int *geno_rows, int geno_cols);

void usage(void);

void calc_iES_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iES);

void calc_EHHS_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, double *EHH);

#endif /* main_h */
