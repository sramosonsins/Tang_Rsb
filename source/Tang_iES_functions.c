//
//  main.c
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 29/01/2019.
//  Copyright © 2019 Sebastian Ramos-Onsins. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

void calc_iES_slow(int **geno, double *lox, long int *geno_rows, int *geno_cols, double *thresh, double *iES)
{
    void calc_EHHS_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, long int *min, long int *max, double *EHH);
    
    double *x;
    double *y;
    double *EHHSa;
    
    long int L;
    long int i,j,k;
    long int min,max;
    
    L = *geno_rows;
    x = (double *)calloc(L-1,sizeof(double));
    y = (double *)calloc(L-1,sizeof(double));
    EHHSa = (double *)calloc(L,sizeof(double));
    //for(i=0;i<L;i++) EHHSa[i]=-1.0;
    
    for(i=0;i<L-1;i++) {x[i] = lox[i+1] - lox[i];}
    for(i=0;i<L;i++) {
        calc_EHHS_pos(&i,geno,geno_rows,geno_cols,thresh,&min,&max,EHHSa);
        for(j=min;j<max;j++) {
            y[j] = EHHSa[j+1] + EHHSa[j];
        }
        for(k=min;k<max;k++) {
            iES[i] += y[k] * x[k] / 2.0;
        }
        for(j=min;j<max;j++) {y[j] = 0.0;}
    }
    free(x);
    free(y);
    free(EHHSa);
    return;
}

void calc_EHHS_pos(long int *i, int **geno, long int *geno_rows, int *geno_cols, double *thresh, long int *min, long int *max, double *EHH)
{
    long int M;
    long int pos;
    long int j,k;
    long int Ii,Ij;
    long int *Ic;
    double cur_EHH;
    
    pos = *i;
    M = *geno_rows;
    Ii = 0;
    for(k=0;k<*geno_cols;k++) {Ii += (geno[k][pos]==0);}
    if(Ii == 0) {
        EHH[pos] = -1;
        *min = *max = pos;
        return;
    }
    EHH[pos] = 1.0;
    
    //left-flank
    Ic = (long int *)calloc(*geno_cols,sizeof(long int));

    for(k=0;k<*geno_cols;k++) {Ic[k] += geno[k][pos];}
    j = pos - 1;
    while( j >= 0) {
        Ij = 0;
        for(k=0;k<*geno_cols;k++) {
            Ic[k] += geno[k][j];
            Ij += (Ic[k]==0);
        }
        cur_EHH = (double)Ij / (double)Ii;
        if(cur_EHH < *thresh) {
            break;
        }
        else {
            EHH[j] = cur_EHH;
        }
        j -= 1;
    }
    *min = j+1;
    
    //right-flank
    for(k=0;k<*geno_cols;k++) {Ic[k] = geno[k][pos];}
    j = pos + 1;
    while( j < M) {
        Ij = 0;
        for(k=0;k<*geno_cols;k++) {
            Ic[k] += geno[k][j];
            Ij += (Ic[k]==0);
        }
        cur_EHH = (double)Ij / (double)Ii;
        if(cur_EHH < *thresh) {
            break;
        }
        else {
            EHH[j] = cur_EHH;
        }
        j += 1;
    }
    *max = j-1;
    free(Ic);
    
    return;
}

