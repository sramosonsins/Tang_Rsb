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

/*
calc_iES.slow <- function(geno,lox) {#function(EHH,lox)
    if( nrow(geno) != length(lox) ) { stop("Number of positions given does not agree with number of markers.\n") }
 
    L <- length(lox)
    iES <- rep(NA, L)
    x = lox[2:L] - lox[1:(L-1)]
    
    for(i in 1:L) {
        EHHSa <- calc_EHHS.pos(i,geno)
        y = EHHSa[1:(L-1)] + EHHSa[2:L]
        if( !all(is.na(y)) ) {
            iES[i] = sum(y*x, na.rm=T) / 2
        }; rm(y)
    }
    iES
}
*/
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
/*
calc_EHHS.pos <- function(pos,geno,thresh=0.1) {
    M <- nrow(geno)
    EHH <- rep(NA,M)
    i <- pos
 
    Ii <- sum(geno[i,]==0)
    EHH[i] = 1
    
    ## left-flank
    j = i-1
    while( j >= 1 ) {
        Ij <- sum(colSums(geno[j:i,])==0)
        cur.EHH = Ij/Ii
        if (is.na(cur.EHH) | cur.EHH < thresh) { break } else {
            EHH[j] <- cur.EHH
        }
        j = j-1
    }
    
    ## right-flank
    j = i + 1
    while( j <= M ) {
        Ij <- sum(colSums(geno[i:j,])==0)
        cur.EHH = Ij/Ii
        if (is.na(cur.EHH) | cur.EHH < thresh) { break } else {
            EHH[j] <- cur.EHH
        }
        j = j+1
    }
    return(EHH)
}

*/

