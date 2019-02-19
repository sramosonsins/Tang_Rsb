//
//  imputation.c
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 31/01/2019.
//
#include "ran1.h"

void impute_genotypes(int **geno, long int *geno_rows, int geno_cols)
{
    int pj[3][3];
    int pjt;
    int gj[3];
    double prob_pj[3];
    long int i,i1;
    int j,k,l;
    int *is_na;
    double ran_val;
    void remove_row(long int i, int **geno, long int *geno_rows, int geno_cols);
    
    is_na = (int *)calloc(geno_cols,sizeof(int));
    
    for(i=0;i<*geno_rows;i++) {
        //init arrays
        for(k=0;k<3;k++) {
            gj[k] = 0;
            for(l=0;l<3;l++) {
                pj[k][l] = 0;
            }
        }
        //estimate probs
        k=0;
        if(i==0) i1=0;
        else i1 = i-1;
        for(j=0;j<geno_cols;j++) {
            gj[0] += (geno[j][i1]==0);
            gj[1] += (geno[j][i1]==1);
            gj[2] += (geno[j][i1]==2);
            pj[0][0] += (geno[j][i1]==0 && geno[j][i]==0);
            pj[0][1] += (geno[j][i1]==0 && geno[j][i]==1);
            pj[0][2] += (geno[j][i1]==0 && geno[j][i]==2);
            pj[1][0] += (geno[j][i1]==1 && geno[j][i]==0);
            pj[1][1] += (geno[j][i1]==1 && geno[j][i]==1);
            pj[1][2] += (geno[j][i1]==1 && geno[j][i]==2);
            pj[2][0] += (geno[j][i1]==2 && geno[j][i]==0);
            pj[2][1] += (geno[j][i1]==2 && geno[j][i]==1);
            pj[2][2] += (geno[j][i1]==2 && geno[j][i]==2);
            if(geno[j][i]==9) is_na[k++] = j;
        }
        pjt = geno_cols - k;
        if(i==0) {gj[0]=gj[1]=gj[2]=1;}
        prob_pj[0] = gj[0]*pj[0][0] + gj[1]*pj[0][1] + gj[2]*pj[0][2];
        prob_pj[1] = gj[0]*pj[1][0] + gj[1]*pj[1][1] + gj[2]*pj[1][2];
        prob_pj[2] = gj[0]*pj[2][0] + gj[1]*pj[2][1] + gj[2]*pj[2][2];
        
        //assign genotypes
        if(pjt==0 || prob_pj[0] + prob_pj[1] + prob_pj[2] == 0) {
            remove_row(i,geno,geno_rows,geno_cols);
        }
        else {
            ran_val = ran1();
            for(j=0;j<k;j++) {
                if((double)prob_pj[0]/(double)pjt > ran_val) {
                    geno[is_na[j]][i] = 0;
                }
                else {
                    if((double)prob_pj[1]/(double)(pjt-prob_pj[0]) > ran_val) {
                        geno[is_na[j]][i] = 1;
                    }
                    else {
                        if((double)prob_pj[2] > 0.0)
                            geno[is_na[j]][i] = 2;
                        else {
                            remove_row(i,geno,geno_rows,geno_cols);
                        }
                    }
                }
            }
        }
    }
    free(is_na);
}

void remove_row(long int row, int **geno, long int *geno_rows, int geno_cols)
{
    long int i;
    int j;
    
    for(j=0;j<geno_cols;j++) {
        for(i=row;i<*geno_rows-1;i++) {
            geno[j][i] = geno[j][i+1];
        }
        geno[j] = (int *)realloc(geno[j],(*geno_rows - 1)*sizeof(int));
    }
    *geno_rows = *geno_rows - 1;
}
