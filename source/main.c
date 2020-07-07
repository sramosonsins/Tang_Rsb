//
//  main.c
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 29/01/2019.
//  Copyright © 2019 Sebastian Ramos-Onsins. All rights reserved.
//

#include "main.h"

int main(int arg, const char *argv[])
{
    double *lox;
    int **geno;
    long int geno_rows;
    int geno_cols;
    int **pop_geno;
    int pop_geno_cols;
    double thresh;
    double **all_iES;
    long int L;
    int N;
    long int i,h;
    int j,k,l;
    int npops;
    int *popsize;
    int popcum;
    
    double *mean;
    double *median;
    double *med;
    double *sd;
    double **all_Rsb;
    long int *nit;
    char **pop_name;

    FILE *plink_file = 0;
    FILE *genot_file = 0;
    FILE *results_file = 0;
    int c;
    long int row;
    char chr_name[11]; //the file must contain the same ID for all positions!
    char file_out[1024];
    char file_in[1024];
    int argc = 1;

    if(arg < 8) {
        printf("\nError introducing arguments\n");
        usage();
        exit(1);
    }
    
    printf(TANG_SOFTW);
    printf("\n\nTang_stats ");
    while(argc < arg) {
        printf("%s ",argv[argc]);
        argc++;
    }
    printf("\n");
    
    //define variables and arrays:
    strcpy( file_in, argv[1]);
    geno_rows = L = atol(argv[2]);
    geno_cols = N = atoi(argv[3]);
    thresh = atof(argv[4]);
    init_seed1(atol(argv[5]));
    npops = atoi(argv[6]);
    popsize = (int *)calloc(npops,sizeof(int));
    for(i=0;i<npops;i++) {
        popsize[i] = atoi(argv[7+i]);
    }
    pop_name = (char **)calloc(npops,sizeof(char *));
    for(i=0;i<npops;i++) {
        pop_name[i] = (char *)calloc(100,sizeof(char));
        strcpy(pop_name[i], argv[7+npops+i]);
    }

    if (!(plink_file = fopen(file_in,"r"))) {
        printf("Error reading the input file %s\n",file_in);
        exit(1);
    }
    
    //geno is transposed to facilitate the analysis
    geno = (int **)calloc(N,sizeof(int *));
    for(j=0;j<N;j++) {geno[j] = (int *)calloc(L,sizeof(int));}
    lox = (double *)calloc(L,sizeof(double));
    all_iES = (double **)calloc(npops,sizeof(double *));
    for(i=0;i<npops;i++) {all_iES[i] = (double *)calloc(L,sizeof(double));}
    
    mean = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    median = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    sd = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    nit = (long int *)calloc(npops*(npops-1)/2,sizeof(long int));

    all_Rsb = (double **)calloc(npops*(npops-1)/2,sizeof(double *));
    for(i=0;i<npops*(npops-1)/2;i++) {all_Rsb[i] = (double *)calloc(L,sizeof(double));}

    //read input data: skip header
    printf("\nReading input file...");
    while((c=getc(plink_file))!='\n' && c!='\r');
    row = 0;
    for(row=0;row<L;row++) {
        if((c = read_row(plink_file,chr_name,lox,geno,geno_cols,row))==0) {
            printf("\nError: input file has less rows or more cols than defined. nrows: %ld\n",row+1);
            exit(1);
        }
    }
    fclose(plink_file);
    
    //imputation
    printf("\nimputing missing values (previously assigned as genotype=9)...\n");
    impute_genotypes(geno, L, geno_cols);

    //writing imputed file
    memset(file_out, '\0', 1024);
    strcat(file_out,file_in);
    strcat(file_out,"_imputed.txt\0");
    printf("\nWriting imputed genotype file %s...\n",file_out);
        //header
    if (!(genot_file = fopen(file_out,"w"))) {
        printf("Error writing the imputed file %s\n",file_in);
        exit(1);
    }
    fprintf(genot_file,"CHR\tPOS");
    for(j=0;j<npops;j++) {
        for(k=0;k<popsize[j];k++) {
            fprintf(genot_file,"\t%s%d_IND%d",pop_name[j],j+1,k+1);
        }
    }
    fprintf(genot_file,"\n");
        //genotypes
    for(i=0;i<L;i++) {
        fprintf(genot_file,"%s\t%f",chr_name,lox[i]);
        for(j=0;j<N;j++) {
            fprintf(genot_file,"\t%d",geno[j][i]);
        }
        fprintf(genot_file,"\n");
    }
    fclose(genot_file);
    
    //converting 2 homozygotes to 0:
    printf("\nConverting all homozygotes to value 0...");
    for(j=0;j<N;j++) {
        for(i=0;i<L;i++) {
            if(geno[j][i]==2) geno[j][i]=0;
        }
    }
    
    //calculation iES
    popcum = 0;
    for(j=0;j<npops;j++) {
        printf("\nComputing iES for pop %d of %d ...",j+1,npops);
        pop_geno_cols = popsize[j];
        pop_geno = &geno[popcum];
        calc_iES_slow(pop_geno, lox, &geno_rows, &pop_geno_cols, &thresh, all_iES[j]);
        popcum += popsize[j];
    }
    
    //calculation Rsb, median, mean and sd:
    printf("\nComputing Rsb...");
    for(i=0;i<L;i++) {
        l=0;
        for(j=0;j<npops-1;j++) {
            for(k=j+1;k<npops;k++) {
                if(all_iES[j][i] > 0.0 && all_iES[k][i] > 0.0) {
                    all_Rsb[l][i] = log(all_iES[j][i])-log(all_iES[k][i]);
                    mean[l] = mean[l] + all_Rsb[l][i];
                    sd[l] = sd[l] + all_Rsb[l][i] * all_Rsb[l][i];
                    nit[l] += 1;
                 }
                else {
                    all_Rsb[l][i] = 1234567890;
                }
                l++;
            }
        }
    }
    // median ...
    l=0;
    med = (double *)calloc(L,sizeof(double));
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            h = 0;
            for(i=0;i<L;i++) {
                if(all_Rsb[l][i] != 1234567890) {
                    med[h++] = all_Rsb[l][i];
                }
            }
            qsort(med,nit[l],sizeof(double),compare_);
            if((double)nit[l]/2.0 == nit[l]/2.0) {
                median[l] = (med[(long int)((double)nit[l]/2.0 - 1.0)] +
                             med[(long int)((double)nit[l]/2.0)]) / 2.0;
            }
            else median[l] = med[nit[l]/2];
            l++;
        }
    }
    // mean, sd ...
    l=0;
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            mean[l] = mean[l]/(double)nit[l];
            sd[l] = sqrt(sd[l]/((double)nit[l]-1.0) - mean[l]*mean[l] * (double)nit[l]/((double)nit[l]-1.0)); //smapled sd
            l++;
        }
    }
    
    //writing output file
    file_out[0] = '\0';
    strcat(file_out,file_in);
    strcat(file_out,"_Results_Tang.txt\0");
    printf("\nWriting results in the output file %s...",file_out);

    if (!(results_file = fopen(file_out,"w"))) {
        printf("Error writing the input file %s\n",file_out);
        exit(1);
    }
    //header
    fprintf(results_file,"Position\t");
    for(j=0;j<npops;j++) fprintf(results_file,"iES_%s\t",pop_name[j]);
    for(j=0;j<npops;j++) fprintf(results_file,"log(iES_%s)\t",pop_name[j]);
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            fprintf(results_file,"Rsb(%s/%s)\t",pop_name[j],pop_name[k]);
        }
    }
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            fprintf(results_file,"RsbN(%s/%s)\t",pop_name[j],pop_name[k]);
        }
    }
    fprintf(results_file,"\n");
    //data
    for(i=0;i<L;i++) {
        fprintf(results_file,"%f\t",lox[i]);
        for(j=0;j<npops;j++) {
            if(all_iES[j][i] > 0.0) fprintf(results_file,"%f\t",all_iES[j][i]);
            else fprintf(results_file,"NA\t");
        }
        for(j=0;j<npops;j++) {
            if(all_iES[j][i] > 0.0) fprintf(results_file,"%f\t",log(all_iES[j][i]));
            else fprintf(results_file,"NA\t");
        }

        l=0;
        for(j=0;j<npops-1;j++) {
            for(k=j+1;k<npops;k++) {
                if(all_Rsb[l][i] != 1234567890) fprintf(results_file,"%f\t",all_Rsb[l][i]);
                else fprintf(results_file,"NA\t");
                l++;
            }
        }
        l=0;
        for(j=0;j<npops-1;j++) {
            for(k=j+1;k<npops;k++) {
                if(all_Rsb[l][i] != 1234567890) fprintf(results_file,"%f\t",(all_Rsb[l][i]-median[l])/sd[l]);
                else fprintf(results_file,"NA\t");
                l++;
            }
        }
        fprintf(results_file,"\n");
    }
    fclose(results_file);
    
    printf("\ndone\n");
    exit(0);
}

int read_row(FILE *plink_file, char *chr_name, double *lox, int **geno, int geno_cols, long int row)
{
    int ch,i;
    int ncol;
    int nfield;
    char field[100];
    
    for(i=0;i<11;i++) {chr_name[i]='\0';}
    for(i=0;i<100;i++) {field[i]='\0';}
    ncol = 0; nfield = 0;
    while((ch = fgetc(plink_file)) != '\n' && ch != EOF && ch != '\r') {
        if(ch != 9 && ch != 32) {
            field[nfield++] = ch;
        }
        else {
            if(nfield) {
                switch(ncol) {
                        case 0:
                        for(i=0;i<nfield;i++)
                        chr_name[i] = field[i];
                        break;
                        case 1:
                        lox[row] = atof(field);
                        break;
                    default:
                        geno[ncol - 2][row] = atoi(field);
                        break;
                }
                ncol++;
                for(i=0;i<nfield;i++) {field[i]='\0';}
            }
            nfield = 0;
            if(ncol-2 > geno_cols) {
                printf("\nError: input file having more columns than defined. row: %ld ncols: %d\n",row,ncol);
                exit(1);
            }
        }
    }
    if(nfield) {
        geno[ncol - 2][row] = atoi(field);
        ncol++;
    }
    if(ncol-2 != geno_cols) {
        printf("\nError: input file having different individuals than defined. row: %ld n_indiv: %d \n",row,ncol);
        exit(1);
    }
    if(ch == EOF)
    return(0);
    return(1);
}

/*compare two double numbers in a long int list*/
int compare_(const void *i,const void *j)
{
    if(*(double *)i < *(double *)j) return -1;
    if(*(double *)i > *(double *)j) return  1;
    return 0;
}

void usage()
{
    printf(TANG_SOFTW);
    printf("\n\nUsage:");
    printf("\nTang_stats [genotype filename (one chrom)] [number of SNPs] [number of indiv] [threshold value (eg=0.1)] [seed (eg=123456)] [number pops] [size pop1] [size pop2] ... [size popN] [name_pop1] ... [name_popN]");
    printf("\n\nOutput file is automatically generated with the input filename plus '_Results_Tang.txt'");
    printf("\n\n");
}

