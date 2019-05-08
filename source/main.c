//
//  main.c
//  Tang_functions
//
//  Created by Sebastian Ramos-Onsins on 29/01/2019.
//  From an original R script of Eva KF Chan
//  Function to calculate the Rsb statistic for a given chromosome as described in:
//  Tang K, Thornton KR, Stoneking M (2007) A New Approach for Using Genome Scans to Detect Recent Positive Selection in the Human Genome . PLoS Biol 5(7): e171 doi:10.1371/journal.pbio.0050171
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
    long int i;
    int j,k,l;
    int npops;
    int *popsize;
    int popcum;
    
    double *mean;
    double *sd;
    double **all_Rsb;
    long int *nit;

    FILE *plink_file = 0;
    int c;
    long int row;
    char chr_name[11]; //the file must contain the same ID for all positions!
    char file_out[1024];
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
    
    
    if (!(plink_file = fopen(argv[1],"r"))) {
        printf("Error reading the input file %s\n",argv[1]);
        exit(1);
    }
    
    //define variables and arrays:
    geno_rows = L = atol(argv[2]);
    geno_cols = N = atoi(argv[3])/* - 2*/;
    thresh = atof(argv[4]);
    init_seed1(atol(argv[5]));
    npops = atoi(argv[6]);
    popsize = (int *)calloc(npops,sizeof(int));
    for(i=0;i<npops;i++) {
        if(i!=npops-1) popsize[i] = atoi(argv[7+i+1]) - atoi(argv[7+i]);
        else popsize[i] = (N+2) - atoi(argv[7+i]);
    }
    
    //geno is transposed to facilitate the analysis
    geno = (int **)calloc(N,sizeof(int *));
    for(j=0;j<N;j++) {geno[j] = (int *)calloc(L,sizeof(int));}
    lox = (double *)calloc(L,sizeof(double));
    all_iES = (double **)calloc(npops,sizeof(double *));
    for(i=0;i<npops;i++) {all_iES[i] = (double *)calloc(L,sizeof(double));}
    
    mean = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    sd = (double *)calloc(npops*(npops-1)/2,sizeof(double));
    nit = (long int *)calloc(npops*(npops-1)/2,sizeof(long int));

    all_Rsb = (double **)calloc(npops*(npops-1)/2,sizeof(double *));
    for(i=0;i<npops*(npops-1)/2;i++) {all_Rsb[i] = (double *)calloc(L,sizeof(double));}

    //read input data: skip header
    printf("\nReading input file...");
    while((c=getc(plink_file))!='\n' && c!='\r'); //skip header
    row = 0;
    for(row=0;row<L;row++) {
        if((c = read_row(plink_file,chr_name,lox,geno,geno_cols,row))==0) {
            printf("\nError: input file having less rows than defined. nrows: %ld\n",row+1);
            exit(1);
        }
    }
    fclose(plink_file);
    
    //imputation
    printf("\nimputing missing values (previously assigned as genotype=9)...");
    impute_genotypes(geno, &geno_rows, geno_cols);
    L = geno_rows;

    //writing imputed file
    strcat(file_out,argv[1]);
    strcat(file_out,"_imputed.txt\0");
    printf("\nWriting imputed genotype file %s...",file_out);
        //header
    if (!(plink_file = fopen(file_out,"w"))) {
        printf("Error reading the input file %s\n",argv[1]);
        exit(1);
    }
    fprintf(plink_file,"CHR\tPOS\t");
    for(j=0;j<npops;j++) {
        for(k=0;k<popsize[j];k++) {
            fprintf(plink_file,"POP%d_IND%d\t",j+1,k+1);
        }
    }
    fprintf(plink_file,"\n");
        //genotypes
    for(i=0;i<L;i++) {
        fprintf(plink_file,"%s\t%f\t",chr_name,lox[i]);
        for(j=0;j<N;j++) {
            fprintf(plink_file,"%d\t",geno[j][i]);
        }
        fprintf(plink_file,"\n");
    }
    fclose(plink_file);
    
    //converting 2 homozygotes to 0:
    printf("\nConverting all homozygotes to value 0...");
    for(j=0;j<npops;j++) {
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
    
    //calculation Rsb, mean and sd:
    printf("\nComputing lnRsb...");
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
    l=0;
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            mean[l] = mean[l]/(double)nit[l];
            sd[l] = sqrt(sd[l]/(double)nit[l] - mean[l]*mean[l]);
            l++;
        }
    }
    
    //writing output file
    file_out[0] = '\0';
    strcat(file_out,argv[1]);
    strcat(file_out,"_Results_Tang.txt\0");
    printf("\nWriting results in the output file %s...",file_out);

    if (!(plink_file = fopen(file_out,"w"))) {
        printf("Error writing the input file %s\n",file_out);
        exit(1);
    }
        //header
    fprintf(plink_file,"Position\t");
    for(j=0;j<npops;j++) fprintf(plink_file,"iES_POP%d\t",j+1);
    for(j=0;j<npops;j++) fprintf(plink_file,"log(iES_POP%d)\t",j+1);
    for(j=0;j<npops;j++) {
        for(k=j+1;k<npops-1;k++) {
            fprintf(plink_file,"lnRsb(POP%d/POP%d)\t",j+1,k+1);
        }
    }
    for(j=0;j<npops-1;j++) {
        for(k=j+1;k<npops;k++) {
            fprintf(plink_file,"lnRsbN(POP%d/POP%d)\t",j+1,k+1);
        }
    }
    fprintf(plink_file,"\n");
        //data
    for(i=0;i<L;i++) {
        fprintf(plink_file,"%f\t",lox[i]);
        for(j=0;j<npops;j++) {
            if(all_iES[j][i] > 0.0) fprintf(plink_file,"%f\t",all_iES[j][i]);
            else fprintf(plink_file,"NA\t");
        }
        for(j=0;j<npops;j++) {
            if(all_iES[j][i] > 0.0) fprintf(plink_file,"%f\t",log(all_iES[j][i]));
            else fprintf(plink_file,"NA\t");
        }

        l=0;
        for(j=0;j<npops-1;j++) {
            for(k=j+1;k<npops;k++) {
                if(all_Rsb[l][i] != 1234567890) fprintf(plink_file,"%f\t",all_Rsb[l][i]);
                else fprintf(plink_file,"NA\t");
                l++;
            }
        }
        l=0;
        for(j=0;j<npops-1;j++) {
            for(k=j+1;k<npops;k++) {
                if(all_Rsb[l][i] != 1234567890) fprintf(plink_file,"%f\t",(all_Rsb[l][i]-mean[l])/sd[l]);
                else fprintf(plink_file,"NA\t");
                l++;
            }
        }
        fprintf(plink_file,"\n");
    }
    fclose(plink_file);
    
    printf("\ndone\n");
}

int read_row(FILE *plink_file, char *chr_name, double *lox, int **geno, int geno_cols, long int row)
{
    int ch,i;
    int ncol;
    int nfield;
    char field[100];
    
    for(i=0;i<11;i++) {chr_name[i]='\0';}
    ncol = 0; nfield = 0;
    while((ch = fgetc(plink_file)) != '\n' && ch != EOF && ch != '\r') {
        if(ch != 9 && ch != 32) {
            field[nfield++] = ch;
        }
        else {
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
            nfield = 0;
            if(ncol-2 > geno_cols) {
                printf("\nError: input file having more columns than defined. row: %ld ncols: %d\n",row,ncol);
                exit(1);
            }
        }
    }
    if(ncol-2+1 != geno_cols) {
        printf("\nError: input file having different columns than defined. row: %ld ncols: %d\n",row,ncol);
        exit(1);
    }
    if(ch == EOF)
        return(0);
    return(1);
}

void usage()
{
    printf(TANG_SOFTW);
    printf("\n\nUsage:");
    printf("\nTang_stats [Plink filename] [number of rows] [number of individuals] [threshold value] [seed] [number pops] [col pop1] [col pop2] ... [col pop N]");
    printf("\n\nOtput file is automatically generated with the input filename plus '_Results_Tang.txt'");
    printf("\n\n");
}

