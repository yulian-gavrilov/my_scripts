#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define total_CA 64
//int total_CA = 64;


struct aframe
{
    float *x;
    float *y;
    float *z;
};


int main()
{
    clock_t begin = clock();

    //int total_CA = 64;

    FILE *fp = NULL;
    int countCA = 0;
    int countFrames = 0;
    float x_test,y_test,z_test;
    char *filename = "md1us_all_no200ns_wt_disp_CA_fit_dt100.txt";
    char *filename_out = "md1us_all_no200ns_wt_disp_CA_fit_dt100_DM.txt";
    char *filename_out_sd = "md1us_all_no200ns_wt_disp_CA_fit_dt100_SDM.txt";

    int nframes = 2400;//100
    
    //struct aframe frames_vec[nframes];
    struct aframe frames;
    struct aframe *pframes;

    pframes = &frames;

    // allocate memory for vectors

    pframes->x = (float*)malloc(total_CA * sizeof(float));
    pframes->y = (float*)malloc(total_CA * sizeof(float));
    pframes->z = (float*)malloc(total_CA * sizeof(float));


    if(pframes == NULL) // check that the memory was allocated successfully
        exit(-1);


    float DM[total_CA][total_CA];
    float DM_sum_avg[total_CA][total_CA];
    float forSDM[total_CA][total_CA];
    float forSDMsum[total_CA][total_CA];
    float SDM[total_CA][total_CA];


    fp = fopen(filename, "r"); 
    if(fp == NULL)
    {
        printf("Failed to open %s.\n", filename);
        return(-1);
    }

    while (countFrames < nframes)
    {
        while (countCA < total_CA)
        {

            fscanf(fp, "%f %f %f", &pframes->x[countCA], &pframes->y[countCA], &pframes->z[countCA]);
            //printf ("%f %f %f\n",frames_vec[countFrames].x[countCA], frames_vec[countFrames].y[countCA], frames_vec[countFrames].z[countCA]);

            //fscanf(fp, "%f %f %f", &frames.x[countCA], &frames.y[countCA], &frames.z[countCA]);

            ++countCA;
        }

        countCA = 0;
        countFrames+=1;

        // now we have coordinates of one frame

        for (int i = 0;i<total_CA;++i)
        {
            for (int j = 0;j<total_CA;++j)
        {
                DM[i][j] = pow(pow(frames.x[i]-frames.x[j],2) + pow(frames.y[i]-frames.y[j],2) + pow(frames.z[i]-frames.z[j],2),0.5);
                //printf("%f ",DM[i][j]);
                DM_sum_avg[i][j] += DM[i][j]; 
            }
            //printf("\n CA %d\n",i);
        }
        // empty DM arry. Probably not required.
        memset(DM, 0, sizeof(DM));

    }

    memset(DM, 0, sizeof(DM));
    fclose(fp); // close the file
    fp = NULL;



    fp = fopen(filename_out, "w"); 
    if(fp == NULL)
    {
        printf("Failed to open %s.\n", filename_out);
        return(-1);
    }

    for(int i = 0;i<total_CA;++i) 
    {
        for(int j = 0;j<total_CA;++j)
        {
            DM_sum_avg[i][j] /= nframes;
            fprintf(fp, "%f ", DM_sum_avg[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp); // close the file
    fp = NULL;

    //////////////////////////////////////

    fp = fopen(filename, "r"); 
    if(fp == NULL)

    {
        printf("Failed to open %s.\n", filename);
        return(-1);
    }

    countFrames = 0;
    countCA=0;
    while (countFrames < nframes)
    {
        while (countCA < total_CA)
        {
            //fscanf(fp, "%f %f %f", &frames.x[countCA], &frames.y[countCA], &frames.z[countCA]);
            fscanf(fp, "%f %f %f", &pframes->x[countCA], &pframes->y[countCA], &pframes->z[countCA]);
            ++countCA;
        }

        countCA = 0;
        countFrames+=1;

        // now we have coordinates of one frame

        for (int i = 0;i<total_CA;++i)
        {
            for (int j = 0;j<total_CA;++j)
            {
                DM[i][j] = pow(pow(frames.x[i]-frames.x[j],2) + pow(frames.y[i]-frames.y[j],2) + pow(frames.z[i]-frames.z[j],2),0.5);
                forSDM[i][j] = pow(DM[i][j]-DM_sum_avg[i][j],2);
                forSDMsum[i][j] += forSDM[i][j]; 
                //printf("DM[i][j] %f SDM[i][j] %f\n",DM[i][j], forSDM[i][j]);

            }
        }
        // empty DM arry. Probably not required.
        memset(DM, 0, sizeof(DM));

    }

    fclose(fp); // close the file
    fp = NULL;



    fp = fopen(filename_out_sd, "w"); 
    if(fp == NULL)
    {
        printf("Failed to open %s.\n", filename_out_sd);
        return(-1);
    }

    for(int i = 0;i<total_CA;++i) 
    {
        for(int j = 0;j<total_CA;++j)
        {
            SDM[i][j] = pow(forSDMsum[i][j]/nframes,0.5);
            //SDM[i][j] = forSDM[i][j];

            fprintf(fp, "%f ", SDM[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp); // close the file
    fp = NULL;
    /////////////////////////////////////

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time_spent = %f s\n",time_spent);


    free(pframes->x);
    free(pframes->y);
    free(pframes->z);

    return 0;
}

