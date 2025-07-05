#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// 数据结构定义
#define MAX_NODES 1000
#define dt 0.1

int main()
{

    // const int N = MAX_NODES; // 假设最大节点数为 2000
    srand(time(NULL));
    // const double beta = 0.04;
    const double miu = 0.1;
    const double beta2 = 0.0 * miu / 4.38;
    const double w1 = 0.8;
    // double rx=0;

    srand(time(NULL)); // 设置随机种子

    for (size_t z = 0; z < 61; z++)
    {

        double beta1 = z * 0.4 * miu / 20.688;
        // double beta1 = z*0.0002;

        double dI = 0, I = 0.6, S = 1 - I, k = 20.688, k2 = 4.38, dS = 0;
        double dSS = 0, dII = 0, dSI = 0;
        double SS = S * S, SI = S * I, II = I * I;
        double fsi = k * k * 998 / (999 * 999) / (k - 1) + 2 * k2 / k / (k - 1);
        double IISS0 = SI * SI * II * SS / (S * S * I * I);
        double IISS1 = SI * SI * SI * SS * II / (S * S * S * I * I * I);
        double IISS2 = SS * pow(SI, 4) * II / (pow(S, 4) * pow(I, 4));
        double IISI0 = pow(SI, 3) * II / (S * S * I * I);
        double IISI1 = pow(SI, 3) * II * II / (S * S * pow(I, 4));
        double IISI2 = pow(SI, 3) * pow(II, 3) / (S * S * pow(I, 6));
        double IISS = pow((1 - fsi), 2) * IISS0 + 2 * (1 - fsi) * fsi * IISS1 + fsi * fsi * IISS2;
        double IISI = pow((1 - fsi), 2) * IISI0 + 2 * (1 - fsi) * fsi * IISI1 + fsi * fsi * IISI2;
        double SSI = (1 - fsi) * SS * SI / S + fsi * SS * SI * SI / (S * S * I);
        double ISI = (1 - fsi) * SI * SI / S + fsi * II * SI * SI / (S * I * I);
        double ISID = I * S * I;
        IISS = SI * SI * II * SS / (S * S * I * I);
        IISI = pow(SI, 3) * II / (S * S * I * I);
        SSI = SS * SI / S;
        ISI = SI * SI / S;
        ISID = I * S * I;
        //printf("%lf %lf %lf %lf %lf %lf %lf\n", SS, SI*k,II,IISS,IISI,SSI,ISI);
       
        for (size_t c = 0; c < 2000000; c++)
        {

            dI = -miu * I + beta1 * k * SI + beta2 * k2 * ISID;
            dS = miu * I - beta1 * k * SI - beta2 * k2 * ISID;
            dSS = 2 * (miu + w1) * SI - 2 * beta1 * (k - 1) * SSI - 2 * beta2 * k2 / k * (k - 2) * IISS;
            // dSS = 2 * (1 - (1 - miu) * (1 - w1)) * SI - 2 * beta1 * (k - 1) * SSI - 2 * beta2 * k2 / k * (k - 2) * IISS;
            // dSI = miu * II - (1 - (1 - miu) * (1 - w1)) * SI + beta1 * (k - 1) * SSI - beta1 * (k - 1) * ISI + beta2 * k2 / k * (k - 2) * IISS - 2 * beta2 * k2 / k * ISID - beta2 * k2 / k * (k - 2) * IISI - beta1 * SI;
            dSI = miu * II - (miu + w1) * SI + beta1 * (k - 1) * SSI - beta1 * (k - 1) * ISI + beta2 * k2 / k * (k - 2) * IISS - 2 * beta2 * k2 / k * ISID - beta2 * k2 / k * (k - 2) * IISI - beta1 * SI;
            dII = -2 * miu * II + 2 * beta1 * (k - 1) * ISI + 2 * beta1 * SI + 4 * beta2 * k2 / k * ISID + 2 * beta2 * k2 / k * (k - 2) * IISI;
            // dSS = 2 * (miu) * SI + w1* SI - 2 * beta1 * (k - 1) * SSI - 2 * beta2 * k2 / k * (k - 2) * IISS;
            // dSI = miu * II - (miu ) * SI- w1* SI + beta1 * (k - 1) * SSI - beta1 * (k - 1) * ISI + beta2 * k2 / k * (k - 2) * IISS - 2 * beta2 * k2 / k * ISID - beta2 * k2 / k * (k - 2) * IISI - beta1 * SI;

            I = I + dI * dt;
            S = S + dS * dt;
            double sun = S + I;
            S = S / sun;
            I = I / sun;
            // S = 1-I;

            SS = SS + dSS * dt;
            II = II + dII * dt;
            // SI = (1 - SS - II) / 2.0;
            SI = SI + dSI * dt;
            //SS = 1-2*SI-II;
            double SUM = SS + II + 2 * SI;
            SS = SS / SUM;
            II = II / SUM;
            SI = SI / SUM;
            //SI = (1 - SS - II) / 2.0;
            //printf("%lf %lf %lf %lf %lf %.16lf\n", beta1 * 20.82 / miu, I, dI*pow(10,9), dSS, dII,dI*dt);
            //break;
            //  printf("%lf %lf %lf\n",I*I,II,I);
            if (fabs(dII*pow(10,14)) <= 1)
            {
                break;
            }
            else
            {

                IISS0 = SI * SI * II * SS / (S * S * I * I);
                IISS1 = SI * SI * SI * SS * II / (S * S * S * I * I * I);
                IISS2 = SS * pow(SI, 4) * II / (pow(S, 4) * pow(I, 4));
                IISI0 = pow(SI, 3) * II / (S * S * I * I);
                IISI1 = pow(SI, 3) * II * II / (S * S * pow(I, 4));
                IISI2 = pow(SI, 3) * pow(II, 3) / (S * S * pow(I, 6));
                IISS = pow((1 - fsi), 2) * IISS0 + 2 * (1 - fsi) * fsi * IISS1 + fsi * fsi * IISS2;
                IISI = pow((1 - fsi), 2) * IISI0 + 2 * (1 - fsi) * fsi * IISI1 + fsi * fsi * IISI2;
                SSI = (1 - fsi) * SS * SI / S + fsi * SS * SI * SI / (S * S * I);
                ISI = (1 - fsi) * SI * SI / S + fsi * II * SI * SI / (S * I * I);
                ISID = I * S * I;

                IISS = SI * SI * II * SS / (S * S * I * I);
                IISI = pow(SI, 3) * II / (S * S * I * I);
                SSI = SS * SI / S;
                ISI = SI * SI / S;
                
                ISID = I * S * I;
                
            }
            // printf("ss%lf %lf %lf %lf %lf %lf\n", IISS0, IISS1, IISS2, IISI0, IISI1, IISI2);
            // printf("aa%lf %lf %lf %lf %lf %lf\n", IISS, IISI, SSI, ISI, SS * SI / S, (SS / S) * (SI / S) * (SI / I));
            // printf("cc%lf %lf %lf %lf %lf %lf\n", (SI / I), (II / I), (1 / S), (SI / S), ISID, fsi);

            // printf("%lf %lf %lf %lf %lf +%lf %lf %lf+%lf %lf %lf+%lf %lf\n", beta1 * 20.82 / miu, I, dI, dSS, dII, IISS0, IISS1, IISS2, IISI0, IISI1, IISI2, IISS, IISI);
            //   printf("%lf %lf %lf %lf %lf=%lf\n",IISS2,pow((SI / S),4),(II / I),(SS / I),(SI / S),I);
            //printf("=%lf %d %lf\n", I, c,dI);
        }

        // printf("%lf %lf\n", beta1 * 20.78 / miu, I);
        //  if(I==0)
        //  {
        //      printf("%lf %lf %lf %lf %lf\n", I,I,I,I,I);
        //  }
        //  else
        //  {
        //      printf("%lf %lf %lf %lf %lf\n", I-((double)rand() / RAND_MAX) * 0.03,I-((double)rand() / RAND_MAX) * 0.03,I-((double)rand() / RAND_MAX) * 0.03,I-((double)rand() / RAND_MAX) * 0.03,I-((double)rand() / RAND_MAX) * 0.03);
        //  }

        if (I >= 0 && I <= 1)
        {
            printf("%lf %lf\n", z * 0.4, I);
        }
        else
        {
            printf("%lf %lf\n", z * 0.4, 0);
        }

        //break;
    }

    system("pause");
    return 0;
}