#include <stdio.h>
#include <cmath>
#include <math.h>
#include <complex>
#include <iostream>
#include <vector>

#define barreras 5
#define anch 100
#define dist 150
#define N (600+anch+(barreras-1)*dist)

using namespace std;
void inicializaphi(complex<double> phi[N+1], double ko);
void calculosgamma(complex<double> gamma[N], complex<double> alpha[N],double V[N+1],double s);
void calculoB(complex<double> b[N], complex<double> phi[N+1], double s );
void calculobeta(complex<double> beta[N], complex<double> gamma[N], complex<double> b[N]);
double calculonorma(complex<double> vector[N+1]);


int main(){
    double s, ko, nciclos, lambda;
    std::complex<double> phi[N+1];
    std::complex<double> alpha[N];
    std::complex<double> beta[N];
    std::complex<double> b[N];
    std::complex<double> Chi[N];
    std::complex<double> gamma[N];
    double V[N+1];
    double normao, norma;
    double inter;
    inter=dist-anch;

    lambda=0.3;

    nciclos=N/5;
    ko=2*M_PI*nciclos/N;
    s=1/(4.0* ko*ko);

    /*
    for (int i = 0; i < N+1; i++)
    {
        if (i>2*N/5 && i<3*N/5)
        {
            V[i]=lambda*ko*ko;
        }
        else
        {
            V[i]=0.0;
        }
        
        
     }
    */
   for (int i = 0; i < 300; i++)
   {
    V[i]=0;
   }
   

    for (int k = 0; k < barreras; k++)
   {
        for (int i = 300+k*anch+k*inter; i < 300+(k+1)*anch+(k+1)*inter; i++)
        {
            if (i>(300+k*anch+k*inter) && i<(300+(k+1)*anch+k*inter))
            {
                V[i]=lambda*ko*ko;
            }
            else
            {
                V[i]=0.0;
            }
        

        }
   }

    //inicializa la funcion en phi
    inicializaphi(phi,ko);
    
    normao=calculonorma(phi);
    for (int i = 0; i < N+1; i++)
    {
        phi[i]=complex(phi[i]);
    }


    //calculo alpha y gamma
    calculosgamma(gamma, alpha, V, s);

    //abro ficherp
    FILE *datosschrodinger;
    datosschrodinger = fopen("datosschrodinger.txt", "w");

    FILE *normaschrodinger;
    normaschrodinger = fopen("norma.txt","w");
    //comienza el ciclo 
    for (int j = 0; j < N*2; j++)
    {
         //rellenamos beta
        calculoB(b,phi,s);
        calculobeta(beta,gamma, b);

        Chi[0]=complex(0.0,0.0);
        for (int i = 0; i < N-1; i++)
        {
            Chi[i+1]=complex(alpha[i]*Chi[i]+beta[i]);
        }
    
        
        for (int i = 0; i < N+1; i++)
        {
            phi[i]=complex(Chi[i]-phi[i]);
            fprintf(datosschrodinger,"%d,%lf,%lf\n", i, norm(phi[i]), V[i]);
        }
        norma=calculonorma(phi);
        fprintf(datosschrodinger,"\n");
        fprintf(normaschrodinger,"%d %lf\n", j, norma/normao);
        printf("%lf\n",norma);
        

    }
     //printf("%lf\n",norma);
   
    

    /*
    for (int i = 0; i < N; i++)
    {
        printf("%lf+i%lf\t", real(phi[i]),imag(phi[i]));
    }
    //*/

    // printf("%lf",s);
    fclose(datosschrodinger);
    fclose(normaschrodinger);



    return 0;
}



void inicializaphi(complex<double> phi[N+1], double ko)
{
    phi[0]=complex(0.0,0.0);
    phi[N]=complex(0.0,0.0);
    std::complex<double> aux1;
    double aux2;
    double Nant=500;

    for (int i = 1; i < N; i++)
    {
        aux1=complex(cos(ko*i), sin(ko*i));
        aux2=exp(-8.0*pow(4.0*i-Nant,2)/(Nant*Nant));
        phi[i]=complex(aux1*aux2);
    }

    return ;
}

void calculosgamma(complex<double> gamma[N], complex<double> alpha[N],double V[N+1],double s)
{
    std::complex<double> a0[N];
    alpha[N-1]=complex(0.0,0.0);
    for (int i = N-1; i > 0; i--)
    {
        a0[i]=complex(-2-V[i],2.0/s);
        gamma[i]=complex(1.0/(a0[i]+alpha[i]));
        alpha[i-1]=complex(-gamma[i]);
    }
    
    return ;
}

void calculoB(complex<double> b[N], complex<double> phi[N+1], double s )
{
    double reales;
    double imagen;
    for (int i = 0; i < N; i++)
    {
        reales=-4*imag(phi[i])/s;
        imagen=4*real(phi[i])/s;
        b[i]=complex(reales,imagen);
        
        

    }

      return ; 
}

void calculobeta(complex<double> beta[N], complex<double> gamma[N], complex<double> b[N])
{
    beta[N-1]=complex(0.0,0.0);
    for (int i = N-1; i > 0; i--)
    {
        beta[i-1]=complex(gamma[i]*(b[i]-beta[i]));
    }
    return ;
}

double calculonorma(complex<double> vector[N+1])
{
    double norma=0.0;
    for (int i = 0; i < N+1; i++)
    {
        norma=norma+norm(vector[i]);
    }
    
    return sqrt(norma);
}

double probderecha(complex<double> vector[N+1])
{
    double suma=0;
    for (int i = N/5; i < N+1; i++)
    {
        suma=suma +abs(vector[i]);
    }
    
    return suma;
}

double probizquierda(complex<double> vector[N+1])
{
    double suma=0;
    for (int i = 0; i < N/5; i++)
    {
        suma=suma +abs(vector[i]);
    }
    
    return suma;
}

