#include <stdio.h>
#include <cmath>
#include <math.h>
#include <complex>
#include <iostream>
#include <vector>

#define N 500
#define T N*10
using namespace std;
void inicializaphi(complex<double> phi[N+1], double ko);
void calculosgamma(complex<double> gamma[N], complex<double> alpha[N],double V[N+1],double s);
void calculoB(complex<double> b[N], complex<double> phi[N+1], double s );
void calculobeta(complex<double> beta[N], complex<double> gamma[N], complex<double> b[N]);
double calculonorma(complex<double> vector[N+1]);
double probderecha(complex<double> vector[N+1]);
void Traspuesta(complex<double> vector[N+1],complex<double> vectort[N+1]);
double ValoresperadoX(complex<double> vector[N+1]);
double desviacion(double vector[T], double media);
void derivadaVecto(complex<double> vector[N+1], complex<double> derivada[N+1]);
double ValoresperadoP(complex<double> derivada[N+1],complex<double> vectort[N+1]);
double ValoresperadoEC(complex<double> derivada2[N+1],complex<double> vectort[N+1]);
double ValoresperadoEP(double V[N+1],complex<double> vector[N+1]);
double ValoresperadoEP(double V[N+1],complex<double> vector[N+1]);



int main(){
    double s, ko, nciclos, lambda, Xmedia, varX, Pmedia, varP, Ec, varEC, Ep, varEp, E, varE;
    std::complex<double> phi[N+1];
    std::complex<double> phitran[N+1];
    std::complex<double> derivada[N+1];
    std::complex<double> derivada2[N+1];

    std::complex<double> alpha[N];
    std::complex<double> beta[N];
    std::complex<double> b[N];
    std::complex<double> Chi[N];
    std::complex<double> gamma[N];
    double V[N+1];
    double xesperado[T];
    double pesperado[T];
    double Ecinetica[T];
    double Epesperado[T];
    double Eesperado[T];
    double normao, norma;
    double Pd, Pdinstan, mt, tmax, p;

    srand(time(NULL));
 //observables
    Xmedia=0;
    Pmedia=0;
    Ec=0;
    Ep=0;
    E=0;


    Pd=0;

    lambda=0.3;

    nciclos=N/5; //N/100;
    ko=2*M_PI*nciclos/N;
    s=1/(4.0* ko*ko);

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
    
    //inicializa la funcion en phi
    inicializaphi(phi,ko);
    

    normao=calculonorma(phi);
    for (int i = 0; i < N+1; i++)
    {
        phi[i]=complex(phi[i]);
    }

    FILE *Xdatos;
    Xdatos = fopen("X.txt", "w");
    FILE *Pdatos;
    Pdatos = fopen("P.txt", "w");
    FILE *ECdatos;
    ECdatos = fopen("EC.txt", "w");
    FILE *Epdatos;
    Epdatos = fopen("Ep.txt", "w");
    FILE *Edatos;
    Edatos = fopen("E.txt", "w");

    //calculo alpha y gamma
    calculosgamma(gamma, alpha, V, s);
    //comienza el ciclo 
    for (int j = 0; j < T; j++)
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
        }
        norma=calculonorma(phi);

        //calculo de X
        xesperado[j]=ValoresperadoX(phi);
        Xmedia=Xmedia+xesperado[j];
        fprintf(Xdatos,"%d\t%lf\n",j,xesperado[j]);
        
        //calculo de P
        Traspuesta(phi, phitran);
        derivadaVecto(phi,derivada);
        pesperado[j]=ValoresperadoP(derivada, phitran);
        Pmedia=Pmedia+pesperado[j];
        fprintf(Pdatos,"%d\t%lf\n",j,pesperado[j]);
        
        //calculo de EC
        derivadaVecto(derivada, derivada2);
        Ecinetica[j]=ValoresperadoEC(derivada2, phitran);
        Ec=Ec+Ecinetica[j];
        fprintf(ECdatos,"%d\t%lf\n",j,Ecinetica[j]);
        
        //calculo Ep
        Epesperado[j]=ValoresperadoEP(V,phi);
        Ep=Ep+Epesperado[j];
        fprintf(Epdatos,"%d\t%lf\n",j,Epesperado[j]);

        //calculo E
        Eesperado[j]=Epesperado[j]+Ecinetica[j];
        E=E+Eesperado[j];
        fprintf(Edatos,"%d\t%lf\n",j,Eesperado[j]);
        


    }
     
    Xmedia=Xmedia/(T);
    varX=desviacion(xesperado,Xmedia);
    printf("%lf\t%lf\n",Xmedia,varX);
    
    Pmedia=Pmedia/(T);
    varP=desviacion(pesperado,Pmedia);
    printf("%lf\t%lf\n",Pmedia,varP);
    
    Ec=Ec/T;
    varEC=desviacion(Ecinetica, Ec);
    printf("%lf\t%lf\n",Ec,varEC);
    
    Ep=Ep/T;
    varEp=desviacion(Epesperado, Ep);
    printf("%lf\t%lf\n",Ep,varEp);

    E=E/T;
    varE=desviacion(Eesperado, E);
    printf("%lf\t%lf\n",E,varE);



    fclose(Xdatos);
    fclose(Pdatos);
    fclose(ECdatos);
    fclose(Epdatos);
    fclose(Edatos);

    return 0;
}



void inicializaphi(complex<double> phi[N+1], double ko)
{
    phi[0]=complex(0.0,0.0);
    phi[N]=complex(0.0,0.0);
    std::complex<double> aux1;
    double aux2;

    for (int i = 1; i < N; i++)
    {
        aux1=complex(cos(ko*i), sin(ko*i));
        aux2=exp(-8.0*pow(4.0*i-N,2)/(N*N));
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
    
    return norma;
}

double probderecha(complex<double> vector[N+1])
{
    double suma=0;
    for (int i = 4*N/5; i < N+1; i++)
    {
        suma=suma +norm(vector[i]);
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

void Traspuesta(complex<double> vector[N+1],complex<double> vectort[N+1])
{
    double preal, pimag;
    for (int i = 0; i < N+1; i++)
    {
        preal=real(vector[i]);
        pimag=-imag(vector[i]);
        vectort[i]=complex(preal,pimag);
    }
    
    return ;
}

double ValoresperadoX(complex<double> vector[N+1])
{
    double suma, norma;
    suma=0;
    norma=0;
    for (int i = 0; i < N+1; i++)
    {
        suma=suma+i*norm(vector[i]);
        norma=norma+norm(vector[i]);
    }
    
    return suma/norma;
}

double desviacion(double vector[T], double media)
{
    double suma=0;
    for (int i = 0; i < 10*N; i++)
    {
        suma=suma+pow((vector[i]-media),2);
    }
    return sqrt(suma/N);

}

void derivadaVecto(complex<double> vector[N+1], complex<double> derivada[N+1])
{
    for (int i = 0; i < N+1; i++)
    {
        if (i==N)
        {
            derivada[i]=0;
        }
        else
        {
            derivada[i]=vector[i+1]-vector[i];
        }
        
        
    }
    return;
}

double ValoresperadoP(complex<double> derivada[N+1],complex<double> vectort[N+1])
{
    double suma=0;
    double norma=0;
    for (int i = 0; i < N+1; i++)
    {
        suma=suma+imag(vectort[i]*derivada[i]);
        norma=norma+norm(vectort[i]);

    }

    return suma/norma;  
} 

double ValoresperadoEC(complex<double> derivada2[N+1],complex<double> vectort[N+1])
{
    double suma=0;
    double norma=0;
    for (int i = 0; i < N+1; i++)
    {
        suma=suma-real(vectort[i]*derivada2[i]);
        norma=norma+norm(vectort[i]);
    }

    return suma/(2*norma); 
}

double ValoresperadoEP(double V[N+1],complex<double> vector[N+1])
{
    double suma=0;
    double norma=0;
    for (int i = 0; i < N+1; i++)
    {
        suma=suma+V[i]*norm(vector[i]);
        norma=norma+norm(vector[i]);
    }
    
    return suma/norma;
}

