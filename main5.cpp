#define _CRT_SECURE_NO_WARNINGS 
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "winbgi2.h"
#include "rk4.h"

void rhs_fun(double t, double *X, double *F);
void veuler(double t, double*X, double h, int n, void(*rhs_fun)(double, double*, double*), double *X1);
double energia(double* X,double m);


double g=9.81;			//przyspieszenie ziemskie
double l=10;			//dlugosc wahadla


void main()
{
	double omega0=0.1;
	double alfa0=0.1;
	double t;
	double h=0.01;			//krok calkowania
	double FE[2];			//tablice prawych stron rownan
	double FR[2];
	double m;			//masa kulki
	



	graphics(800,600);
	//scale(-0.5,-0.5,0.5,0.5);
	scale(0,0,10,10);

	printf("podaj mase kulki m: \n");
	scanf("%lf", &m);
	

	double XE[2]={omega0,alfa0};			//tablice wartosci zmiennych zaleznych
	double XR[2]={omega0,alfa0};


	for (t=0; t<=10;t+=h)
	{

		veuler(t, XE, h, 2, rhs_fun, FE);
		vrk4(t,XR,h,2,rhs_fun,FR);
		printf("t = %lf\t , energia Euler = %lf\t energia RK4 = %lf\t omega Euler = %lf\t alfa Euler = %lf\t omega RK4 = %lf\t alfa RK4 = %lf\n\n", t, energia(XE,m), energia(XR,m), XE[0], XE[1], XR[0], XR[1]); 

		setcolor(RED);
		//point(XE[1],XE[0]);					//wykres trajektorii ukladu w przestrzeni fazowej (alfa - omega) dla metody Eulera
		//point(XR[1],XR[0]);					//wykres trajektorii ukladu w przestrzeni fazowej (alfa - omega) dla metody RK4
		//point(t,energia(XE,m));			//wykres zaleznosci energii wahadla od t dla metody Eulera
		point(t,energia(XR,m));			//wykres zaleznosci energii wahadla od t dla metody RK4

		XE[0]=FE[0];
		XE[1]=FE[1];
		XR[0]=FR[0];
		XR[1]=FR[1];
		
		
		
	}
	
	
wait();
}

void rhs_fun(double t, double *X, double *F)
{
	double prawa_str[2];
	prawa_str[0]= -(g/l)*sin(X[1]);
	prawa_str[1]= X[0];
	F[0]=prawa_str[0];
	F[1]=prawa_str[1];
}

void veuler(double t, double*X, double h, int n, void(*rhs_fun)(double, double *, double *), double *X1)
{

		double F[2];
		rhs_fun(t, X, F);
		X1[0]=X[0]+h*F[0];
		X1[1]=X[1]+h*F[1];
		
}
double energia(double* X,double m)
{
	return ((m*l*l)/2)*X[0]*X[0]+m*g*l*(1-cos(X[1]));
}