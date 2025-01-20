#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#define PI 3.141592653589
#define X 81
#define Y 81
#define T 300
#define MD X*Y // MD = X*Y
#define P 1000

double complex psi_0[MD];

double complex A[MD][MD];
double complex M[MD][MD];

double complex b[MD];
double complex sum_GS1[MD];
double complex sum_GS2[MD];

double complex psi_seg[MD];
double complex psi_ant[MD];


double complex psi_t[T][MD];
double complex psi_guess[P][MD];

double V(double x, double y, double L);
double f(double x, double y, double L);


int main(){
    double x[X];
    double y[Y];
    double t[T];

    int t_print = T; //iter temporal que es vol printejar
    double norm[T];

	int j, i, n, k;    // j = x, i = y, n = t
	int p, l; // dummy index
	int p_conv;
	double sum_dist, dist;

    double dt = 0.0001;
    double dx = 0.01;

	double L = (X-1)*dx;
	printf("L : %lf\n", L);
    double sigma = 0.04;
    double kx = 100;
    double ky = 0;
    double x0 = -L/3;
    double y0 = 0;
    int N = X*Y; // MD

    double complex rx = I*dt/(2*dx*dx);
    double complex ry = I*dt/(2*dx*dx);

    // definim vectors x, y, t

    for (j=0; j<X; j++) {
        x[j] = -0.5*L + j*dx;
    }
    for (i=0; i<X; i++) {
        y[i] = -0.5*L + i*dx;
    }
    for (n=0; n<T; n++) {
        t[n] = n*dt;
    }

    // definim psi_0
    for (k = 0; k<MD; k++){
        j = k%(X);
        i = k/(X);
        psi_t[0][k] = 1/(sqrt(2*PI*sigma*sigma)) * exp(-0.25* ( (x[j]-x0)*(x[j]-x0) + (y[i]-y0)*(y[i]-y0) ) /(sigma*sigma)) * cexp( I* (kx*(x[j]-x0)+ky*(y[i]-y0)) );
        }


    // definim CC
    for (k = 0; k<MD; k++) {
        j = k%(X);
        i = k/(X);

        if (i == 0 ) {
            psi_t[0][k] = 0;
        }
        else if (j == 0 ) {
            psi_t[0][k] = 0;
        }
        else if (i == X-1 ) {
            psi_t[0][k] = 0;
        }
        else if (j == X-1) {
            psi_t[0][k] = 0;
        }
    }


    // Crank-Nicolson

    for (p=0; p<MD; p++) { //zero matrix
        for(k=0; k<MD; k++) {
            A[p][k] = 0;
            M[p][k] = 0;
        }
    }


    for (p=0; p<MD; p++) {
        for(k=0; k<MD; k++) {
            j = k%(X);
            i = k/(X);

            if (k == p) {
                A[p][k] = 1 + 2*rx + 2*ry + I*dt/2 * (V(x[j], y[i], L) +f(x[j],y[i], L));
                M[p][k] = 1 - 2*rx - 2*ry - I*dt/2 * (V(x[j], y[i], L) +f(x[j],y[i], L));

                if (k != 0) {
                    A[p][k-1] = -rx;
                    M[p][k-1] = rx;
                }

                if (k != MD-1) {
                    A[p][k+1] = -rx;
                    M[p][k+1] = rx;
                }

                if (k >= X) {
                    A[p][k-X] = -ry;
                    M[p][k-X] = ry;
                }

                if (k < MD-X) {
                    A[p][k+X] = -ry;
                    M[p][k+X] = ry;
                }
            }
        }
    }


    // Gauss-Seidel Ax = b = My

    for (n=1; n<T+1; n++) {

        for (k=0; k<MD; k++) {
            psi_guess[0][k] = psi_t[n-1][k];
            b[k] = 0;

            for (p=0; p<MD; p++) {
                b[k] += M[k][p] * psi_t[n-1][p]; // = (My)_k = b_k
            }

        }

        for (l=1; l<P; l++) {
            for(k=0; k<MD; k++) {

                sum_GS1[k] = 0;
                sum_GS2[k] = 0;

                for (p=1; p<k; p++) {
                    sum_GS1[k] += A[k][p] * psi_guess[l][p];
                }

                for (p=k+1; p<MD-1; p++) {
                    sum_GS2[k] += A[k][p] * psi_guess[l-1][p];
                }

                psi_guess[l][k] = 1/A[k][k] * (b[k] - sum_GS1[k] - sum_GS2[k]);
            }


            sum_dist = 0;
            dist = 0;

            for (k=0; k<MD; k++) {
                sum_dist += pow( (cabs(conj(psi_guess[l][k])*psi_guess[l][k])-cabs(conj(psi_guess[l-1][k])*psi_guess[l-1][k])) ,2);
            }
            dist = sqrt(sum_dist);


            if (dist < 10E-6) {
                printf("%d de %d \n", n, T);
                p_conv = l;
                break;
            }
        }

        for (k=0; k<MD; k++) {
            psi_t[n][k] = psi_guess[p_conv][k];
        }
    }



    // NORMALITZACIO

    int k1, k2, k3, k4;
    double complex sum_norm;

    for (n=0; n< T+1; n++) {
        sum_norm = 0;

        for (j=0; j<X-1; j++) {
            for (i=0; i<Y-1; i++) {
                k = i * (X) + (j);

                // k1 k2
                // k3 k4

                k1 = i * (X) + (j); // = k
                k2 = i * (X) + (j+1);
                k3 = (i+1) * (X) + (j);
                k4 = (i+1) * (X) + (j+1);


                sum_norm += 0.25 *( pow(cabs(psi_t[n][k1]),2) + pow(cabs(psi_t[n][k2]),2) + pow(cabs(psi_t[n][k3]),2) + pow(cabs(psi_t[n][k4]),2) ) ;
            }
        }
    }

    // Simpson (trapezis 2d)

    double sum_1, sum_2, sum_3, sum_4, sum_5, sum_6;
    double sum_norm1;
    double a,b,c,d; // int (a-->b) dx int (c-->d) dy f(x,y)
    int k_ac,k_bc,k_ad,k_bd;

    k_ad = 0;
    k_ac = (X-1) * X;
    k_bc = (X-1) * X + (X-1);
    k_bd = 0 * X + (X-1);

    sum_1 = pow(cabs(psi_t[n][k_ac]),2) + pow(cabs(psi_t[n][k_bc]),2) + pow(cabs(psi_t[n][k_ad]),2) + pow(cabs(psi_t[n][k_bd]),2);

    for (n=0; n< T+1; n++) {

        sum_norm1 = 0;

        sum_2 = 0;
        sum_3 = 0;
        sum_4 = 0;
        sum_5 = 0;
        sum_6 = 0;

        k = i * (X) + (j);

        for (j=0; j<X; j++) {

            sum_4 += pow(cabs(psi_t[n][0 * (X) + (j)]),2);
            sum_5 += pow(cabs(psi_t[n][(X-1) * (X) + (j)]),2);
        }

        for (i=0; i<Y; i++) {

            sum_2 += pow(cabs(psi_t[n][i * (X) + (0)]),2);
            sum_3 += pow(cabs(psi_t[n][i * (X) + (X-1)]),2);
        }

        for (i=1; i<Y-1; i++) {
            for (j=1; j<X-1; j++) {
                k = i * (X) + (j);

                sum_6 += pow(cabs(psi_t[n][k]),2);
            }
        }

        sum_norm1 = (dx*dx)/4 * ( sum_1 + 2*sum_2 + 2*sum_3 + 2*sum_4 + 2*sum_5 + 4*sum_6 );
        norm[n] = sum_norm1;
    }

    // ERROR

    double error[T];

    for (n=0; n<T+1; n++) {
        error[n] = fabs(1-norm[n])/1.0 *100;
        printf("norm Simp %d: %lf \t error: %lf %%\n", n, norm[n], error[n]);

    }



    // Obrim fitxers
    FILE* mides;
    mides = fopen("mides.txt", "w");

    FILE* gaussT;
    gaussT = fopen("gaussT.txt", "w");



    for (n=0; n<T; n++) {
        for (k=0; k<MD; k++) {
            fprintf(gaussT, "%lf\t", cabs(conj(psi_t[n][k]) * psi_t[n][k]));
            if ((k+1)%X == 0) {
                fprintf(gaussT, "\n");
            }
        }
        fprintf(gaussT, "\n\n\n");
    }
    fprintf(mides, "%d \t %d \t %d \t", T, X, Y);

    fclose(mides);

    fclose(gaussT);

    return 0;
}

// POTENCIAL POU INFINIT

/*
double V(double x, double y, double L) {
    double V0;
    V0 = 10E6;


    if (x >= -0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (y >=-0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (x <= -0.5*L) {
        return V0;
    }
    else if (y <= -0.5*L) {
        return V0;
    }
    else {
        return 0;
    }
}
*/


// POTENCIAL HARM�NIC

/*
double V(double x, double y, double L) {
    double m = 1;
    double omega = 1000;
    double V0;
    V0 = 10E6;

    if (x >= -0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (y >=-0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (x <= -0.5*L) {
        return V0;
    }
    else if (y <= -0.5*L) {
        return V0;
    }
    else {
        return 0.5*m*omega*omega*(x*x + y*y);
    }
}
*/

// POTENCIAL PERI�DIC

/*
double V(double x, double y, double L) {
    double V0;
    V0 = 10E6;

    double R=0.1;
    double x_new, y_new;

    if (x >= -0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (y >= -0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (x <= -0.5*L) {
        return V0;
    }
    else if (y <= -0.5*L) {
        return V0;
    }

    //

    if (x >= 0) {
        x_new = x - 2*R* ceil((x-R)/(2*R));
    }
    else {
        x_new = x + 2*R * ceil((-x-R)/(2*R));
    }
    if (y >= 0) {
        y_new = y - 2*R* ceil((y-R)/(2*R));
    }
    else {
        y_new = y + 2*R * ceil((-y-R)/(2*R));
    }

    if ( (pow(x_new,2)+pow(y_new,2)) <= pow(R,2) && pow(x_new,2)+pow(y_new,2) != 0 ) {
        return -1/sqrt(pow(x_new,2)+pow(y_new,2));
    }
    else {
        return V0;
    }
}
*/


// POTENCIAL DOBLE ESCLETXA


double V(double x, double y, double L) {
    double V0 = 10E6;
    double d = 2*0.01; // (gruix)
    double x_esc = L/6; // posici� de la pared escletxa
    double e = 3*0.01; //escletxa

    // escletxa

    if (x >= x_esc && x <= x_esc+d && y >= -L/2 && y <= -L/2+2*L/5) {
        return V0;
    }
    else if (x >= x_esc && x <= x_esc+d && y >= -L/2+2*L/5+e && y <= -L/2+3*L/5-e) {
        return V0;
    }
    else if (x >= x_esc && x <= x_esc+d && y >= -L/2+3*L/5 && y <= L/2) {
        return V0;
    }

    // pou infinit

    else if (x >= -0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (y >=-0.5*L + (X-1)*0.01) {
        return V0;
    }
    else if (x <= -0.5*L) {
        return V0;
    }
    else if (y <= -0.5*L) {
        return V0;
    }
    else {
        return 0;
    }
}


// POTENCIAL ARREGLITO

double f(double x, double y, double L) {

    double value;
    value = pow( (x*x + y*y)/(0.5*L*L), 5);

    return value;
}
