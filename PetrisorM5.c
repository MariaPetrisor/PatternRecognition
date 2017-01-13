#include <stdio.h>
#include <stdlib.h>
#include "mersenne.c"
#include <time.h>
#include <math.h>
#include <ctype.h>

void eroare() //Afiseaza mesaj de eroare in cazul in care alocarea dinamica nu este efectuata corect
{
    printf( "Eroare alocare dinamica a memoriei. \n" );
    exit(1);
}

// functie care simuleaza distributia discreta de probabilitate
int simVarDiscr( double *p )
{
    int k = 0;
    double F = p[0];
    double u = genrand_real3();
    while ( u > F )
    {
        k++;
        F = F + p[k];
    }
    return k;
}

void normal_polar( double *z1, double *z2, double m1, double sigma1 )
{
    double x1, x2, f, u;
    do
    {
        x1 = -1 + 2 * genrand_real3();
        x2 = -1 + 2 * genrand_real3();
        u = x1 * x1 + x2 * x2;
    }while ( u >= 1 );
    f = sqrt( (-2 * log(u)) / u );
    (*z1) = f * x1;
    (*z2) = f * x2;
}

void simulare( double z1, double z2, double *x1, double *x2, double m1, double m2, double sigma1, double sigma2, double ro )
{
    double zz1, zz2;
    normal_polar( &z1, &z2, 0, 1 );
    zz1 = z1;
    normal_polar( &z1, &z2, 0, 1 );
    zz2 = z2;
    (*x1) = sigma1 * zz1 + m1;
    (*x2) = ro * sigma2 * zz1 + (double)sqrt( 1.0 - ro * ro ) * sigma2 * zz2 + m2;
}

void meniu()
{
    printf( "N: normalpolar\n" );
    printf( "V: vector\n" );
    printf( "M: mixtura\n" );
    printf( "\n" );
}

int main()
{
    time_t secunde;
    secunde=time(NULL);
    init_genrand(secunde);

    char c;
    int n = 2020, i, contor = 0, k, K = 4;
    double m[4][2] = { {1.4, -0.3}, {1.3, 1.9}, {-2.3, 0.8}, {-1.1, 0.6} }, sigma[4][2] = { {0.9, 0.3}, {0.65, -0.4}, {-0.52, 0.8}, {0.65, -0.8}};
    double ro[] = { 0.0, 0.64, -0.15, 0.4 }, p[] = { 0.2, 0.4, 0.15, 0.25 };
    double y1, y2, procent, z1, z2, x1, x2;

    FILE *pf, *v1, *v2, *v3, *mix;
    pf = fopen( "normal.txt", "w" );
    if ( pf == NULL )
    {
        printf( "Eroare afisare in fisier.\n" );
        exit(1);
    }
    v1 = fopen( "vector1.txt", "w" );
    if ( v1 == NULL )
    {
        printf( "Eroare afisare in fisier.\n" );
        exit(1);
    }
    v2 = fopen( "vector2.txt", "w" );
    if ( v2 == NULL )
    {
        printf( "Eroare afisare in fisier.\n" );
        exit(1);
    }
    v3 = fopen( "vector3.txt", "w" );
    if ( v3 == NULL )
    {
        printf( "Eroare afisare in fisier.\n" );
        exit(1);
    }
    mix = fopen( "mixtura.txt", "w" );
    if ( mix == NULL )
    {
        printf( "Eroare afisare in fisier.\n" );
        exit(1);
    }

    meniu();
    c = getchar();
    switch( toupper(c) )
        {
            case 'N':
                {
                    for ( i = 0; i < n/2; i++ )
                    {
                        normal_polar( &z1, &z2, m[0][0], sigma[0][0] );
                        y1 = m[0][0] + sigma[0][0] * z1;
                        y2 = m[0][0] + sigma[0][0] * z2;
                        fprintf( pf, "%lf\n%lf\n", y1, y2 );
                        if ( y1 > ( m[0][0] -  2 * sigma[0][0] ) && y1 < ( m[0][0] + 2 * sigma[0][0] ) )
                            contor++;
                        if ( y2 > ( m[0][0] -  2 * sigma[0][0] ) && y2 < ( m[0][0] + 2 * sigma[0][0] ) )
                            contor++;
                    }
                    procent = (double)contor / n;
                    printf( "Procent: %.2lf\n", procent );
                    break;
                }
            case 'V':
                {
                    for ( i = 0; i < n; i++ )
                    {
                        simulare( z1, z2, &x1, &x2, m[0][0], m[0][1], sigma[0][0], sigma[0][1], ro[0] ); // ro nul
                        fprintf( v1, "%lf %lf\n", x1, x2 );
                    }
                    for ( i = 0; i < n; i++ )
                    {
                        simulare( z1, z2, &x1, &x2, m[0][0], m[0][1], sigma[0][0], sigma[0][1], ro[1] ); // ro pozitiv
                        fprintf( v2, "%lf %lf\n", x1, x2 );
                    }
                    for ( i = 0; i < n; i++ )
                    {
                        simulare( z1, z2, &x1, &x2, m[0][0], m[0][1], sigma[0][0], sigma[0][1], ro[2] ); // ro negativ
                        fprintf( v3, "%lf %lf\n", x1, x2 );
                    }
                    break;
                }
            case 'M':
                {
                    for ( i = 0; i < n; i++ )
                    {
                        k = simVarDiscr( p );
                        simulare( z1, z2, &x1, &x2, m[k][0], m[k][1], sigma[k][0], sigma[k][1], ro[k] );
                        fprintf( mix, "%lf %lf\n", x1, x2 );
                    }
                    break;
                }
            default:
                printf( "Eroare. Ati introdus o comanda gresita.\n" );

        }
    fclose(pf);
    fclose(v1);
    fclose(v2);
    fclose(v3);
    fclose(mix);
    return 0;
}
