#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define Nt 300  //broj koraka unutar bloka
#define Nw 500  //broj šetača
#define Nb 200  //ukupan broj blokova
#define NbSkip 10
#define V0 10

double a;

double Psi(double x)
{   //prima veijednost x i vraća vrijednost probne valne funkcije
    if ((x > 0.0) && (x < 1.0)) {   //provjera jeli čest unutar jame
        return x * (x - 1.0) * exp(-a * x);  //ako je unutar jame vraća funkciju
    }
    else {
        return 0.0;
    }
}
//za određene položaje x vraća valnu funkciju
double lokalna_energija(double x)
{
    return -0.5 * a * a
           - (1.0 + a - 2.0 * a * x) / (x * (x - 1.0))
           + V0 * x;
}
//glavna funkcija koja izvodi jednu VMC simulaciju za trenutačnu vrijednost parametra a
void VMC(long *idum, double *Eprosjek, double *greska, int spremanje)
{  //Eprosjek - adresa za zapis konačne srednje energije
    int ib, it, iw, n = 0;  //cijeli brojevi

    //trenutni položaj čestica, gustoća vjerojatnosti za svaki položaj, 
    //lokalna energija svakog šetaća
    double x[Nw], P[Nw], E[Nw];   //definiraju se 3 polja
    double dx_max = 0.5;  //max pomak šetača

    double suma = 0.0;   //suma blokova
    double suma2 = 0.0;  //suma kvadrata blokova(energije)

    double accept = 0.0;   //ukupan broj prihvaćenih Metropolisovih pomaka
    double acc_ib;

    FILE *f = NULL;

    if (spremanje) {
        f = fopen("Energija.dat", "w");

        if (f == NULL) {
            printf("Greska pri otvaranju Energija.dat.\n");
            exit(EXIT_FAILURE);
        }
        //upisuje zaglavlje datoteke
        fprintf(f, "# blok E_bloka tekuci_prosjek\n");
    }

    // Inicijalizacija setača 
    for (iw = 0; iw < Nw; iw++) {
        x[iw] = ran1(idum);
        P[iw] = Psi(x[iw]) * Psi(x[iw]);
        E[iw] = lokalna_energija(x[iw]);
    }

    // Petlja po blokovima
    for (ib = 1; ib <= Nb; ib++) {
        //na pocetku svakjog bloka energ je postavljena 0
        double Eblok = 0.0; 

        for (it = 0; it < Nt; it++) {
            //Na početku svakog MC koraka postavlja sumu energija šetača na nulu
            double Esrednja = 0.0;
            //U svakom MC koraku pokušava se pomaknuti svaki od 500 šetača
            for (iw = 0; iw < Nw; iw++) {
                //Generira probni položaj xp
                double xp = x[iw]
                            + 2.0 * (ran1(idum) - 0.5) * dx_max;

                double Pp = Psi(xp) * Psi(xp);

                //uspoređuje koliko je probni položaj vjerojatan u odnosu na trenutačni položaj
                double T = P[iw] > 0.0
                           ? Pp / P[iw]
                           : 0.0;

                if (T > 1.0)
                    T = 1.0;

                if (ran1(idum) <= T) {
                    x[iw] = xp;
                    P[iw] = Pp;
                    E[iw] = lokalna_energija(xp);
                    //Povećava ukupan broj prihvaćenih pomaka za 1
                    accept++;
                }
                //Dodaje lokalnu energiju trenutačnog šetača u sumu
                Esrednja += E[iw];
            }
            //srednja energija svih šetača u tom koraku
            Eblok += Esrednja / Nw;
        }
        //dodaje u sumu energija bloka i dijwli s brojem koraka
        Eblok /= Nt; //srednja energija jednog bloka

        // Ukupni udio prihvacenih pomaka do trenutnog bloka 
        acc_ib = accept / ((double)ib * Nw * Nt);

         //dx_max se prilagodava tijekom cijele simulacije
         //kako bi prihvacanje bilo oko 50 %
         //Kad su pomaci mali, novi položaj je vrlo blizu starom položaju i često se prihvaća
        if (acc_ib > 0.5)
            dx_max *= 1.05;
        else if (acc_ib < 0.5)
            dx_max *= 0.95;
 
        //Prvih NbSkip blokova sluzi za stabilizaciju.
        //Oni se ne koriste u racunanju energije i pogreske. 
        if (ib > NbSkip) {

            double tekuci;
            //Povećava broj korisnih blokova
            n++;
            //energiju trenutačnog bloka
            suma += Eblok;
            suma2 += Eblok * Eblok;

            //srednju energiju svih do sada prihvaćenih blokova -tekuća
            tekuci = suma / n;

            /*
             * Zapisuju se samo blokovi nakon stabilizacije.
             * Prvi stupac je efektivni broj bloka:
             * 1, 2, ..., Nb - NbSkip.
             */
            if (spremanje) {
                fprintf(f, "%d %.10f %.10f\n",
                        n, Eblok, tekuci);
            }
        }
    }
    //konačni prosjek energija blokova
    *Eprosjek = suma / n;

    double srednji_kvadrat;
    double varijanca;

    srednji_kvadrat = suma2 / n;
    varijanca = srednji_kvadrat - pow(*Eprosjek, 2);

    *greska = sqrt(varijanca / (n - 1));

    if (f != NULL)
        fclose(f);
}

int main(void)
{
    long idum = -212;

    double E, greska;
    double Emin = 1.0e30;
    double amin = 0.0;
    double greska_min = 0.0;

    FILE *f = fopen("E_vs_a.dat", "w");

    if (f == NULL) {
        printf("Greska pri otvaranju E_vs_a.dat.\n");
        return EXIT_FAILURE;
    }

    fprintf(f, "# a energija pogreska\n");

    for (a = -2.0; a <= 2.0001; a += 0.1) {

        VMC(&idum, &E, &greska, 0);

        fprintf(f, "%.2f %.10f %.10f\n",
                a, E, greska);

        printf("a = %.2f, E = %.8f +- %.8f\n",
               a, E, greska);

        if (E < Emin) {
            Emin = E;
            amin = a;
            greska_min = greska;
        }
    }

    fclose(f);

    printf("\nMinimum: a = %.2f, E = %.8f +- %.8f\n",
           amin, Emin, greska_min);

    a = amin;

    VMC(&idum, &E, &greska, 1);

    printf("Zavrsna procjena: E = %.8f +- %.8f\n",
           E, greska);

    return EXIT_SUCCESS;
}