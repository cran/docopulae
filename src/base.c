#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>


#define DLESS(d) (d < ((double)0))
#define DZERO(d) (d == ((double)0))


void rowmatch_double(double* x, int* _xrn, int* _xcn,
                     double* tab, int* _tabrn,
                     int* r) {
    /*
     * in R:
     * - check shapes
     * - rowsort x, tab
     */
    int xrn = *_xrn;
    int xcn = *_xcn;
    int tabrn = *_tabrn;
    int i1, i2, j, j1, j2;
    double d;
    bool nextX = false;
    bool nextTab = false;

    i2 = 0;
    for (i1 = 0; i1 < xrn; i1++) { // for each row in x
        r[i1] = R_NaInt;
        for (; i2 < tabrn; i2++) { // find row in tab
            j1 = i1;
            j2 = i2;
            for (j = 0; j < xcn; j++) { // for each column
                d = x[j1] - tab[j2];

                if (DLESS(d)) {
                    nextX = true;
                    break;
                } else if (!DZERO(d)) {
                    nextTab = true;
                    break;
                }

                j1 += xrn;
                j2 += tabrn;
            }

            if (nextX) {
                nextX = false;
                break;
            }

            if (nextTab) {
                nextTab = false;
                continue;
            }

            r[i1] = i2;
            break;
        }
    }
}


void rowsduplicated_double(double* x, int* _nrow, int* _ncol, Rboolean* r) {
    int nrow = *_nrow;
    int ncol = *_ncol;
    int i, j, k;
    double d;
    Rboolean duplicate;

    if (nrow != 0) {
        r[0] = FALSE;
    }

    for (i = 1; i < nrow; i++) {
        k = i;
        duplicate = TRUE;

        for (j = 0; j < ncol; j++) {
            d = x[k] - x[k - 1];
            if (!DZERO(d)) {
                duplicate = FALSE;
                break;
            }

            k += nrow;
        }

        r[i] = duplicate;
    }
}
