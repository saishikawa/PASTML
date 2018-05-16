#include "pastml.h"
#include "logger.h"
#include "eigen.h"

#define NUM_AA 20
#define SQNUM_AA 400
#define CUNUM_AA 8000
#define NUM_AA_REL_RATES 190

double aaFreq[NUM_AA];
double aaRelativeRate[NUM_AA_REL_RATES];
static double Qij[SQNUM_AA], Cijk[CUNUM_AA], Root[NUM_AA];

static double jttRelativeRates[NUM_AA_REL_RATES] = {
        0.531678, 0.557967, 0.827445, 0.574478, 0.556725, 1.066681, 1.740159, 0.219970, 0.361684, 0.310007, 0.369437,
        0.469395, 0.138293, 1.959599, 3.887095, 4.582565, 0.084329, 0.139492, 2.924161,
        0.451095, 0.154899, 1.019843, 3.021995, 0.318483, 1.359652, 3.210671, 0.239195, 0.372261, 6.529255, 0.431045,
        0.065314, 0.710489, 1.001551, 0.650282, 1.257961, 0.235601, 0.171995,
        5.549530, 0.313311, 0.768834, 0.578115, 0.773313, 4.025778, 0.491003, 0.137289, 2.529517, 0.330720, 0.073481,
        0.121804, 5.057964, 2.351311, 0.027700, 0.700693, 0.164525,
        0.105625, 0.521646, 7.766557, 1.272434, 1.032342, 0.115968, 0.061486, 0.282466, 0.190001, 0.032522, 0.127164,
        0.589268, 0.425159, 0.057466, 0.453952, 0.315261,
        0.091304, 0.053907, 0.546389, 0.724998, 0.150559, 0.164593, 0.049009, 0.409202, 0.678335, 0.123653, 2.155331,
        0.469823, 1.104181, 2.114852, 0.621323,
        3.417706, 0.231294, 5.684080, 0.078270, 0.709004, 2.966732, 0.456901, 0.045683, 1.608126, 0.548807, 0.523825,
        0.172206, 0.254745, 0.179771,
        1.115632, 0.243768, 0.111773, 0.097485, 1.731684, 0.175084, 0.043829, 0.191994, 0.312449, 0.331584, 0.114381,
        0.063452, 0.465271,
        0.201696, 0.053769, 0.069492, 0.269840, 0.130379, 0.050212, 0.208081, 1.874296, 0.316862, 0.544180, 0.052500,
        0.470140,
        0.181788, 0.540571, 0.525096, 0.329660, 0.453428, 1.141961, 0.743458, 0.477355, 0.128193, 5.848400, 0.121827,
        2.335139, 0.202562, 4.831666, 0.777090, 0.098580, 0.405119, 2.553806, 0.134510, 0.303445, 9.533943,
        0.146481, 3.856906, 2.500294, 1.060504, 0.592511, 0.272514, 0.530324, 0.241094, 1.761439,
        0.624581, 0.024521, 0.216345, 0.474478, 0.965641, 0.089134, 0.087904, 0.124066,
        0.436181, 0.164215, 0.285564, 2.114728, 0.201334, 0.189870, 3.038533,
        0.148483, 0.943971, 0.138904, 0.537922, 5.484236, 0.593478,
        2.788406, 1.176961, 0.069965, 0.113850, 0.211561,
        4.777647, 0.310927, 0.628608, 0.408532,
        0.080556, 0.201094, 1.14398,
        0.747889, 0.239697,
        0.165473
};

static double jttFrequencies[NUM_AA] = {
        0.076862, 0.051057, 0.042546, 0.051269, 0.020279, 0.041061, 0.061820, 0.074714, 0.022983, 0.052569, 0.091111,
        0.059498, 0.023414, 0.040530, 0.050532, 0.068225, 0.058518, 0.014336, 0.032303, 0.066374
};

void exchange_params(size_t num_tips, int *states, char **character, char *model) {
    size_t i;

    if (strcmp(model, HKY) == 0) {
        /*put 4 characters in this order : TCAG*/
        for (i = 0; i < num_tips; i++) {
            if (strcmp(character[states[i]], "T") == 0) {
                states[i] = 0;
            } else if (strcmp(character[states[i]], "C") == 0) {
                states[i] = 1;
            } else if (strcmp(character[states[i]], "A") == 0) {
                states[i] = 2;
            } else if (strcmp(character[states[i]], "G") == 0) {
                states[i] = 3;
            } else {
                states[i] == -1;
            }
        }
        strcpy(character[0], "T");
        strcpy(character[1], "C");
        strcpy(character[2], "A");
        strcpy(character[3], "G");
    }

    if (strcmp(model, JTT) == 0) {
        /*put 20 characters in this order : ARNDCQEGHILKMFPSTWYV*/
        for (i = 0; i < num_tips; i++) {
            if (strcmp(character[states[i]], "A") == 0) {
                states[i] = 0;
            } else if (strcmp(character[states[i]], "R") == 0) {
                states[i] = 1;
            } else if (strcmp(character[states[i]], "N") == 0) {
                states[i] = 2;
            } else if (strcmp(character[states[i]], "D") == 0) {
                states[i] = 3;
            } else if (strcmp(character[states[i]], "C") == 0) {
                states[i] = 4;
            } else if (strcmp(character[states[i]], "Q") == 0) {
                states[i] = 5;
            } else if (strcmp(character[states[i]], "E") == 0) {
                states[i] = 6;
            } else if (strcmp(character[states[i]], "G") == 0) {
                states[i] = 7;
            } else if (strcmp(character[states[i]], "H") == 0) {
                states[i] = 8;
            } else if (strcmp(character[states[i]], "I") == 0) {
                states[i] = 9;
            } else if (strcmp(character[states[i]], "L") == 0) {
                states[i] = 10;
            } else if (strcmp(character[states[i]], "K") == 0) {
                states[i] = 11;
            } else if (strcmp(character[states[i]], "M") == 0) {
                states[i] = 12;
            } else if (strcmp(character[states[i]], "F") == 0) {
                states[i] = 13;
            } else if (strcmp(character[states[i]], "P") == 0) {
                states[i] = 14;
            } else if (strcmp(character[states[i]], "S") == 0) {
                states[i] = 15;
            } else if (strcmp(character[states[i]], "T") == 0) {
                states[i] = 16;
            } else if (strcmp(character[states[i]], "W") == 0) {
                states[i] = 17;
            } else if (strcmp(character[states[i]], "Y") == 0) {
                states[i] = 18;
            } else if (strcmp(character[states[i]], "V") == 0) {
                states[i] = 19;
            } else {
                states[i] == -1;
            }
        }
        strcpy(character[0], "A");
        strcpy(character[1], "R");
        strcpy(character[2], "N");
        strcpy(character[3], "D");
        strcpy(character[4], "C");
        strcpy(character[5], "Q");
        strcpy(character[6], "E");
        strcpy(character[7], "G");
        strcpy(character[8], "H");
        strcpy(character[9], "I");
        strcpy(character[10], "L");
        strcpy(character[11], "K");
        strcpy(character[12], "M");
        strcpy(character[13], "F");
        strcpy(character[14], "P");
        strcpy(character[15], "S");
        strcpy(character[16], "T");
        strcpy(character[17], "W");
        strcpy(character[18], "Y");
        strcpy(character[19], "V");
    }
}

void SetRelativeRates(double *inRelativeRate) {
    int i;
    for (i = 0; i < NUM_AA_REL_RATES; i++) {
        aaRelativeRate[i] = inRelativeRate[i];
    }
}

void SetFrequencies(double *inFrequencies) {
    int i;
    for (i = 0; i < NUM_AA; i++) {
        aaFreq[i] = inFrequencies[i];
    }
}

void setJTTFrequencies(double* parameters, size_t n) {
    size_t i;
    for(i = 0; i < n; i++) {
        parameters[i] = jttFrequencies[i];
    }
}

void setHKYFrequencies(double* parameters) {
    parameters[0] = 0.1;
    parameters[1] = 0.4;
    parameters[2] = 0.2;
    parameters[3] = 0.3;
}

void SetupJTTMatrix() {
    int i, j, k;
    double mr;
    double sum;
    double U[SQNUM_AA], V[SQNUM_AA], T1[SQNUM_AA], T2[SQNUM_AA];

    SetRelativeRates(jttRelativeRates);
    SetFrequencies(jttFrequencies);
    k = 0;
    for (i = 0; i < NUM_AA - 1; i++) {
        for (j = i + 1; j < NUM_AA; j++) {
            Qij[i * NUM_AA + j] = Qij[j * NUM_AA + i] = aaRelativeRate[k++];
        }
    }

    for (i = 0; i < NUM_AA; i++) {
        for (j = 0; j < NUM_AA; j++) {
            Qij[i * NUM_AA + j] *= aaFreq[j];
        }
    }

    mr = 0;
    for (i = 0; i < NUM_AA; i++) {
        sum = 0;
        Qij[i * NUM_AA + i] = 0;
        for (j = 0; j < NUM_AA; j++) {
            sum += Qij[i * NUM_AA + j];
        }
        Qij[i * NUM_AA + i] = -sum;
        mr += aaFreq[i] * sum;
    }

    abyx(1.0 / mr, Qij, SQNUM_AA);

    if ((k = eigen(1, Qij, NUM_AA, Root, T1, U, V, T2)) != 0) {
        fprintf(stderr, "\ncomplex roots in SetupAAMatrix");
        exit(0);
    }
    xtoy(U, V, SQNUM_AA);
    matinv(V, NUM_AA, NUM_AA, T1);
    for (i = 0; i < NUM_AA; i++) {
        for (j = 0; j < NUM_AA; j++) {
            for (k = 0; k < NUM_AA; k++) {
                Cijk[i * SQNUM_AA + j * NUM_AA + k] = U[i * NUM_AA + k] * V[k * NUM_AA + j];
            }
        }
    }
}

void SetJTTMatrix(double *matrix, double len) {
    int i, j, k;
    double expt[NUM_AA];
    double *P;

    SetupJTTMatrix();
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
    P = matrix;
    if (len < 1e-6) {
        for (i = 0; i < NUM_AA; i++) {
            for (j = 0; j < NUM_AA; j++) {
                if (i == j)
                    *P = 1.0;
                else
                    *P = 0.0;
                P++;
            }
        }
        return;
    }
    for (k = 1; k < NUM_AA; k++) {
        expt[k] = exp(len * Root[k]);
    }
    for (i = 0; i < NUM_AA; i++) {
        for (j = 0; j < NUM_AA; j++) {
            (*P) = Cijk[i * SQNUM_AA + j * NUM_AA];
            for (k = 1; k < NUM_AA; k++) {
                (*P) += Cijk[i * SQNUM_AA + j * NUM_AA + k] * expt[k];
            }
            P++;
        }
    }
}

void get_pij_hky(const Node *nd, size_t num_frequencies, const double *frequency, double bl) {
    size_t i, j;
    double ts = 8.0, beta, freqTC, freqAG, mul, mul2, rig, lef;

    freqTC = frequency[0] + frequency[1];
    freqAG = frequency[2] + frequency[3];
    beta = 0.5 * 1 / (freqAG * freqTC + ts * (frequency[0] * frequency[1] + frequency[2] * frequency[3]));
    for (i = 0; i < num_frequencies; i++) {
        for (j = 0; j < num_frequencies; j++) {
            //T->T
            if (i == 0 && j == 0) {
                mul = -1.0 * bl * beta;
                lef = frequency[0] * (freqTC + (freqAG) * exp(mul)) / (freqTC);
                mul2 = mul * (1.0 + (freqTC) * (ts - 1.0));
                rig = frequency[1] / freqTC * exp(mul2);
                nd->pij[i][j] = lef + rig;
                //T->C
            } else if (i == 0 && j == 1) {
                mul = -1.0 * bl * beta;
                lef = frequency[1] * (freqTC + (freqAG) * exp(mul)) / (freqTC);
                mul2 = mul * (1.0 + (freqTC) * (ts - 1.0));
                rig = frequency[1] / freqTC * exp(mul2);
                nd->pij[i][j] = lef - rig;
                //T->A
            } else if (i == 0 && j == 2) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[2] * (1 - exp(mul));
                //T->G
            } else if (i == 0 && j == 3) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[3] * (1 - exp(mul));
                //C->T
            } else if (i == 1 && j == 0) {
                mul = -1.0 * bl * beta;
                lef = frequency[0] * (freqTC + (freqAG) * exp(mul)) / (freqTC);
                mul2 = mul * (1.0 + (freqTC) * (ts - 1.0));
                rig = frequency[0] / freqTC * exp(mul2);
                nd->pij[i][j] = lef - rig;
                //C->C
            } else if (i == 1 && j == 1) {
                mul = -1.0 * bl * beta;
                lef = frequency[1] * (freqTC + (freqAG) * exp(mul)) / (freqTC);
                mul2 = mul * (1.0 + (freqTC) * (ts - 1.0));
                rig = frequency[0] / freqTC * exp(mul2);
                nd->pij[i][j] = lef + rig;
                //C->A
            } else if (i == 1 && j == 2) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[2] * (1 - exp(mul));
                //C->G
            } else if (i == 1 && j == 3) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[3] * (1 - exp(mul));
                //A->T
            } else if (i == 2 && j == 0) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[0] * (1 - exp(mul));
                //A->C
            } else if (i == 2 && j == 1) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[1] * (1 - exp(mul));
                //A->A
            } else if (i == 2 && j == 2) {
                mul = -1.0 * bl * beta;
                lef = frequency[2] * (freqAG + (freqTC) * exp(mul)) / (freqAG);
                mul2 = mul * (1.0 + (freqAG) * (ts - 1.0));
                rig = frequency[3] / freqAG * exp(mul2);
                nd->pij[i][j] = lef + rig;
                //A->G
            } else if (i == 2 && j == 3) {
                mul = -1.0 * bl * beta;
                lef = frequency[3] * (freqAG + (freqTC) * exp(mul)) / (freqAG);
                mul2 = mul * (1.0 + (freqAG) * (ts - 1.0));
                rig = frequency[3] / freqAG * exp(mul2);
                nd->pij[i][j] = lef - rig;
                //G->T
            } else if (i == 3 && j == 0) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[0] * (1 - exp(mul));
                //G->C
            } else if (i == 3 && j == 1) {
                mul = -1.0 * bl * beta;
                nd->pij[i][j] = frequency[1] * (1 - exp(mul));
                //G->A
            } else if (i == 3 && j == 2) {
                mul = -1.0 * bl * beta;
                lef = frequency[2] * (freqAG + (freqTC) * exp(mul)) / (freqAG);
                mul2 = mul * (1.0 + (freqAG) * (ts - 1.0));
                rig = frequency[2] / freqAG * exp(mul2);
                nd->pij[i][j] = lef - rig;
                //G->G
            } else if (i == 3 && j == 3) {
                mul = -1.0 * bl * beta;
                lef = frequency[3] * (freqAG + (freqTC) * exp(mul)) / (freqAG);
                mul2 = mul * (1.0 + (freqAG) * (ts - 1.0));
                rig = frequency[2] / freqAG * exp(mul2);
                nd->pij[i][j] = lef + rig;
            }
        }
    }
}
