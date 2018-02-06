//
// Created by azhukova on 1/24/18.
//

#ifndef PASTML_OUTPUT_STATES_H
#define PASTML_OUTPUT_STATES_H

#include "pastml.h"

void output_state_anc_PP(Node *nd, int nb, int nbanno, char **character, FILE *outfile);

void output_state_tip_PP(Node *nd, int nb, int nbanno, char **character, FILE *outfile);

#endif //PASTML_OUTPUT_STATES_H
