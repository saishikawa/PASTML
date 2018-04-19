//
// Created by azhukova on 2/20/18.
//

#ifndef PARSIMONY_H
#define PARSIMONY_H

#include "pastml.h"

void parsimony(Tree *tree, unsigned long num_annotations, char* method);
void select_parsimonious_states(Tree *tree, size_t num_annotations);
#endif //PARSIMONY_H