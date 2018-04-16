/**
* Created by azhukova on 1/24/18.
**/

#ifndef PASTML_MARGINAL_APPROXI_H
#define PASTML_MARGINAL_APPROXI_H

void choose_likely_states(Tree *tree, size_t n);
void normalize_result_probabilities(Tree *tree, size_t n);
void set_id_best_states(Tree *tree, size_t n);
void choose_best_marginal_states(Tree *tree, size_t n);

#endif //PASTML_MARGINAL_APPROXI_H
