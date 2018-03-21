//
// Created by azhukova on 2/20/18.
//

#ifndef PASTML_PARAM_MINIMIZATION_H
#define PASTML_PARAM_MINIMIZATION_H
double minimize_params(Tree* s_tree, size_t num_annotations, double *parameters, char **character, char *model,
                       double scale_low, double scale_up, double epsilon_low, double epsilon_up);
#endif //PASTML_PARAM_MINIMIZATION_H
