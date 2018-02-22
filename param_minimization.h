//
// Created by azhukova on 2/20/18.
//

#ifndef PASTML_PARAM_MINIMIZATION_H
#define PASTML_PARAM_MINIMIZATION_H
double minimize_params(char **tip_array, const int *state_array, int num_tips, int num_annotations, double mu,
                       double *parameters, char **character, char *model_name, double scale_low, double scale_up,
                       double epsilon_low, double epsilon_up);
#endif //PASTML_PARAM_MINIMIZATION_H
