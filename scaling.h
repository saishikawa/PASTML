//
// Created by azhukova on 2/22/18.
//

#ifndef PASTML_SCALING_H
#define PASTML_SCALING_H
int upscale_node_probs(double* array, int n);
void rescale(double *array, int i, int curr_scaler_pow);
int get_scaling_pow(double value);
#endif //PASTML_SCALING_H
