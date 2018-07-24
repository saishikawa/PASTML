/**
* Created by azhukova on 1/24/18.
**/

#ifndef PASTML_PASTML_H
#define PASTML_PASTML_H
int runpastml(char *annotation_name, char *tree_name, char *out_annotation_name, char *out_tree_name,
              char *out_parameter_name, char *model, char *prob_method, char *parameter_name, char *out_mp_name);
int is_valid_model(char* model);
int is_valid_prediction_method(char *prob_method);
int is_parsimonious_method(char *prob_method);
int is_marginal_method(char *prob_method);
#endif //PASTML_PASTML_H
