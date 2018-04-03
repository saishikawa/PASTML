/**
* Created by azhukova on 1/24/18.
**/

#ifndef PASTML_MAKE_TREE_H
#define PASTML_MAKE_TREE_H

#include "pastml.h"

Tree *read_tree(char *nwk, size_t num_anno);
void free_tree(Tree *tree, size_t num_anno);
int write_nh_tree(Tree *s_tree, char *output_filepath);

#endif //PASTML_MAKE_TREE_H
