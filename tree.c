#include <errno.h>
#include <stdbool.h>
#include "pastml.h"
#include "logger.h"


bool isTip(const Node *nd) {
    return nd->nb_neigh == 1;
}

int index_toplevel_colon(const char *in_str, int begin, int end) {
    /* returns the index of the (first) toplevel colon only, -1 if not found */
    int level = 0, i;
    for (i = end; i >= begin; i--) {/* more efficient to proceed from the end in this case */
        switch (in_str[i]) {
            case ')':
                level++;
                break;
            case '(':
                level--;
                break;
            case ':':
                if (level == 0) return i;
        } /* endswitch */
    } /* endfor */
    return -1;
} /* end index_toplevel_colon */

int parse_double(char *in_str, int begin, int end, double *location) {
    /* this function parses a numerical value and puts it into location. Meant to be used for branch lengths. */
    if (end < begin) {
        fprintf(stderr, "Missing branch length at offset %d in the New Hampshire string. Branch length set to 0.\n",
                begin);
        sscanf("0.0", "%lg", location);
        return EXIT_FAILURE;
    }
    char numerical_string[52] = {'\0'};
    strncpy(numerical_string, in_str + begin, end - begin + 1);
    int n_matches = sscanf(numerical_string, "%lg", location);
    if (n_matches != 1) {
        fprintf(stderr, "Fatal error in parse_double: unable to parse a number out of \"%s\". Aborting.\n",
                numerical_string);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
} /* end parse_double */

int strip_toplevel_parentheses(const char *in_str, int begin, int end, int *pair) {
    /* returns the new (begin,end) pair comprising all chars found strictly inside the toplevel parentheses.
       The input "pair" is an array of two integers, we are passing the output values through it.
       It is intended that here, in_str[pair[0]-1] == '(' and in_str[pair[1]+1] == ')'.
       In case no matching parentheses are simply return begin and end in pair[0] and pair[1]. It is NOT an error. */
    /* This function also tests the correctness of the NH syntax: if no balanced pars, then return an error and abort. */
    int i, found_par = 0;

    pair[0] = end + 1;
    pair[1] = -1; /* to ensure termination if no parentheses are found */

    /* first seach opening par from the beginning of the string */
    for (i = begin; i <= end; i++)
        if (in_str[i] == '(') {
            pair[0] = i + 1;
            found_par += 1;
            break;
        }

    /* and then search the closing par from the end of the string */
    for (i = end; i >= begin; i--)
        if (in_str[i] == ')') {
            pair[1] = i - 1;
            found_par += 1;
            break;
        }

    switch (found_par) {
        case 0:
            pair[0] = begin;
            pair[1] = end;
            break;
        case 1:
            fprintf(stderr,
                    "Syntax error in NH tree: unbalanced parentheses between string indices %d and %d. Aborting.\n",
                    begin, end);
            return EXIT_FAILURE;
    } /* end of switch: nothing to do in case 2 (as pair[0] and pair[1] correctly set), and found_par can never be > 2 */
    return EXIT_SUCCESS;
}

int count_outer_commas(const char *in_str, int begin, int end) {
    /* returns the number of toplevel commas found, from position begin included, up to position end. */
    int count = 0, level = 0, i;
    for (i = begin; i <= end; i++) {
        switch (in_str[i]) {
            case '(':
                level++;
                break;
            case ')':
                level--;
                break;
            case ',':
                if (level == 0) count++;
        } /* endswitch */
    } /* endfor */
    return count;
} /* end count_outer_commas */

int process_name_and_brlen(Node *son_node, char *in_str, int begin, int end) {
    /* looks into in_str[begin..end] for the branch length of the "father" edge
       and updates the edge and node structures accordingly */
    int colon = index_toplevel_colon(in_str, begin, end);
    int closing_par = -1, opening_bracket = -1;
    int i, ignore_mode, name_begin, name_end, name_length, effective_length;
    double brlen = .0;

    /* processing the optional BRANCH LENGTH... */
    if (colon == -1) {
        son_node->branch_len = 0.0;
    } else {
        if (EXIT_SUCCESS != parse_double(in_str, colon + 1, end, &brlen)) {
            return EXIT_FAILURE;
        }
        son_node->branch_len = brlen;
    }

    /* then scan backwards from the colon (or from the end if no branch length) to get the NODE NAME,
       not going further than the first closing par */
    /* we ignore the NHX-style comments for the moment, hence the detection of the brackets, which can contain anything but nested brackets */
    ignore_mode = 0;
    for (i = (colon == -1 ? end : colon - 1); i >= begin; i--) {
        if (in_str[i] == ']' && ignore_mode == 0) { ignore_mode = 1; }
        else if (in_str[i] == ')' && ignore_mode == 0) {
            closing_par = i;
            break;
        }
        else if (in_str[i] == '[' && ignore_mode) {
            ignore_mode = 0;
            opening_bracket = i;
        }
    } /* endfor */

    name_begin = (closing_par == -1 ? begin : closing_par + 1);
    name_end = (opening_bracket != -1) ? (opening_bracket - 1): (colon == -1 ? end : colon - 1);
    /* but now if the name starts and ends with single or double quotes, remove them */
    if (in_str[name_begin] == in_str[name_end] && (in_str[name_begin] == '"' || in_str[name_begin] == '\'')) {
        name_begin++;
        name_end--;
    }
    name_length = name_end - name_begin + 1;
    effective_length = (name_length > MAX_NAMELENGTH ? MAX_NAMELENGTH : name_length);
    son_node->name = (char *) malloc((MAX_NAMELENGTH) * sizeof(char));
    if (name_length >= 1) {
        strncpy(son_node->name, in_str + name_begin, (size_t) effective_length);
        son_node->name[effective_length] = '\0'; /* terminating the string */
    }
    return EXIT_SUCCESS;
} /* end of process_name_and_brlen */


Node *create_son_and_connect_to_father(Node *current_node, Tree *current_tree, int direction, char *in_str, int begin,
                                       int end) {
    /* This function creates (allocates) the son node in the given direction from the current node.
       It also creates a new branch to connect the son to the father.
       The array structures in the tree (nodes) are updated accordingly.
       Branch length and node name are processed.
       The input string given between the begin and end indices (included) is of the type:
       (...)node_name:length
       OR
       leaf_name:length
       OR
       a:1,b:0.31,c:1.03
       In both cases the length is optional, and replaced by MIN_BR_LENGTH if absent.
       In case of a failure the function returns NULL*/

    if (direction < 0) {
        fprintf(stderr, "Error in the direction given to create a son! Aborting.\n");
        return NULL;
    }

    Node *son = (Node *) malloc(sizeof(Node));
    son->id = current_tree->next_avail_node_id++;
    current_tree->nodes[son->id] = son;
    current_tree->nb_nodes++;

    son->name = NULL;

    current_tree->nb_edges++;
    current_node->neigh[direction] = son;

    /* process node name (of the son) and branch length (of the edge we just created)... */
    if (EXIT_SUCCESS != process_name_and_brlen(son, in_str, begin, end)) {
        return NULL;
    };

    return son;
} /* end of create_son_and_connect_to_father */



int index_next_toplevel_comma(const char *in_str, int begin, int end) {
    /* returns the index of the next toplevel comma, from position begin included, up to position end.
       the result is -1 if none is found. */
    int level = 0, i;
    for (i = begin; i <= end; i++) {
        switch (in_str[i]) {
            case '(':
                level++;
                break;
            case ')':
                level--;
                break;
            case ',':
                if (level == 0) return i;
        } /* endswitch */
    } /* endfor */
    return -1; /* reached if no outer comma found */
} /* end index_next_toplevel_comma */



int parse_substring_into_node(char *in_str, int begin, int end, Node *current_node, int has_father, Tree *current_tree,
                               size_t nbanno) {
    /* this function supposes that current_node is already allocated, but not the data structures in there.
       It reads starting from character of in_str at index begin and stops at character at index end.
       It is supposed that the input to this function is what has been seen immediately within a set of parentheses.
       The outer parentheses themselves are not included in the range [begin, end].
       So we expect in_str[begin, end] to contain something like:
       MyTaxa:1.2e-3
       OR
       (A:4,B:6)Archae:0.45,Ctax:0.004
       OR
       MyTaxa
       OR
       (A:4,B:6),Ctax
       OR
       A,B,C,D,E,etc (we allow large multifurcations, with no limit on the number of sons)
     */

    /* When called, the current node has just been created but doesn't know yet its number of neighbours. We are going to discover
       this when counting the number of outer commas in the substring. This function:
       (1) checks how many outer commas are here: this is the number of "sons" of this node. Add one to it if the node has a father.
       (2) creates the stuctures (array of node pointers and array of edge pointers) accordingly (+1 for the father)
       (3) fills them. index 0 corresponds to the "father", the other to the "sons". */

    if (begin > end) {
        fprintf(stderr, "Error in parse_substring_into_node: begin > end. Aborting.\n");
        return EXIT_FAILURE;
    }

    int i;
    int pair[2]; /* to be the beginning and end points of the substrings describing the various nodes */
    int inner_pair[2]; /* to be the beginning and end points of the substrings after removing the name and branch length */
    int nb_commas = count_outer_commas(in_str, begin, end);
    int comma_index = begin - 1;
    int direction;
    Node *son;

    /* allocating the data structures for the current node */
    /* FIXME: in case the tree has an inner node with just one child, nb_neigh will become 1 instead of 2.*/
    current_node->nb_neigh = (nb_commas == 0) ? 1 : (nb_commas + 1 + has_father);
    current_node->neigh = malloc(current_node->nb_neigh * sizeof(Node *));

    current_node->bottom_up_likelihood = calloc(nbanno, sizeof(double));
    // do not allocate parsimony_states nor up_parsimony_states as we'll do it when we use them
    current_node->result_probs = calloc(nbanno, sizeof(double));
    current_node->pij = calloc(nbanno, sizeof(double *));
    for (i = 0; i < nbanno; i++) {
        current_node->pij[i] = calloc(nbanno, sizeof(double));
    }
    current_node->best_states = calloc(nbanno, sizeof(size_t));
    current_node->top_down_likelihood = calloc(nbanno, sizeof(double));
    current_node->scaling_factor_down = calloc(1, sizeof(int));
    current_node->scaling_factor_up = calloc(1, sizeof(int));

    if (nb_commas != 0) { /* at least one comma, so at least two sons: */
        for (i = 0; i <= nb_commas; i++) { /* e.g. three iterations for two commas */
            direction = i + has_father;
            pair[0] = comma_index + 1; /* == begin at first iteration */
            comma_index = (i == nb_commas ? end + 1 : index_next_toplevel_comma(in_str, pair[0], end));
            pair[1] = comma_index - 1;

            son = create_son_and_connect_to_father(current_node, current_tree, direction /* dir from current */,
                                                   in_str, pair[0], pair[1]);
            if (son == NULL) {
                return EXIT_FAILURE;
            }
            /* RECURSIVE TREATMENT OF THE SON */
            if (EXIT_SUCCESS != strip_toplevel_parentheses(in_str, pair[0], pair[1],
                                       inner_pair)) {
                return EXIT_FAILURE;
            } /* because name and branch_len already processed by create_son */
            if (EXIT_SUCCESS != parse_substring_into_node(in_str, inner_pair[0], inner_pair[1], son, 1, current_tree,
                                      nbanno)){
                return EXIT_FAILURE;
            } /* recursive treatment */
            /* after the recursive treatment of the son, the data structures of the son have been created, so now we can write
               in it the data corresponding to its direction0 (father) */
            son->neigh[0] = current_node;
        } /* end for i (treatment of the various sons) */
    } /* end if/else on the number of commas */

    return EXIT_SUCCESS;
} /* end parse_substring_into_node */


Tree *parse_nh_string(char *in_str, size_t nbanno) {
    /* this function allocates, populates and returns a new tree. */
    /* returns NULL if the file doesn't correspond to NH format */
    int in_length = (int) strlen(in_str);
    int i; /* loop counter */
    int begin, end; /* to delimitate the string to further process */
    int nodecount = 0;
    size_t n_otu = 0;
    size_t maxpoly;
    double tip_branch_len_sum=0.0;
    Node *cur_node;

    /* SYNTACTIC CHECKS on the input string */
    i = 0;
    while (isspace(in_str[i])) i++;
    if (in_str[i] != '(') {
        fprintf(stderr, "Error: tree doesn't start with an opening parenthesis.\n");
        return NULL;
    }
    else begin = i + 1;
    /* begin: AFTER the very first parenthesis */

    i = in_length - 1;
    while (isspace(in_str[i])) i--;
    if (in_str[i] != ';') {
        fprintf(stderr, "Error: tree doesn't end with a semicolon.\n");
        return NULL;
    }
    while (in_str[--i] != ')');
    end = i - 1;
    /* end: BEFORE the very last parenthesis, discarding optional name for the root and uncanny branch length for its "father" branch */

    /* we make a first pass on the string to discover the number of taxa. */
    /* there are as many OTUs as commas plus 1 in the nh string */
    for (i = 0; i < in_length; i++) if (in_str[i] == ',') n_otu++;
    n_otu++;

    /************************************
    initialisation of the tree structure
    *************************************/
    Tree *t = (Tree *) malloc(sizeof(Tree));
    /* in a rooted binary tree with n taxa, (2n-2) branches and (2n-1) nodes in total.
      this is the maximum we can have. multifurcations will reduce the number of nodes and branches, so set the data structures to the max size */
    t->nb_taxa = n_otu;

    t->nodes = (Node **) calloc(2 * n_otu - 1, sizeof(Node *));
    t->nb_nodes = 1; /* for the moment we only have the root node. */

    t->nb_edges = 0; /* none at the moment */

    t->root = (Node *) malloc(sizeof(Node));
    t->nodes[0] = t->root;

    t->root->id = 0;
    t->root->name = "ROOT";
    t->root->branch_len = 0.0;

    t->next_avail_node_id = 1; /* root node has id 0 */

    /* ACTUALLY READING THE TREE... */
    if (EXIT_SUCCESS != parse_substring_into_node(in_str, begin, end, t->root, 0 /* no father node */, t, nbanno)) {
        return NULL;
    }

    /* SANITY CHECKS AFTER READING THE TREE */

    /* name tree nodes if needed */
    for (i = 0; i < t->nb_nodes; i++) {
        cur_node = t->nodes[i];
        if (cur_node->nb_neigh > 1 && i > 0) {
            nodecount++;
            if (!cur_node->name || strcmp(cur_node->name, "") == 0 || strcmp(cur_node->name, "\0") == 0) {
                sprintf(cur_node->name, "Pastml_Node_%d", nodecount);
            }
        }
    }
    t->min_branch_len = -1.0;
    t->max_branch_len = 0.0;
    t->min_tip_branch_len = -1.0;
    t->num_zero_tip_branches = 0;
    maxpoly=0;
    int num_nonzero_tips = 0;
    int num_nonzero_inner_branches = 0;

    double branch_len_sum = 0.;
    for (i = 0; i < t->nb_nodes; i++) {
        cur_node = t->nodes[i];
        if(isTip(cur_node)) {
          if (cur_node->branch_len > 0.0) {
              tip_branch_len_sum += cur_node->branch_len;
              num_nonzero_tips += 1;
              if ((t->min_tip_branch_len < 0) || (t->min_tip_branch_len > cur_node->branch_len)) {
                  t->min_tip_branch_len = cur_node->branch_len;
              }
          } else {
              t->num_zero_tip_branches += 1;
          }
        }
        if(maxpoly < cur_node->nb_neigh)  {
            maxpoly = cur_node->nb_neigh;
        }
        branch_len_sum += cur_node->branch_len;
        if ((cur_node != t->root) && (cur_node->branch_len > 0.0)) {
            num_nonzero_inner_branches += 1;
            if ((t->min_branch_len < 0) || (t->min_branch_len > cur_node->branch_len)) {
                t->min_branch_len = cur_node->branch_len;
            }
        }
        if (t->max_branch_len < cur_node->branch_len) {
            t->max_branch_len = cur_node->branch_len;
        }
    }
    t->avg_tip_branch_len = tip_branch_len_sum / num_nonzero_tips;
    t->avg_branch_len = branch_len_sum / (num_nonzero_inner_branches + num_nonzero_tips);

    log_info("BASIC TREE STATISTICS:\n\n");
    log_info("\tNumber of taxa:\t%zd\n", t->nb_taxa);
    if (t->nb_taxa > MAXNSP) {
        fprintf(stderr, "Fatal error: too many taxa: more than %d.\n", MAXNSP);
        return NULL;
    }
    log_info("\tNumber of nodes:\t%zd\n", t->nb_nodes - t->nb_taxa);
    log_info("\tNumber of edges:\t%d\n", t->nb_edges);
    log_info("\tMax branch length:\t%e\n", t->max_branch_len);
    log_info("\tAvg branch length:\t%e\n", t->avg_branch_len);
    log_info("\tAvg tip branch length:\t%e\n", t->avg_tip_branch_len);
    log_info("\tNumber of zero tip branches:\t%d\n", t->num_zero_tip_branches);
    log_info("\tMin non-zero inner branch length:\t%e\n", t->min_branch_len);
    log_info("\tMin non-zero tip branch length:\t%e\n", t->min_tip_branch_len);
    log_info("\tMax number of children per node:\t%zd\n", maxpoly-1);
    log_info("\n");

    return t;

} /* end parse_nh_string */


Tree *complete_parse_nh(char *big_string, size_t nbanno) {
    Tree *mytree = parse_nh_string(big_string, nbanno);
    if (mytree == NULL) {
        fprintf(stderr, "Not a syntactically correct NH tree.\n");
        return NULL;
    }

    return mytree;
}


size_t tell_size_of_one_tree(char *filename) {
    /* the only purpose of this is to know about the size of a treefile (NH format)
     * in order to save memspace in allocating the string later on */
    size_t mysize = 0;
    int u;
    FILE *myfile = fopen(filename, "r");
    if (myfile) {
        while ((u = fgetc(myfile)) != ';') { /* termination character of the tree */
            if (feof(myfile)) {
                break;
            } /* shouldn't happen anyway */
            if (!isspace(u)) {
                mysize++;
            }
        }
        fclose(myfile);
    } /* end if(myfile) */
    return mysize + 1;
}

int copy_nh_stream_into_str(FILE *nh_stream, char *big_string) {
    int index_in_string = 0;
    int u;
    /* rewind(nh_stream);
     * DO NOT go to the beginning of the stream
     * if we want to make this flexible enough to read several trees per file */
    while ((u = fgetc(nh_stream)) != ';') { /* termination character of the tree */
        if (feof(nh_stream)) {
            big_string[index_in_string] = '\0';
            return EXIT_FAILURE;
        } /* error code telling that no tree has been read properly */
        if (index_in_string == MAX_TREELENGTH - 1) {
            fprintf(stderr, "Fatal error: tree file seems too big, are you sure it is a newick tree file?\n");
            return EXIT_FAILURE;
        }
        if (!isspace(u)) {
            big_string[index_in_string++] = (char) u;
        }
    }
    big_string[index_in_string++] = ';';
    big_string[index_in_string] = '\0';
    return EXIT_SUCCESS; /* leaves the stream right after the terminal ';' */
} /*end copy_nh_stream_into_str */

void free_node(Node *node, size_t num_anno) {
    size_t j;

    if (node == NULL) return;
    free(node->neigh);
    if (strcmp(node->name, "ROOT") != 0 && node->name) {
        free(node->name);
    }
    free(node->bottom_up_likelihood);
    // do not free parsimony_states as we'll do it as soon as we do not need them
    free(node->result_probs);
    free(node->best_states);
    free(node->top_down_likelihood);
    free(node->scaling_factor_down);
    free(node->scaling_factor_up);
    for (j = 0; j < num_anno; j++) {
        free(node->pij[j]);
    }
    free(node->pij);
    free(node);
}

void free_tree(Tree *tree, size_t num_anno) {
    size_t i;
    if (tree == NULL) return;
    for (i = 0; i < tree->nb_nodes; i++) {
        free_node(tree->nodes[i], num_anno);
    }
    free(tree->nodes);
    free(tree);
}


Tree *read_tree(char *nwk, size_t num_anno) {
    /**
     * Read a tree from newick file
     */

    Tree *s_tree;

    FILE *tree_file = fopen(nwk, "r");
    if (tree_file == NULL) {
        fprintf(stderr, "Tree file %s is not found or is impossible to access.\n", nwk);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return NULL;
    }

    size_t tree_file_size = 3 * tell_size_of_one_tree(nwk);
    if (tree_file_size > MAX_TREELENGTH) {
        fprintf(stderr, "Tree filesize for %s is more than %d bytes: are you sure it's a valid newick tree?\n",
                nwk, MAX_TREELENGTH / 3);
        return NULL;
    }

    void *retval;
    if ((retval = calloc(tree_file_size + 1, sizeof(char))) == NULL) {
        fprintf(stderr, "Not enough memory\n");
        return NULL;
    }
    char *c_tree = (char *) retval;

    if (EXIT_SUCCESS != copy_nh_stream_into_str(tree_file, c_tree)) {
        fprintf(stderr, "A problem occurred while parsing the reference tree.\n");
        return NULL;
    }
    fclose(tree_file);

    /*Make Tree structure*/
    s_tree = complete_parse_nh(c_tree, num_anno);
    if (NULL == s_tree) {
        fprintf(stderr, "A problem occurred while parsing the reference tree.\n");
        return NULL;
    }
    return s_tree;
}


int dir_a_to_b(Node *a, Node *b) {
    /* this returns the direction from a to b when a and b are two neighbours, otherwise return -1 */
    int i, n = a->nb_neigh;
    for (i = 0; i < n; i++) if (a->neigh[i] == b) break;
    if (i < n) return i;
    else {
        return -1;
    }
} /* end dir_a_to_b */

int write_subtree_to_stream(Node *node, Node *node_from, FILE *stream) {
    int i, direction_to_exclude, n = node->nb_neigh;
    if (node_from == NULL) {
        return EXIT_SUCCESS;
    }
    // internal node => write its children first
    if (n != 1) {
        direction_to_exclude = dir_a_to_b(node, node_from);
        if (-1 == direction_to_exclude) {
            fprintf(stderr, "Fatal error : nodes are not neighbours.\n");
            return EXIT_FAILURE;
        }

        putc('(', stream);
        /* we have to write (n-1) subtrees in total. The last print is not followed by a comma */
        for (i = 1; i < n - 1; i++) {
            if (EXIT_SUCCESS !=
                write_subtree_to_stream(node->neigh[(direction_to_exclude + i) % n], node, stream)) {
                return EXIT_FAILURE;
            } /* a son */
            putc(',', stream);
        }
        if (EXIT_SUCCESS != write_subtree_to_stream(node->neigh[(direction_to_exclude + i) % n], node, stream)) {
            return EXIT_FAILURE;
        } /* last son */
        putc(')', stream);
    }
    // write node's name and dist to father
    fprintf(stream, "%s:%f", (node->name ? node->name : ""), node->branch_len);
    return EXIT_SUCCESS;
} /* end write_subtree_to_stream */

int write_nh_tree(Tree *s_tree, char *output_filepath) {

    FILE* output_file = fopen(output_filepath, "w");
    if (!output_file) {
        fprintf(stderr, "Output tree file %s is impossible to access.", output_filepath);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    /* writing the tree from the current position in the file */
    int i, n = s_tree->root->nb_neigh;
    putc('(', output_file);
    for (i = 0; i < n - 1; i++) {
        if (EXIT_SUCCESS != write_subtree_to_stream(s_tree->root->neigh[i], s_tree->root, output_file)) {
            return EXIT_FAILURE;
        } /* a son */
        putc(',', output_file);
    }
    if (EXIT_SUCCESS != write_subtree_to_stream(s_tree->root->neigh[i], s_tree->root, output_file)) {
        return EXIT_FAILURE;
    } /* last son */
    putc(')', output_file);

    if (s_tree->root->name) {
        fprintf(output_file, "%s", s_tree->root->name);
    }
    /* terminate with a semicol AND and end of line */
    putc(';', output_file);
    putc('\n', output_file);

    fclose(output_file);
    return EXIT_SUCCESS;
}
