#include "pastml.h"
#include <assert.h>

double BL_sum = 0., BL_avg = 0.;

void free_node(Node *node, int count, int num_anno);


int index_toplevel_colon(char *in_str, int begin, int end) {
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

int strip_toplevel_parentheses(char *in_str, int begin, int end, int *pair) {
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

int count_outer_commas(char *in_str, int begin, int end) {
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
        son_node->brlen = 0.0;
    } else {
        if (EXIT_SUCCESS != parse_double(in_str, colon + 1, end, &brlen)) {
            return EXIT_FAILURE;
        }
        son_node->brlen = brlen;
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
    if (opening_bracket != -1) name_end = opening_bracket - 1; else name_end = (colon == -1 ? end : colon - 1);
    /* but now if the name starts and ends with single or double quotes, remove them */
    if (in_str[name_begin] == in_str[name_end] && (in_str[name_begin] == '"' || in_str[name_begin] == '\'')) {
        name_begin++;
        name_end--;
    }
    name_length = name_end - name_begin + 1;
    effective_length = (name_length > MAX_NAMELENGTH ? MAX_NAMELENGTH : name_length);
    son_node->name = (char *) malloc((MAX_NAMELENGTH) * sizeof(char));
    if (name_length >= 1) {
        strncpy(son_node->name, in_str + name_begin, effective_length);
        son_node->name[effective_length] = '\0'; /* terminating the string */
    }
    return EXIT_SUCCESS;
} /* end of process_name_and_brlen */


Node *create_son_and_connect_to_father(Node *current_node, Tree *current_tree, int direction, char *in_str, int begin,
                                       int end) {
    /* This function creates (allocates) the son node in the given direction from the current node.
       It also creates a new branch to connect the son to the father.
       The array structures in the tree (a_nodes) are updated accordingly.
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
    current_tree->a_nodes[son->id] = son;
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



int index_next_toplevel_comma(char *in_str, int begin, int end) {
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
                               int nbanno) {
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

    int i, j;
    int pair[2]; /* to be the beginning and end points of the substrings describing the various nodes */
    int inner_pair[2]; /* to be the beginning and end points of the substrings after removing the name and branch length */
    int nb_commas = count_outer_commas(in_str, begin, end);
    int comma_index = begin - 1;
    int direction;
    Node *son;

    /* allocating the data structures for the current node */
    /* FIXME: in case the tree has an inner node with just one child, nneigh will become 1 instead of 2.*/
    current_node->nneigh = (nb_commas == 0 ? 1 : nb_commas + 1 + has_father);
    current_node->neigh = malloc(current_node->nneigh * sizeof(Node *));

    size_t nbanno_size_t = (size_t) nbanno;
    current_node->bottom_up_likelihood = calloc(nbanno_size_t, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->bottom_up_likelihood[i] = 0.0;
    current_node->condlike_mar = calloc(nbanno_size_t, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->condlike_mar[i] = 0.0;
    current_node->pij = calloc(nbanno_size_t, sizeof(double *));
    for (i = 0; i < nbanno; i++) current_node->pij[i] = calloc(nbanno_size_t, sizeof(double));
    for (i = 0; i < nbanno; i++) {
      for (j = 0; j < nbanno; j++) current_node->pij[i][j] = 0.0;
    }

    current_node->marginal = calloc(nbanno_size_t, sizeof(double));
    for (i = 0; i < nbanno; i++) {
        current_node->marginal[i] = 0.0;
    }
    current_node->best_states = calloc(nbanno_size_t, sizeof(int));
    for (i = 0; i < nbanno; i++) {
        current_node->best_states[i] = 0;
    }
    current_node->top_down_likelihood = calloc(nbanno_size_t, sizeof(double));
    for (i = 0; i < nbanno; i++) {
        current_node->top_down_likelihood[i] = 0.0;
    }

    if (has_father == 0) { /* root */
    }
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
            } /* because name and brlen already processed by create_son */
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


Tree *parse_nh_string(char *in_str, int nbanno) {
    /* this function allocates, populates and returns a new tree. */
    /* returns NULL if the file doesn't correspond to NH format */
    int in_length = (int) strlen(in_str);
    int i; /* loop counter */
    int begin, end; /* to delimitate the string to further process */
    int n_otu = 0, nodecount = 0;
    int maxpoly;
    double tip_branch_len_sum=0.0;

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

    t->a_nodes = (Node **) calloc(2 * n_otu - 1, sizeof(Node *));
    t->nb_nodes = 1; /* for the moment we only have the node0 node. */

    t->nb_edges = 0; /* none at the moment */

    t->node0 = (Node *) malloc(sizeof(Node));
    t->a_nodes[0] = t->node0;

    t->node0->id = 0;
    t->node0->name = "ROOT";

    t->next_avail_node_id = 1; /* root node has id 0 */

    /* ACTUALLY READING THE TREE... */
    if (EXIT_SUCCESS != parse_substring_into_node(in_str, begin, end, t->node0, 0 /* no father node */, t, nbanno)) {
        return NULL;
    }

    /* SANITY CHECKS AFTER READING THE TREE */

    /* name tree nodes if needed */
    for (i = 0; i < t->nb_nodes; i++) {
        if (t->a_nodes[i]->nneigh > 1 && i > 0) {
            nodecount++;
            if (!t->a_nodes[i]->name || strcmp(t->a_nodes[i]->name, "") == 0 || strcmp(t->a_nodes[i]->name, "\0") == 0) {
                sprintf(t->a_nodes[i]->name, "Pastml_Node_%d", nodecount);
            }
        }
    }
    t->min_bl = DBL_MAX;
    maxpoly=0;
    for (i = 0; i < t->nb_nodes; i++) {
        if(t->a_nodes[i]->nneigh == 1){ //tips
          tip_branch_len_sum += t->a_nodes[i]->brlen;
        }
        if (t->a_nodes[i]->nneigh - 1 > MAXPOLY) {
            fprintf(stderr, "Fatal error: too many polytomy more than %d at the node %s.\n", MAXPOLY,
                    t->a_nodes[i]->name);
            return NULL;
        }
        if(maxpoly < t->a_nodes[i]->nneigh) maxpoly = t->a_nodes[i]->nneigh;
        if (t->min_bl > t->a_nodes[i]->brlen && t->a_nodes[i]->brlen != 0.0)
            t->min_bl = t->a_nodes[i]->brlen;
        if (i != 0) {
            BL_sum += t->a_nodes[i]->brlen;
        }
    }
    t->tip_avg_branch_len = tip_branch_len_sum / (double) t->nb_taxa;
    BL_avg = BL_sum / (double) t->nb_edges;
    t->avg_branch_len = BL_avg;

    printf("\n*** BASIC STATISTICS ***\n\n", in_str);
    printf("Number of taxa in the tree read: %d\n", t->nb_taxa);
    if (t->nb_taxa > MAXNSP) {
        fprintf(stderr, "Fatal error: too many taxa: more than %d.\n", MAXNSP);
        return NULL;
    }
    printf("Number of nodes in the tree read: %d\n", t->nb_nodes - t->nb_taxa);
    printf("Number of edges in the tree read: %d\n", t->nb_edges);
    printf("Average branch lengths of the tree: %lf\n", t->avg_branch_len);
    printf("Minimum branch length in the tree: %lf\n", t->min_bl);
    printf("The maximum number of multifurcations at a single node of the input tree: %d\n", maxpoly-1);

    return t;

} /* end parse_nh_string */


Tree *complete_parse_nh(char *big_string, int nbanno) {
    Tree *mytree = parse_nh_string(big_string, nbanno);
    if (mytree == NULL) {
        fprintf(stderr, "Not a syntactically correct NH tree.\n");
        return NULL;
    }

    return mytree;
}
