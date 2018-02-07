#include "pastml.h"
#include <assert.h>

int ntax;
double BL_sum = 0., BL_avg = 0.;

void free_node(Node *node, int count, int num_anno);

void free_edge(Edge *edge);

int min_int(int a, int b) {
    return (a < b ? a : b);
}


/* collapsing a branch */
void collapse_branch(Edge *branch, Tree *tree, int count, int nbanno) {
    /* this function collapses the said branch and creates a higher-order multifurcation (n1 + n2 - 2 neighbours for the resulting node).
       We also have to remove the extra node from tree->a_nodes and the extra edge from t->a_edges.
       to be done:
       (1) create a new node with n1+n2-2 neighbours. Ultimately we will destroy the original node.
       (2) populate its list of neighbours from the lists of neighbours corresponding to the two original nodes
       (3) populate its list of neighbouring edges form the br lists of the two original nodes
       (4) for each of the neighbours, set the info regarding their new neighbour (that is, our new node)
       (5) for each of the neighbouring branches, set the info regarding their new side (that is, our new node)
       (6) destroy the two original nodes and commit this info to a_nodes. Modify tree->nb_nodes
       (7) destroy the original edge and commit this info to a_edges. Modify tree->nb_edges */

    /* WARNING: this function won't accept to collapse terminal edges */
    Node *node1 = branch->left, *node2 = branch->right;
    int i, j, n1 = node1->nneigh, n2 = node2->nneigh;
    if (n1 == 1 || n2 == 1) {
        fprintf(stderr, "Warning: %s() won't collapse terminal edges.\n", __FUNCTION__);
        return;
    }
    int degree = n1 + n2 - 2;
    /* (1) */
    /* Node* new = new_node("collapsed", tree, n1 + n2 - 2); */ /* we cannot use that because we want to reuse n1's spot in tree->a_nodes */
    //printf("start collapse");
    Node *new = (Node *) malloc(sizeof(Node));
    new->nneigh = degree;
    new->neigh = malloc(degree * sizeof(Node *));
    new->br = malloc(degree * sizeof(Edge *));
    new->id = node1->id; /* because we are going to store the node at this index in tree->a_nodes */
    new->name = strdup("collapsed");
    new->comment = NULL;
    new->depth = min_int(node1->depth, node2->depth);

    /* very important: set tree->node0 to new in case it was either node1 or node2 */
    if (tree->node0 == node1 || tree->node0 == node2) tree->node0 = new;


    int ind = 0; /* index in the data structures in new */
    /* (2) and (3) and (4) and (5) */
    for (i = 0; i < n1; i++) {
        if (node1->neigh[i] == node2) continue;
        new->neigh[ind] = node1->neigh[i];
        /*  then change one of the neighbours of that neighbour to be the new node... */
        for (j = 0; j < new->neigh[ind]->nneigh; j++) {
            if (new->neigh[ind]->neigh[j] == node1) {
                new->neigh[ind]->neigh[j] = new;
                break;
            }
        } /* end for j */

        new->br[ind] = node1->br[i];
        /* then change one of the two ends of that branch to be the new node... */
        if (new->neigh[ind] == new->br[ind]->right) new->br[ind]->left = new; else new->br[ind]->right = new;
        ind++;
    }

    for (i = 0; i < n2; i++) {
        if (node2->neigh[i] == node1) continue;
        new->neigh[ind] = node2->neigh[i];
        /*  then change one of the neighbours of that neighbour to be the new node... */
        for (j = 0; j < new->neigh[ind]->nneigh; j++) {
            if (new->neigh[ind]->neigh[j] == node2) {
                new->neigh[ind]->neigh[j] = new;
                break;
            }
        } /* end for j */

        new->br[ind] = node2->br[i];
        /* then change one of the two ends of that branch to be the new node... */
        if (new->neigh[ind] == new->br[ind]->right) new->br[ind]->left = new; else new->br[ind]->right = new;
        ind++;
    }

    /* (6) tidy up tree->a_nodes and destroy old nodes */
    assert(tree->a_nodes[new->id] == node1);
    tree->a_nodes[new->id] = new;
    /* current last node in tree->a_edges changes id and is now placed at the position were node2 was */
    int id2 = node2->id;
    assert(tree->a_nodes[id2] == node2);
    tree->a_nodes[id2] = tree->a_nodes[--tree->next_avail_node_id]; /* moving the last node into the spot occupied by node2... */
    tree->a_nodes[id2]->id = id2;                    /* and changing its id accordingly */
    tree->a_nodes[tree->next_avail_node_id] = NULL; /* not strictly necessary, but... */
    tree->nb_nodes--;
    //printf("finished before free");
    free_node(node1, 1, nbanno);
    free_node(node2, 1, nbanno);

    /* (7) tidy up tree->a_edges and destroy the old branch */
    assert(tree->a_edges[branch->id] == branch);
    tree->a_edges[branch->id] = tree->a_edges[--tree->next_avail_edge_id]; /* moving the last branch into the spot occupied by 'branch' */
    tree->a_edges[branch->id]->id = branch->id;                /* ... and changing its id accordingly */
    tree->a_edges[tree->next_avail_edge_id] = NULL; /* not strictly necessary, but... */
    tree->nb_edges--;
    free_edge(branch);

} /* end collapse_branch */

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
        return;
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

int count_zero_length_branches(Tree *tree) {
    int count = 0;
    int i, n = tree->nb_edges;
    for (i = 0; i < n; i++) if (tree->a_edges[i]->had_zero_length) count++;
    return count;
}

int count_leaves(Tree *tree) {
    int count = 0;
    int i, n = tree->nb_nodes;
    for (i = 0; i < n; i++) if (tree->a_nodes[i]->nneigh == 1) count++;
    return count;
}

int count_roots(Tree *tree) { /* to ensure there is exactly zero or one root */
    int count = 0;
    int i, n = tree->nb_nodes;
    for (i = 0; i < n; i++) {
        if (tree->a_nodes[i]->nneigh > 1) {
            if (strcmp(tree->a_nodes[i]->neigh[1]->name, "Node2") == 0) count++;
        }
    }
    return count;
}

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

int process_name_and_brlen(Node *son_node, Edge *edge, Tree *current_tree, char *in_str, int begin, int end) {
    /* looks into in_str[begin..end] for the branch length of the "father" edge
       and updates the edge and node structures accordingly */
    int colon = index_toplevel_colon(in_str, begin, end);
    int closing_par = -1, opening_bracket = -1;
    int i, ignore_mode, name_begin, name_end, name_length, effective_length;
    double brlen = .0;

    /* processing the optional BRANCH LENGTH... */
    if (colon == -1) {
        edge->had_zero_length = TRUE;
        edge->brlen = 0.0;
    } else {
        if (EXIT_SUCCESS != parse_double(in_str, colon + 1, end, &brlen)) {
            return EXIT_FAILURE;
        }
        edge->had_zero_length = (brlen == 0.0);
        edge->brlen = brlen;
        //edge->brlen = (brlen < MIN_BRLEN ? MIN_BRLEN : brlen);
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
       The array structures in the tree (a_nodes and a_edges) are updated accordingly.
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

    int i;
    Node *son = (Node *) malloc(sizeof(Node));
    son->id = current_tree->next_avail_node_id++;
    current_tree->a_nodes[son->id] = son;
    current_tree->nb_nodes++;

    son->name = son->comment = NULL;
    son->depth = MAX_NODE_DEPTH;

    Edge *edge = (Edge *) malloc(sizeof(Edge));
    edge->id = current_tree->next_avail_edge_id++;
    current_tree->a_edges[edge->id] = edge;
    current_tree->nb_edges++;

    // for (i=0; i<2; i++) edge->subtype_counts[i] = (int*) calloc(NUM_SUBTYPES, sizeof(int));
    for (i = 0; i < 2; i++) edge->subtype_counts[i] = NULL; /* subtypes.c will have to create that space */

    edge->right = son;
    edge->left = current_node;

    current_node->neigh[direction] = son;
    current_node->br[direction] = edge;

    /* process node name (of the son) and branch length (of the edge we just created)... */
    if (EXIT_SUCCESS != process_name_and_brlen(son, edge, current_tree, in_str, begin, end)) {
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
    current_node->nneigh = (nb_commas == 0 ? 1 : nb_commas + 1 + has_father);
    current_node->neigh = malloc(current_node->nneigh * sizeof(Node *));
    current_node->br = malloc(current_node->nneigh * sizeof(Edge *));

    current_node->condlike = calloc(nbanno, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->condlike[i] = 0.0;
    current_node->condlike_mar = calloc(nbanno, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->condlike_mar[i] = 0.0;
    current_node->up_like = calloc(nbanno, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->up_like[i] = 0.0;
    current_node->mar_prob = calloc(nbanno, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->mar_prob[i] = 0.0;
    current_node->pij = calloc(nbanno, sizeof(double *));
    for (i = 0; i < nbanno; i++) current_node->pij[i] = calloc(nbanno, sizeof(double));
    for (i = 0; i < nbanno; i++) {
      for (j = 0; j < nbanno; j++) current_node->pij[i][j] = 0.0;
    }

    current_node->marginal = calloc(nbanno, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->marginal[i] = 0.0;
    current_node->tmp_best = calloc(nbanno, sizeof(int)); for (i = 0; i < nbanno; i++) current_node->tmp_best[i] = 0;
    current_node->sum_down = calloc(nbanno, sizeof(double)); for (i = 0; i < nbanno; i++) current_node->sum_down[i] = 0.0;

    current_node->rootpij = calloc(nbanno, sizeof(double *));
    for (i = 0; i < nbanno; i++) current_node->rootpij[i] = calloc(nbanno, sizeof(double));
    for (i = 0; i < nbanno; i++) {
      for (j = 0; j < nbanno; j++) current_node->rootpij[i][j] = 0.0;
    }

    current_node->local_flag = calloc(nbanno, sizeof(int));
    for (i = 0; i < nbanno; i++) current_node->local_flag[i] = 1;


    if (has_father == 0) { /* root */
    }
    if (nb_commas == 0) { /* leaf: no recursive call */
        /* this means there is no split here, terminal node: we know that the current node is a leaf.
           Its name is already there in node->name, we just have to update the taxname table and all info related
           to the fact that we have a taxon here. */
        /* that's also the moment when we check that there are no two identical taxa on different leaves of the tree */
        for (i = 0; i < current_tree->next_avail_taxon_id; i++) {
            if (!strcmp(current_node->name, current_tree->taxa_names[i])) {
                fprintf(stderr, "Fatal error: duplicate taxon %s.\n", current_node->name);
                return EXIT_FAILURE;
            } /* end if */
        } /* end for */
        current_tree->taxa_names[current_tree->next_avail_taxon_id++] = strdup(current_node->name);
    } else { /* at least one comma, so at least two sons: */
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
            son->br[0] = current_node->br[direction];
        } /* end for i (treatment of the various sons) */
    } /* end if/else on the number of commas */

    return EXIT_SUCCESS;
} /* end parse_substring_into_node */




Tree *parse_nh_string(char *in_str, int nbanno, char *keep_ID) {
    /* this function allocates, populates and returns a new tree. */
    /* returns NULL if the file doesn't correspond to NH format */
    int in_length = (int) strlen(in_str);
    int i, j; /* loop counter */
    int begin, end; /* to delimitate the string to further process */
    int n_otu = 0, nodecount = 0;
    char str[MAXLNAME];
    int collapsed_one = 0, uncollapsed_terminal = 0, collapsed_internal = 0, maxpoly;

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

    /* immediately, we set the global variable ntax. TODO: see if we can simply get rid of this global var. */
    ntax = n_otu;



    /************************************
    initialisation of the tree structure
    *************************************/
    Tree *t = (Tree *) malloc(sizeof(Tree));
    /* in a rooted binary tree with n taxa, (2n-2) branches and (2n-1) nodes in total.
      this is the maximum we can have. multifurcations will reduce the number of nodes and branches, so set the data structures to the max size */
    t->nb_taxa = n_otu;

    t->a_nodes = (Node **) calloc(2 * n_otu - 1, sizeof(Node *));
    t->nb_nodes = 1; /* for the moment we only have the node0 node. */

    t->a_edges = (Edge **) calloc(2 * n_otu - 2, sizeof(Edge *));
    t->nb_edges = 0; /* none at the moment */

    t->node0 = (Node *) malloc(sizeof(Node));
    t->a_nodes[0] = t->node0;

    t->node0->id = 0;
    t->node0->name = "ROOT";
    t->node0->comment = NULL;

    t->node0->depth = MAX_NODE_DEPTH;
    t->taxa_names = (char **) malloc(n_otu * sizeof(char *));

    t->next_avail_node_id = 1; /* root node has id 0 */
    t->next_avail_edge_id = 0; /* no branch added so far */
    t->next_avail_taxon_id = 0; /* no taxon added so far */

    /* ACTUALLY READING THE TREE... */
    if (EXIT_SUCCESS != parse_substring_into_node(in_str, begin, end, t->node0, 0 /* no father node */, t, nbanno)) {
        return NULL;
    }

    /* SANITY CHECKS AFTER READING THE TREE */

    if (strcmp(keep_ID, "F") == 0) {
        for (i = 0; i < t->nb_nodes; i++) {
            if (t->a_nodes[i]->nneigh > 1 && i > 0) {
                nodecount++;
                sprintf(t->a_nodes[i]->name, "%s%d", "Node", nodecount);
            }
        }
    }
    t->min_bl = DBL_MAX;
    maxpoly=0;
    for (i = 0; i < t->nb_nodes; i++) {
        if (t->a_nodes[i]->nneigh > MAXPOLY) {
            fprintf(stderr, "Fatal error: too many polytomy more than %d at the node %s.\n", MAXPOLY,
                    t->a_nodes[i]->name);
            return NULL;
        }
        if(maxpoly < t->a_nodes[i]->nneigh) maxpoly = t->a_nodes[i]->nneigh;
        if (t->min_bl > t->a_nodes[i]->br[0]->brlen && t->a_nodes[i]->br[0]->brlen != 0.0)
            t->min_bl = t->a_nodes[i]->br[0]->brlen;
        for (j = 0; j < t->a_nodes[i]->nneigh; j++) {
            if (j == 0) {
                if (i == 0) {
                } else {
                    BL_sum += t->a_nodes[i]->br[j]->brlen;
                }
            } else {
            }
        }
    }
    BL_avg = (double) BL_sum / (double) t->nb_edges;
    t->avgbl = BL_avg;
    //printf("Tree Length = %lf\n", BL_sum);

    printf("\n*** BASIC STATISTICS ***\n\n", in_str);
    printf("Number of taxa in the tree read: %d\n", t->nb_taxa);
    if (t->nb_taxa > MAXNSP) {
        fprintf(stderr, "Fatal error: too many taxa: more than %d.\n", MAXNSP);
        return NULL;
    }
    //printf("Number of leaves according to the tree structure: %d\n", count_leaves(t));
    printf("Number of nodes in the tree read: %d\n", t->nb_nodes - t->nb_taxa);
    printf("Number of edges in the tree read: %d\n", t->nb_edges);
    printf("Average branch lengths of the tree: %lf\n", t->avgbl);
    printf("Minimum branch length in the tree: %lf\n", t->min_bl);
    printf("Number of edges with zero length: %d\n", count_zero_length_branches(t));
    printf("The maximum number of multifurcations at a single node of the input tree: %d\n", maxpoly-1);

    return t;

} /* end parse_nh_string */




Tree *complete_parse_nh(char *big_string, int nbanno, char *keep_ID) {
    int i;
    Tree *mytree = parse_nh_string(big_string, nbanno, keep_ID);
    if (mytree == NULL) {
        fprintf(stderr, "Not a syntactically correct NH tree.\n");
        return NULL;
    }

    return mytree;
}
