#include "pastml.h"
#include "tree.h"


bool intersect_with_states2(long *states1, const long *states2, size_t n) {
    bool intersection_is_empty = true;
    size_t k;

    for (k = 0; k < n; k++) {
        states1[k] &= states2[k];
        if (states1[k] != 0l) {
            intersection_is_empty = false;
        }
    }
    return !intersection_is_empty;
}

long* get_neighbour_state_intersection(Node *node, size_t first_index, size_t last_index,
                                       size_t mask_array_len, int include_self,
                                       long *(*get_neighbour_states)(Node *)) {
    /**
     * Calculates child state intersection and if it's non-empty sets the node->parsimony_states to it.
     * Returns true is the intersection is non-empty else false.
     */
    size_t i, k;
    bool intersection_is_empty = true;
    long* neighbour_states;

    long *state_intersection = calloc(mask_array_len, sizeof(long));
    if (include_self) {
        // set intersection mask to current node mask if it should be included
        memcpy(state_intersection, get_neighbour_states(node), mask_array_len * sizeof(long));
    } else {
        // set all intersection bits to 1
        for (k = 0; k < mask_array_len; k++) {
            state_intersection[k] = ~0l;
        }
    }

    for (i = first_index; i < last_index; i++) {
        neighbour_states = get_neighbour_states(node->neigh[i]);
        if (neighbour_states != NULL) {
            intersection_is_empty = !intersect_with_states2(state_intersection, neighbour_states, mask_array_len);
            // if it's already empty no need to go further
            if (intersection_is_empty) {
                break;
            }
        }
    }
    if (intersection_is_empty) {
        free(state_intersection);
    }
    return intersection_is_empty ? NULL: state_intersection;
}
size_t update_counts(Node* nd, size_t mc);

long * get_most_common_neighbour_states(Node *node, size_t first_index, size_t last_index, size_t num_annotations,
                                        int include_self, long *(*get_neighbour_states)(Node *)) {
    size_t i, k;
    long cur_state = 1;
    size_t max_count;
    size_t *state2count;
    long* most_common_states;
    size_t mask_array_len = (size_t) ceil(num_annotations / 64.0);

    most_common_states = get_neighbour_state_intersection(node, first_index, last_index, mask_array_len, include_self,
                                                          get_neighbour_states);
    if (most_common_states != NULL) {
        return most_common_states;
    }

    state2count = calloc(num_annotations, sizeof(size_t));
    most_common_states = calloc(mask_array_len, sizeof(long));
    max_count = 0;

    size_t update_counts(Node* nd, size_t mc) {
        long cur_mask = 1;
        long* neightbour_states = get_neighbour_states(nd);
        if (neightbour_states != NULL) {
            for (i = 0; i < num_annotations; i++) {
                if (i % 64 == 0) {
                    cur_mask = 1;
                } else {
                    cur_mask <<= 1;
                }
                if ((neightbour_states[i / 64] & cur_mask) != 0l) {
                    state2count[i] += 1;
                    if (state2count[i] > mc) {
                        mc = state2count[i];
                    }
                }
            }
        }
        return mc;
    }

    if (include_self) {
        // set counts to current state counts
        max_count = update_counts(node, max_count);
    }
    for (k = first_index; k < last_index; k++) {
        max_count = update_counts(node->neigh[k], max_count);
    }

    for (i = 0; i < num_annotations; i++) {
        if (i % 64 == 0) {
            cur_state = 1;
        } else {
            cur_state <<= 1;
        }
        // if this state is (one of) the most frequent, set the corresponding bit to 1
        if (state2count[i] == max_count) {
            most_common_states[i / 64] |= cur_state;
        }
    }
    free(state2count);
    return most_common_states;
}


long* get_down_parsimony_states(Node *nd) {
    return nd->down_parsimony_states;
}

long* get_parsimony_states(Node *nd) {
    return nd->parsimony_states;
}

void _uppass_node(Node* node, Node* root, size_t num_annotations) {
    size_t i;
    long cur_state;
    size_t mask_array_len = (size_t) ceil(num_annotations / 64.0);

    if (isTip(node)) {
        // initialize the bit mask for the tips: put 1 in the position of its state(s).
        node->down_parsimony_states = calloc(mask_array_len, sizeof(long));
        cur_state = 1;
        for (i = 0; i < num_annotations; i++) {
            if (i % 64 == 0) {
                cur_state = 1;
            } else {
                cur_state <<= 1;
            }
            if (node->bottom_up_likelihood[i] == 1.0) {
                node->down_parsimony_states[i / 64] |= cur_state;
            }
        }
        return;
    }
    // process children
    for (i = (node == root) ? 0: 1; i < node->nb_neigh; i++) {
        _uppass_node(node->neigh[i], root, num_annotations);
    }
    node->down_parsimony_states = get_most_common_neighbour_states(node, (node == root) ? 0 : 1, node->nb_neigh,
                                                                       num_annotations, false,
                                                                       get_down_parsimony_states);
}

long *get_states(Node *nd);
long *get_up_down_states(Node *nd);

void _downpass_node(Node* node, Node* root, size_t num_annotations) {
    /**
     * if N is not a tip:
     *      L, R <- left and right children of N
     *      if N is root:
     *          UP_S(N) <- union of all states
     *      else:
     *          P <- parent of N
     *          B <- brother of N
     *          UP_S(N) <- most_common_states(UP_S(P), S(B))
     *          S(N) <- most_common_states(UP_S(N), S(L), S(R))
     *     DOWNPASS(L)
     *     DOWNPASS(R)
     */
    size_t i;
    size_t mask_array_len = (size_t) ceil(num_annotations / 64.0);
    Node* parent;

    if (node == root) {
        node->up_parsimony_states = calloc(mask_array_len, sizeof(long));
        // set up_parsimony_states to all 1
        for (i = 0; i < mask_array_len; i++) {
            node->up_parsimony_states[i] = ~0l;
        }
        node->parsimony_states = node->down_parsimony_states;
    } else {
        parent = node->neigh[0];

        long *get_states(Node *nd) {
            // parent
            if (nd == parent) {
                return nd->up_parsimony_states;
            }
            // do not include this node
            if (nd == node) {
                return NULL;
            }
            // brother
            return nd->down_parsimony_states;
        }

        node->up_parsimony_states = get_most_common_neighbour_states(parent, (parent == root) ? 0 : 1, parent->nb_neigh,
                                                                     num_annotations, true, get_states);
        if (isTip(node)) {
            /* if we were hesitating between several states for this tip before,
             * we might be able to reduce the number of possibilities
             * by intersecting them with the up state */
            if (intersect_with_states2(node->up_parsimony_states, node->down_parsimony_states, mask_array_len)) {
                node->parsimony_states = node->up_parsimony_states;
            } else {
                node->parsimony_states = node->down_parsimony_states;
            }
        } else {
            long *get_up_down_states(Node *nd) {
                // this node
                if (nd == node) {
                    return nd->up_parsimony_states;
                }
                // child
                return nd->down_parsimony_states;
            }

            /* find most common states among node->up_parsimony_states and its children parsimony_states,
             * and set node->parsimony_states to it. */
            node->parsimony_states = get_most_common_neighbour_states(node, (node == root) ? 0 : 1, node->nb_neigh,
                                                                      num_annotations, true, get_up_down_states);
        }
    }
    // process children
    for (i = (node == root) ? 0: 1; i < node->nb_neigh; i++) {
        _downpass_node(node->neigh[i], root, num_annotations);
    }
}


void _acctran_node(Node* node, Node* root, size_t num_annotations) {
    size_t i;
    Node *child;
    size_t mask_array_len = (size_t) ceil(num_annotations / 64.0);
    long* state_intersection;

    for (i = (node == root) ? 0: 1; i < node->nb_neigh; i++) {
        child = node->neigh[i];
        // set child state to its intersection with its parent (this node) state, if non-empty
        state_intersection = get_neighbour_state_intersection(child, 0, 1, mask_array_len, true,
                                                              get_down_parsimony_states);
        if (state_intersection != NULL) {
            free(child->down_parsimony_states);
            child->down_parsimony_states = state_intersection;
        }
        // acctran child
        _acctran_node(child, root, num_annotations);
    }
    node->parsimony_states = node->down_parsimony_states;
}


void _deltran_node(Node* node, Node* root, size_t num_annotations) {
    size_t i;
    size_t mask_array_len = (size_t) ceil(num_annotations / 64.0);
    long* state_intersection;

    if (node != root) {
        // set node state to its intersection with its parent, if non-empty
        state_intersection = get_neighbour_state_intersection(node, 0, 1, mask_array_len, true, get_parsimony_states);
        if (state_intersection != NULL) {
            free(node->parsimony_states);
            node->parsimony_states = state_intersection;
        }
    }
    // deltran children
    for (i = (node == root) ? 0: 1; i < node->nb_neigh; i++) {
        _deltran_node(node->neigh[i], root, num_annotations);
    }
}


void uppass(Tree* tree, size_t num_annotations) {
    /**
     * UPPASS traverses the tree starting from the tips and going up till the root,
     * and assigns to each parent node a state based on the states of its child nodes.
     *
     *   if N is a tip:
     *   S(N) <- state of N
     *   else:
     *       L, R <- left and right children of N
     *       UPPASS(L)
     *       UPPASS(R)
     *       if S(L) intersects with S(R):
     *          S(N) <- intersection(S(L), S(R))
     *       else:
     *          S(N) <- union(S(L), S(R))
     */
    _uppass_node(tree->root, tree->root, num_annotations);
}

void downpass(Tree* tree, size_t num_annotations) {
    /**
     * DOWNPASS traverses the tree starting from the root and going down till the tips,
     * and for each node combines the state information from its supertree and its subtree (calculated at UPPASS).
     * As the root state was already the most parsimonious after the UPPASS,
     * we skip it and start directly with the root children.
     *
     * if N is not a tip:
     *      L, R <- left and right children of N
     *      if N is root:
     *          UP_S(N) <- union of all states
     *      else:
     *          P <- parent of N
     *          B <- brother of N
     *          UP_S(N) <- most_common_states(UP_S(P), S(B))
     *          S(N) <- most_common_states(UP_S(N), S(L), S(R))
     *     DOWNPASS(L)
     *     DOWNPASS(R)
     */
    size_t i;
    Node * node;
    _downpass_node(tree->root, tree->root, num_annotations);
    for (i = 0; i < tree->nb_nodes; i++) {
        node = tree->nodes[i];
        // do not need our up_parsimony_states nor down_parsimony_states anymore
        if (node->parsimony_states != node->up_parsimony_states) {
            free(node->up_parsimony_states);
        }
        if (node->parsimony_states != node->down_parsimony_states) {
            free(node->down_parsimony_states);
        }
    }
}

void acctran(Tree* tree, size_t num_annotations) {
    /**
     * ACCTRAN (accelerated transformation) (Farris, 1970) aims at reducing the number of ambiguities
     * in the parsimonious result. ACCTRAN forces the state changes to be performed as close to the root as possible,
     * and therefore prioritises the reverse mutations.
     *
     * if N is not a tip:
     *     L, R <- left and right children of N
     *     if intersection(S(N), S(L)) is not empty:
     *         S(L) <- intersection(S(N), S(L))
     *     if intersection(S(N), S(R)) is not empty:
     *         S(R) <- intersection(S(N), S(R))
     *     ACCTRAN(L)
     *     ACCTRAN(R)
     */
    _acctran_node(tree->root, tree->root, num_annotations);
}

void deltran(Tree* tree, size_t num_annotations) {
    /**
     * DELTRAN (delayed transformation) (Swofford & Maddison, 1987) aims at reducing the number of ambiguities
     * in the parsimonious result. DELTRAN makes the changes as close as possible to the leaves,
     * hence prioritizing parallel mutations. DELTRAN is performed after DOWNPASS.
     *
     * if N is not a root:
     *     P <- parent(N)
     *     if intersection(S(N), S(P)) is not empty:
     *         S(N) <- intersection(S(N), S(P))
     * if N is not a tip:
     *     L, R <- left and right children of N
     *     DELTRAN(L)
     *     DELTRAN(R)
     */
    _deltran_node(tree->root, tree->root, num_annotations);
}

void parsimony(Tree *tree, size_t num_annotations, char* method) {
    /**
     * Calculates parsimonious states on the tree and stores them in the node->parsimony_states mask arrays.
     *
     * We use a bit mask for node states to speed up state intersection and union operations.
     * A set of states {0, 2, 3} would look like a mask ...00001101,
     * a set of states {1, 2} would look like ...00000110.
     *
     * The bit mask is stored in node->parsimony_states.
     * As a long type in C is guaranteed to be 64 bits,
     * we need to add longs to our masks as the number of states exceeds 64.
     * Therefore node->parsimony_states is an array, where the states 0-63 are stored in the first element,
     * the states 64-127 are stored in the second element, etc.
     */
    uppass(tree, num_annotations);
    if (strcmp(method, ACCTRAN) == 0) {
        acctran(tree, num_annotations);
    } else {
        downpass(tree, num_annotations);
        if (strcmp(method, DELTRAN) == 0) {
            deltran(tree, num_annotations);
        }
    }
}

void select_parsimonious_states(Tree *tree, size_t num_annotations) {
    size_t i, j;
    Node *node;
    long cur_state;

    for (j = 0; j < tree->nb_nodes; j++) {
        node = tree->nodes[j];
        // initialize the bit mask.
        cur_state = 1;
        for (i = 0; i < num_annotations; i++) {
            if (i % 64 == 0) {
                cur_state = 1;
            } else {
                cur_state <<= 1;
            }
            if ((node->parsimony_states[i / 64] & cur_state) != 0l) {
                node->result_probs[i] = 1.0;
            } else {
                node->result_probs[i] = 0.0;
            }
        }
        free(node->parsimony_states);
    }
}