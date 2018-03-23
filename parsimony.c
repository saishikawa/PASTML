#include "pastml.h"


int set_state_to_neighbour_state_intersection(Node *node, size_t first_index, size_t mask_array_len, int include_self) {
    /**
     * Calculates child state intersection and if it's non-empty sets the node->parsimony_states to it.
     * Returns TRUE is the intersection is non-empty else FALSE.
     */
    size_t i, k;
    Node *child;
    int intersection_is_empty = TRUE;

    long *state_intersection = calloc(mask_array_len, sizeof(long));
    if (TRUE == include_self) {
        // set intersection mask to current node mask if it should be included
        memcpy(state_intersection, node->parsimony_states, mask_array_len * sizeof(long));
    } else {
        // set all intersection bits to 1
        for (k = 0; k < mask_array_len; k++) {
            state_intersection[k] = ~0l;
        }
    }

    for (i = first_index; i < node->nb_neigh; i++) {
        child = node->neigh[i];
        intersection_is_empty = TRUE;
        for (k = 0; k < mask_array_len; k++) {
            state_intersection[k] &= child->parsimony_states[k];
            if (state_intersection[k] != 0l) {
                intersection_is_empty = FALSE;
            }
        }
        // if it's already empty no need to go further
        if (intersection_is_empty == TRUE) {
            break;
        }
    }
    if (intersection_is_empty == FALSE) {
        memcpy(node->parsimony_states, state_intersection, mask_array_len * sizeof(long));
    }
    free(state_intersection);
    return intersection_is_empty ? FALSE: TRUE;
}

void set_state_to_most_common_neighbour_states(Node *node, size_t first_index, size_t num_annotations) {
    size_t i, k;
    Node *child;
    long cur_state = 1;
    size_t max_count = 0;
    size_t *state2count = calloc(num_annotations, sizeof(size_t));

    for (k = first_index; k < node->nb_neigh; k++) {
        child = node->neigh[k];
        cur_state = 1;
        for (i = 0; i < num_annotations; i++) {
            if (i % 64 == 0) {
                cur_state = 1;
            } else {
                cur_state <<= 1;
            }
            if ((child->parsimony_states[i / 64] & cur_state) != 0l) {
                state2count[i] += 1;
                if (state2count[i] > max_count) {
                    max_count = state2count[i];
                }
            }
        }
    }
    for (i = 0; i < num_annotations; i++) {
        if (i % 64 == 0) {
            cur_state = 1;
        } else {
            cur_state <<= 1;
        }
        // if this state is (one of) the most frequent, set the corresponding bit to 1, else to 0
        if (state2count[i] == max_count) {
            node->parsimony_states[i / 64] |= cur_state;
        } else {
            node->parsimony_states[i / 64] &= ~cur_state;
        }
    }
    free(state2count);
}


void _uppass_node(Node* node, Node* root, size_t num_annotations) {
    /**
     * We use a bit mask for node states to speed up state intersection and union operations.
     * A set of states {0, 2, 3} would look like a mask ...00001101,
     * a set of states {1, 2} would look like ...00000110.
     *
     * The bit mask is stored in node->parsimony_states.
     * As a long in C is guaranteed to be 64 bits, we need to add longs to our masks as the number of states exceeds 64.
     * Therefore node->parsimony_states is an array, where the states 0-63 are stored in the first element,
     * the states 64-127 are stored in the second element, etc.
     */
    size_t i;
    long cur_state;
    size_t mask_array_len;

    // a tip
    if (node->nb_neigh == 1) {
        // initialize the bit mask for the tips: put 1 in the position of its state(s).
        cur_state = 1;
        for (i = 0; i < num_annotations; i++) {
            if (i % 64 == 0) {
                cur_state = 1;
            } else {
                cur_state <<= 1;
            }
            if (node->bottom_up_likelihood[i] == 1.0) {
                node->parsimony_states[i / 64] |= cur_state;
            }
        }
        return;
    }
    // process children
    for (i = (node == root) ? 0: 1; i < node->nb_neigh; i++) {
        _uppass_node(node->neigh[i], root, num_annotations);
    }
    mask_array_len = (size_t) ceil(num_annotations / 64.0);
    /* quickly calculate state intersection for all children.
     * if it's non-empty, the method sets our state to it,
     * otherwise we need to find the most common states. */
    if (FALSE == set_state_to_neighbour_state_intersection(node, (node == root) ? 0 : 1, mask_array_len, FALSE)) {
        set_state_to_most_common_neighbour_states(node, (node == root) ? 0 : 1, num_annotations);
    }
}

void _downpass_node(Node* node, Node* root, size_t num_annotations) {
    /**
     * We use a bit mask for node states to speed up state intersection and union operations.
     * A set of states {0, 2, 3} would look like a mask ...00001101,
     * a set of states {1, 2} would look like ...00000110.
     *
     * The bit mask is stored in node->parsimony_states.
     * As a long in C is guaranteed to be 64 bits, we need to add longs to our masks as the number of states exceeds 64.
     * Therefore node->parsimony_states is an array, where the states 0-63 are stored in the first element,
     * the states 64-127 are stored in the second element, etc.
     */
    size_t i;
    size_t mask_array_len = (size_t) ceil(num_annotations / 64.0);

    // a tip
    if (node->nb_neigh == 1) {
        /* if we were hesitating between several states for this tip before,
         * we might be able to reduce the number of possibilities by intersecting them with the parent states */
        set_state_to_neighbour_state_intersection(node, 0, mask_array_len, TRUE);
    // the root states are already based on all the tree, nothing else to be done for it
    } else if (node != root) {
        /* quickly calculate state intersection for all the neighbours.
         * if it's non-empty, the method sets our state to it,
         * otherwise we need to find the most common states. */
        if (FALSE == set_state_to_neighbour_state_intersection(node, 0, mask_array_len, FALSE)) {
            set_state_to_most_common_neighbour_states(node, 0, num_annotations);
        }
    }
    // process children
    for (i = (node == root) ? 0: 1; i < node->nb_neigh; i++) {
        _downpass_node(node->neigh[i], root, num_annotations);
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
     * and for each node combines the state information from its children and its parent.
     * As the root state was already the most parsimonious after the UPPASS,
     * we skip it and start directly with the root children.
     *
     * if N is not a tip:
     *     L, R <- left and right children of N
     *     if N is not the root:
     *        P <- parent of N
     *        V <- intersection(S(L), S(R), S(P))
     *        if V is empty:
     *            V <- union(intersection(S(L), S(R)), intersection(S(L), S(P)), intersection(S(P), S(R)))
     *        if V is empty:
     *            V <- union(S(L), S(R), S(P))
     *        S(N) <- V
     *     DOWNPASS(L)
     *     DOWNPASS(R)
     */
    _downpass_node(tree->root, tree->root, num_annotations);
}

void parsimony(Tree *tree, size_t num_annotations) {
    uppass(tree, num_annotations);
    downpass(tree, num_annotations);
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
    }
}