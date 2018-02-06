#include "pastml.h"

int dir_a_to_b(Node *a, Node *b) {
    /* this returns the direction from a to b when a and b are two neighbours, otherwise return -1 */
    int i, n = a->nneigh;
    for (i = 0; i < n; i++) if (a->neigh[i] == b) break;
    if (i < n) return i;
    else {
        return -1;
    }
} /* end dir_a_to_b */

int write_subtree_to_stream(Node *node, Node *node_from, FILE *stream, double epsilon, double scaling) {
    int i, direction_to_exclude, n = node->nneigh;
    if (node == NULL || node_from == NULL) return EXIT_SUCCESS;

    if (n == 1) {
        /* terminal node */
        if (node->br[0]->brlen == 0.0) {
            fprintf(stream, "%s:%f", (node->name ? node->name : ""),
                    node->br[0]->brlen + epsilon); /* distance to father */
        } else {
            fprintf(stream, "%s:%f", (node->name ? node->name : ""),
                    node->br[0]->brlen * scaling); /* distance to father */
        }
    } else {
        direction_to_exclude = dir_a_to_b(node, node_from);
        if (-1 == direction_to_exclude) {
            fprintf(stderr, "Fatal error : nodes are not neighbours.\n");
            return EXIT_FAILURE;
        }

        putc('(', stream);
        /* we have to write (n-1) subtrees in total. The last print is not followed by a comma */
        for (i = 1; i < n - 1; i++) {
            if (EXIT_SUCCESS != write_subtree_to_stream(node->neigh[(direction_to_exclude + i) % n], node, stream, epsilon,
                                    scaling)) {
                return EXIT_FAILURE;
            } /* a son */
            putc(',', stream);
        }
        if (EXIT_SUCCESS != write_subtree_to_stream(node->neigh[(direction_to_exclude + i) % n], node, stream, epsilon,
                                scaling)) {
            return EXIT_FAILURE;
        } /* last son */
        putc(')', stream);
        if (node->br[0]->brlen == 0.0) {
            fprintf(stream, "%s:%f", (node->name ? node->name : ""),
                    node->br[0]->brlen + epsilon); /* distance to father */
        } else {
            fprintf(stream, "%s:%f", (node->name ? node->name : ""),
                    node->br[0]->brlen * scaling); /* distance to father */
        }
    }
    return EXIT_SUCCESS;

} /* end write_subtree_to_stream */

int write_nh_tree(Tree *tree, FILE *stream, double epsilon, double scaling) {
    /* writing the tree from the current position in the stream */
    if (!tree) return EXIT_SUCCESS;
    Node *node = tree->node0; /* root or pseudoroot node */
    int i, n = node->nneigh;
    putc('(', stream);
    for (i = 0; i < n - 1; i++) {
        if (EXIT_SUCCESS != write_subtree_to_stream(node->neigh[i], node, stream, epsilon, scaling)) {
            return EXIT_FAILURE;
        } /* a son */
        putc(',', stream);
    }
    if (EXIT_SUCCESS != write_subtree_to_stream(node->neigh[i], node, stream, epsilon, scaling)) {
        return EXIT_FAILURE;
    } /* last son */
    putc(')', stream);

    if (node->name) fprintf(stream, "%s", node->name);
    /* terminate with a semicol AND and end of line */
    putc(';', stream);
    putc('\n', stream);
    return EXIT_SUCCESS;
}
