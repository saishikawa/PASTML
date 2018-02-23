#include <errno.h>
#include "pastml.h"
#include "lik.h"

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
                write_subtree_to_stream(node->neigh[(direction_to_exclude + i) % n], node, stream, epsilon, scaling)) {
                return EXIT_FAILURE;
            } /* a son */
            putc(',', stream);
        }
        if (EXIT_SUCCESS != write_subtree_to_stream(node->neigh[(direction_to_exclude + i) % n], node, stream, epsilon,
                                                    scaling)) {
            return EXIT_FAILURE;
        } /* last son */
        putc(')', stream);
    }
    // write node's name and dist to father
    fprintf(stream, "%s:%f", (node->name ? node->name : ""), node->brlen);
    return EXIT_SUCCESS;
} /* end write_subtree_to_stream */

int write_nh_tree(Node *root, char *output_filepath, double epsilon, double scaling) {

    FILE* output_file = fopen(output_filepath, "w");
    if (!output_file) {
        fprintf(stderr, "Output tree file %s is impossible to access.", output_filepath);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }
    /* writing the tree from the current position in the file */
    int i, n = root->nneigh;
    putc('(', output_file);
    for (i = 0; i < n - 1; i++) {
        if (EXIT_SUCCESS != write_subtree_to_stream(root->neigh[i], root, output_file, epsilon, scaling)) {
            return EXIT_FAILURE;
        } /* a son */
        putc(',', output_file);
    }
    if (EXIT_SUCCESS != write_subtree_to_stream(root->neigh[i], root, output_file, epsilon, scaling)) {
        return EXIT_FAILURE;
    } /* last son */
    putc(')', output_file);

    if (root->name) {
        fprintf(output_file, "%s", root->name);
    }
    /* terminate with a semicol AND and end of line */
    putc(';', output_file);
    putc('\n', output_file);

    fclose(output_file);
    return EXIT_SUCCESS;
}
