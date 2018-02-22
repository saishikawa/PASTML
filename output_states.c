#include <errno.h>
#include "pastml.h"


void output_states(Node *nd, Node *root, int num_annotations, char **character, FILE *outfile) {
    int i, j, num_nonzero_states = 0;

    // process children
    for (i = (nd == root) ? 0: 1; i < nd->nneigh; i++) {
        output_states(nd->neigh[i], root, num_annotations, character, outfile);
    }

    fprintf(outfile, "%s", nd->name);

    // a tip
    if (nd->nneigh == 1) {
        for (i = 0; i < num_annotations; i++) {
            fprintf(outfile, ", %.5f", (i == nd->pupko_state) ? 1.0 : 0.0);
        }
    // an internal node
    } else {
        for (i = 0; i < num_annotations; i++) {
            if (nd->local_flag[i] == 1) {
                num_nonzero_states++;
            }
        }
        double marginal_val = 1.0 / num_nonzero_states;
        for (i = 0; i < num_annotations; i++) {
            for (j = 0; j < num_annotations; j++) {
                if (strcmp(character[i], character[nd->tmp_best[j]]) == 0) {
                    fprintf(outfile, ", %.5f", (nd->local_flag[j] == 1) ? marginal_val : 0.0);
                }
            }
        }
    }

    fprintf(outfile, "\n");
}

int output_state_ancestral_states(Node *root, int num_annotations, char **character, char *output_filepath) {
    FILE* outfile = fopen(output_filepath, "w");
    if (!outfile) {
        fprintf(stderr, "Output annotation file %s is impossible to access.", output_filepath);
        fprintf(stderr, "Value of errno: %d\n", errno);
        fprintf(stderr, "Error opening the file: %s\n", strerror(errno));
        return ENOENT;
    }

    // print the header
    fprintf(outfile, "node ID");
    for (int i = 0; i < num_annotations; i++) {
        fprintf(outfile, ", %s", character[i]);
    }
    fprintf(outfile, "\n");

    output_states(root, root, num_annotations, character, outfile);
    fclose(outfile);
    return EXIT_SUCCESS;
}
