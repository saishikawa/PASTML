#include "pastml.h"
#include "runpastml.h"
#include <getopt.h>
#include <errno.h>

extern QUIET;

int main(int argc, char **argv) {
    char *model = JC, *prob_method = MARGINAL_APPROXIMATION;
    char *annotation_name = NULL;
    char *tree_name = NULL;
    char *out_annotation_name = NULL;
    char *out_tree_name = NULL;
    struct timespec;
    int opt;
    char *arg_error_string = malloc(sizeof(char) * 1024);

    opterr = 0;

    const char *help_string = "usage: PASTML -a ANNOTATION_FILE -t TREE_NWK [-m MODEL] "
            "[-o OUTPUT_ANNOTATION_FILE] [-n OUTPUT_TREE_NWK] [-p PROBABILITY_CALCULATION_METHOD] [-q]\n"
            "\n"
            "required arguments:\n"
            "   -a ANNOTATION_FILE                  path to the annotation csv file containing tip states\n"
            "   -t TREE_NWK                         path to the tree file (in newick format)\n"
            "\n"
            "optional arguments:\n"
            "   -o OUTPUT_ANNOTATION_FILE           path where the output annotation csv file containing node states will be created\n"
            "   -n OUTPUT_TREE_NWK                  path where the output tree file will be created (in newick format)\n"
            "   -m MODEL                            state evolution model (\"JC\" (default) or \"F81\")\n"
            "   -p PROBABILITY_CALCULATION_METHOD   probability calculation method (\"marginal_approx\" (default), \"marginal\", \"max_posteriori\", or \"joint\")\n"
            "   -q                                  quiet, do not print progress information\n";

    opt = getopt(argc, argv, "a:t:o:m:n:p:q");
    do {
        switch (opt) {
            case -1:
                printf(help_string);
                return EINVAL;

            case 'a':
                annotation_name = optarg;
                break;

            case 'o':
                out_annotation_name = optarg;
                break;

            case 'n':
                out_tree_name = optarg;
                break;

            case 't':
                tree_name = optarg;
                break;

            case 'm':
                model = optarg;
                break;

            case 'q':
                QUIET = TRUE;
                break;

            case 'p':
                prob_method = optarg;
                break;

            default: /* '?' */
                snprintf(arg_error_string, 1024, "Unknown arguments...\n\n%s", help_string);
                printf(arg_error_string);
                free(arg_error_string);
                return EINVAL;
        }
    } while ((opt = getopt(argc, argv, "a:t:o:m:n:p:q")) != -1);
    /* Make sure that the required arguments are set correctly */
    if (annotation_name == NULL) {
        snprintf(arg_error_string, 1024, "Annotation file (-a) must be specified.\n\n%s", help_string);
        printf(arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if (tree_name == NULL) {
        snprintf(arg_error_string, 1024, "Tree file (-t) must be specified.\n\n%s", help_string);
        printf(arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if ((strcmp(model, JC) != 0) && (strcmp(model, F81) != 0)) {
        snprintf(arg_error_string, 1024, "Model (-m) must be either %s or %s.\n\n%s", JC, F81, help_string);
        printf(arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if ((strcmp(prob_method, MARGINAL) != 0) && (strcmp(prob_method, MARGINAL_APPROXIMATION) != 0)
        && (strcmp(prob_method, MAX_POSTERIORI) != 0) && (strcmp(prob_method, JOINT) != 0)) {
        snprintf(arg_error_string, 1024, "Probability calculation method (-p) must be one of the following: %s, %s, %s, %s.\n\n%s",
                 MARGINAL_APPROXIMATION, MARGINAL, MAX_POSTERIORI, JOINT, help_string);
        printf(arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    /* No error in arguments */
    free(arg_error_string);

    if (out_annotation_name == NULL) {
        out_annotation_name = calloc(256, sizeof(char));
        sprintf(out_annotation_name, "%s.%s.%s.pastml.out.csv", annotation_name, model, prob_method);
    }
    if (out_tree_name == NULL) {
        out_tree_name = calloc(256, sizeof(char));
        sprintf(out_tree_name, "%s.%s.%s.pastml.out.nwk", tree_name, model, prob_method);
    }
    return runpastml(annotation_name, tree_name, out_annotation_name, out_tree_name, model, prob_method);
}
