#include "pastml.h"
#include "runpastml.h"
#include "logger.h"
#include <getopt.h>
#include <errno.h>

extern int QUIET;

int main(int argc, char **argv) {
    char *model = JC, *prob_method = MARGINAL_APPROXIMATION;
    char *annotation_name = NULL;
    char *tree_name = NULL;
    char *param_name = NULL;
    char *out_annotation_name = NULL;
    char *out_parameter_name = NULL;
    char *out_tree_name = NULL;
    int opt, is_parsimonious;
    char *arg_error_string = malloc(sizeof(char) * 2048);

    opterr = 0;

    const char *help_string = "usage: PASTML -a ANNOTATION_FILE -t TREE_NWK [-m MODEL] "
            "[-o OUTPUT_ANNOTATION_CSV] [-n OUTPUT_TREE_NWK] [-r OUTPUT_PARAMETERS_CSV] [-p PREDICTION_METHOD] [-q]\n"
            "\n"
            "required arguments:\n"
            "   -a ANNOTATION_FILE                  path to the annotation file containing tip states (in csv format)\n"
            "   -t TREE_NWK                         path to the tree file (in newick format)\n"
            "\n"
            "optional arguments:\n"
            "   -o OUTPUT_ANNOTATION_CSV            path where the output annotation file containing node states will be created (in csv format)\n"
            "   -n OUTPUT_TREE_NWK                  path where the output tree file will be created (in newick format)\n"
            "   -r OUTPUT_PARAMETERS_CSV            path where the output parameters file will be created (in csv format)\n"
            "   -i INPUT_PARAMETERS_CSV             path to the parameters file (in csv format)\n"
            "   -m MODEL                            state evolution model for max likelihood prediction methods: "
            "\"JC\" (default) or \"F81\"\n"
            "   -p PREDICTION_METHOD                ancestral state prediction method: \"marginal_approx\" (default), "
            "\"marginal\", \"max_posteriori\", \"joint\", \"downpass\", \"acctran\", or \"deltran\"\n"
            "(\"marginal_approx\", \"marginal\", \"max_posteriori\", and \"joint\" are max likelihood methods, "
            "while \"downpass\", \"acctran\", and \"deltran\" are parsimonious ones)\n"
            "   -q                                  quiet, do not print progress information\n";

    opt = getopt(argc, argv, "a:t:o:r:m:n:i:p:q");
    do {
        switch (opt) {
            case -1:
                printf("%s", help_string);
                return EINVAL;

            case 'a':
                annotation_name = optarg;
                break;

            case 'o':
                out_annotation_name = optarg;
                break;

            case 'i':
                param_name = optarg;
                break;

            case 'r':
                out_parameter_name = optarg;
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
                snprintf(arg_error_string, 2048, "Unknown arguments...\n\n%s", help_string);
                printf("%s", arg_error_string);
                free(arg_error_string);
                return EINVAL;
        }
    } while ((opt = getopt(argc, argv, "a:t:o:r:m:n:i:p:q")) != -1);
    /* Make sure that the required arguments are set correctly */
    if (annotation_name == NULL) {
        snprintf(arg_error_string, 2048, "Annotation file (-a) must be specified.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if (tree_name == NULL) {
        snprintf(arg_error_string, 2048, "Tree file (-t) must be specified.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if (!is_valid_prediction_method(prob_method)) {
        snprintf(arg_error_string, 2048, "Ancestral state prediction method (-p) is not valid.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    is_parsimonious = is_parsimonious_method(prob_method);
    if (!is_parsimonious && !is_valid_model(model)) {
        snprintf(arg_error_string, 2048, "Model (-m) is not valid.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    /* No error in arguments */
    free(arg_error_string);
    return runpastml(annotation_name, tree_name, out_annotation_name, out_tree_name, out_parameter_name,
                     model, prob_method, param_name);
}
