#include "pastml.h"
#include "runpastml.h"
#include "logger.h"
#include <getopt.h>
#include <errno.h>

static const int HELP_STRING_LEN = 4000;
extern bool QUIET;

int main(int argc, char **argv) {
    char *model = F81, *prob_method = MARGINAL_APPROXIMATION;
    char *annotation_name = NULL;
    char *tree_name = NULL;
    char *param_name = NULL;
    char *out_annotation_name = NULL;
    char *out_parameter_name = NULL;
    char *out_mp_name = NULL;
    char *out_tree_name = NULL;
    int opt, is_parsimonious;
    char *arg_error_string = malloc(sizeof(char) * HELP_STRING_LEN);

    opterr = 0;

    const char *help_string = "usage: PASTML -a ANNOTATION_FILE -t TREE_NWK [-o OUTPUT_ANNOTATION_CSV] [-n OUTPUT_TREE_NWK] [-q] "
                              "[-p PREDICTION_METHOD] [-m MODEL] [-r OUTPUT_PARAMETERS_CSV] [-i INPUT_PARAMETERS_CSV] "
                              "[-f OUTPUT_MARGINAL_PROBAB_CSV]\n"
            "\n"
            "required arguments:\n"
            "   -a ANNOTATION_FILE                  path to the annotation file containing tip states (in csv format: <tip_id>,<state>);"
            "                                       should not contain any header, must include all tips (leave <state> blank when unknown)\n"
            "   -t TREE_NWK                         path to the rooted tree file (in newick format)\n"
            "\n"
            "optional arguments:\n"
            "   -o OUTPUT_ANNOTATION_CSV            path where the output annotation file containing node states will be created (in csv format),\n"
            "                                       if not specified, it will be generated by adding a suffix to the input annotation file name\n"
            "   -n OUTPUT_TREE_NWK                  if specified, the tree with named internal nodes will be saved to this file (in newick format)\n"
            "   -q                                  quiet, do not print progress information\n"
            "   -p PREDICTION_METHOD                ancestral state prediction method;\n"
            "                                       can be one of the ML methods: marginal_approx (default, marginal posterior probabilities approximation (MPPA)),\n"
            "                                       marginal, max_posteriori (maximum a posteriori MAP), joint;\n"
            "                                       or one of the parsimonious methods: downpass, acctran (accelerated transformation), or deltran (delayed transformation)\n"
            "\n"
            "optional arguments for ML methods only:\n"
            "   -m MODEL                            state evolution model: F81 (default, Felsenstein 1981-like, character frequencies are optimised, or supplied via -i option),\n"
            "                                       JC (Jukes-Cantor-like, all character frequencies and change rates are equal),\n"
            "                                       or EFT (character frequencies are Estimated From the Tip state distribution)\n"
            "   -r OUTPUT_PARAMETERS_CSV            if specified, the model parameters (state frequencies, scaling factor)\n"
            "                                       will be saved to this file (in csv format)\n"
            "   -i INPUT_PARAMETERS_CSV             if specified, the parameters found in this file will be fixed to the corresponding values;\n"
            "                                       the file should be in csv format: <parameter>,<value>\n"
            "                                       where <parameter> can be scaling factor,\n"
            "                                       or (for F81 model only) state (<value> should contain its frequency, and all the states should be supplied);\n"
            "                                       all other parameters will be ignored; see the output of -r option for the file format example\n"
            "   -f OUTPUT_MARGINAL_PROBAB_CSV       (for marginal ML methods (MAP, MPPA, marginal) only) if specified, \n"
            "                                       the marginal probabilities of node states will be saved to this file (in csv format).\n";

    while ((opt = getopt(argc, argv, "a:t:o:r:m:n:i:f:p:hq")) != -1) {
        switch (opt) {
            case 'h':
                printf("%s", help_string);
                return EXIT_SUCCESS;

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

            case 'f':
                out_mp_name = optarg;
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
                QUIET = true;
                break;

            case 'p':
                prob_method = optarg;
                break;

            default: /* '?' */
                snprintf(arg_error_string, HELP_STRING_LEN, "Unknown arguments...\n\n%s", help_string);
                printf("%s", arg_error_string);
                free(arg_error_string);
                return EINVAL;
        }
    }
    /* Make sure that the required arguments are set correctly */
    if (annotation_name == NULL) {
        snprintf(arg_error_string, HELP_STRING_LEN, "Annotation file (-a) must be specified.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if (tree_name == NULL) {
        snprintf(arg_error_string, HELP_STRING_LEN, "Tree file (-t) must be specified.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if (!is_valid_prediction_method(prob_method)) {
        snprintf(arg_error_string, HELP_STRING_LEN, "Ancestral state prediction method (-p) is not valid.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    is_parsimonious = is_parsimonious_method(prob_method);
    if (!is_parsimonious && !is_valid_model(model)) {
        snprintf(arg_error_string, HELP_STRING_LEN, "Model (-m) is not valid.\n\n%s", help_string);
        printf("%s", arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    /* No error in arguments */
    free(arg_error_string);
    return runpastml(annotation_name, tree_name, out_annotation_name, out_tree_name, out_parameter_name,
                     model, prob_method, param_name, out_mp_name);
}
