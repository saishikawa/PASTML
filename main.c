#include "asrml.h"
#include "runpastml.h"
#include <getopt.h>
#include <errno.h>

void *check_alloc(int nbrelt, int sizelt){
    void *retval;
    if( (retval=calloc(nbrelt,sizelt)) != NULL ) {
        return retval;
    }
    sprintf(stderr, "Not enough memory\n");
    return NULL;
}

int main(int argc, char **argv) {
    char *model = "JC";
    char *annotation_name = NULL;
    char *tree_name = NULL;
    char *out_annotation_name = NULL;
    char *out_tree_name = NULL;
    int i, argnum;
    struct timespec;
    double *frequency;
    int opt;
    int check_freq = 0;
    char *scaling = "T", *keep_ID = "T";
    char *arg_error_string = malloc(sizeof(char) * 1024);

    opterr = 0;

    const char *help_string = "usage: PASTML -a ANNOTATION_FILE -t TREE_NWK [-m MODEL] "
            "[-o OUTPUT_ANNOTATION_FILE] [-n OUTPUT_TREE_NWK] "
            "[-s SCALING_ON_OFF] [-I KEEP_INTERNAL_NODE_IDS_ON_OFF] [-f USER_DEFINED_CHAR_FREQUENCES]\n"
            "\n"
            "required arguments:\n"
            "   -a ANNOTATION_FILE                  path to the annotation csv file containing tip states\n"
            "   -t TREE_NWK                         path to the tree file (in newick format)\n"
            "\n"
            "optional arguments:\n"
            "   -o OUTPUT_ANNOTATION_FILE           path where the output annotation csv file containing node states will be created\n"
            "   -n OUTPUT_TREE_NWK                  path where the output tree file will be created (in newick format)\n"
            "   -m MODEL                            state evolution model (JC or F81)\n"
            "   -f USER_DEFINED_FREQUENCY           state frequencies, the sum must be 1.0 (e.g. -f 0.1 0.2 0.3 ... 0.4)\n"
            "   -s SCALING_ON_OFF                   branch length scaling on (T, by default) or off (F)\n"
            "   -I KEEP_INTERNAL_NODE_IDS_ON_OFF    keep internal node ids from the input tree: T (true) or F (false)\n";

    opt = getopt(argc, argv, "a:t:o:m:n:s:f:I:");
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

            case 'f':
                argnum = atoi(argv[optind-1]);
                frequency = check_alloc(argnum, sizeof(double));
                if (frequency == NULL) {
                    return ENOMEM;
                }
                for (i = 1; i <= argnum; i++) {
                    frequency[i - 1] = atof(argv[optind - 1 + i]);
                    printf("Input Frequency No.%d = %lf\n", i, frequency[i]);
                }
                check_freq = 1;
                break;

            case 's':
                scaling=optarg;
                break;

            case 'I':
                keep_ID=optarg;
                break;

            default: /* '?' */
                snprintf(arg_error_string, 1024, "%s%s", "Unknown arguments...\n\n", help_string);
                printf(arg_error_string);
                free(arg_error_string);
                return EINVAL;
        }
    } while ((opt = getopt(argc, argv, "a:t:o:m:n:s:f:I:")) != -1);

    /* Make sure that the required arguments are set correctly */
    if (annotation_name == NULL) {
        snprintf(arg_error_string, 1024, "%s%s", "Annotation file (-a) must be specified.\n\n", help_string);
        printf(arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if (tree_name == NULL) {
        snprintf(arg_error_string, 1024, "%s%s", "Tree file (-t) must be specified.\n\n", help_string);
        printf(arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    if ((strcmp(model, "JC") != 0) && (strcmp(model, "F81") != 0)) {
        snprintf(arg_error_string, 1024, "%s%s", "Model (-m) must be either JC or F81.\n\n", help_string);
        printf(arg_error_string);
        free(arg_error_string);
        return EINVAL;
    }
    /* No error in arguments */
    free(arg_error_string);

    /* Initialise optional arguments is needed */
    if (check_freq == 0) {
        frequency = calloc(MAXCHAR, sizeof(double));
        if (frequency == NULL) {
            return ENOMEM;
        }
    }
    if (out_annotation_name == NULL) {
        out_annotation_name = calloc(256, sizeof(char));
        sprintf(out_annotation_name, "%s.pastml.out.csv", annotation_name);
    }
    if (out_tree_name == NULL) {
        out_tree_name = calloc(256, sizeof(char));
        sprintf(out_tree_name, "%s.pastml.out.nwk", tree_name);
    }
    int res = runpastml(annotation_name, tree_name, out_annotation_name, out_tree_name, model, frequency, scaling, keep_ID);

    free(frequency);
    return res;
}
