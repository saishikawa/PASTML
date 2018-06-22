#include <Python.h>
#include "runpastml.h"
#include "pastml.h"

extern int QUIET;

/*  wrapped pastml function */
static PyObject *infer_ancestral_states(PyObject *self, PyObject *args) {
    char *annotation_name;
    char *tree_name;
    char *out_annotation_name = NULL;
    char *out_tree_name = NULL;
    char *out_param_name = NULL;
    char *in_param_name = NULL;
    char *model = JC;
    char *prob_method = MARGINAL_APPROXIMATION;
    int quiet = FALSE;
    int sts;

    if (!PyArg_ParseTuple(args, "ss|ssssssi", &annotation_name, &tree_name, &in_param_name, &out_annotation_name,
                          &out_tree_name, &out_param_name, &model, &prob_method, &quiet)) {
        printf("Could not parse parameters O_o\n");
        return NULL;
    }
    if (quiet != FALSE) {
        QUIET = TRUE;
    }
    if (strcmp(in_param_name, "") == 0) {
        in_param_name = NULL;
    }
    sts = runpastml(annotation_name, tree_name, out_annotation_name, out_tree_name, out_param_name, model, prob_method,
            in_param_name);
    if (sts != EXIT_SUCCESS) {
        if (errno) {
            return PyErr_SetFromErrno(PyErr_NewException("pastml.error", NULL, NULL));
        } else {
            PyErr_SetString(PyErr_NewException("pastml.error", NULL, NULL), strerror(sts));
            return NULL;
        }
    }
    return PyLong_FromLong(sts);
}

/*  define functions in module */
static PyMethodDef PastmlMethods[] =
        {
                {"infer_ancestral_states", infer_ancestral_states, METH_VARARGS,
                        "Infer tree ancestral states with PASTML.\n"
                        "   :param annotation_file: str, path to the csv file containing two (unnamed) columns: tree tip ids and their states.\n"
                        "   :param tree_file: str, path to the tree in newick format.\n"
                        "   :param in_param_file: str, path to the csv file containing parameters.\n"
                        "   :param out_annotation_file: str, path where the csv file with the inferred annotations will be stored.\n"
                        "   :param out_tree_file: str, path where the output tree (with named internal nodes) in newick format will be stored.\n"
                        "   :param out_param_file: str, path where the output parameter file in csv format will be stored.\n"
                        "   :param model: str, the model of state evolution: pastml.JC or pastml.F81.\n"
                        "   :param prediction_method: str, ancestral state prediction method: "
                                "pastml.MARGINAL_APPROXIMATION, pastml.MARGINAL, pastml.MAX_POSTERIORI, pastml.JOINT, "
                                "pastml.DOWNPASS, pastml.DELTRAN, or pastml.ACCTRAN.\n"
                        "   :param quiet: int, set to non-zero value to prevent PASTML from printing log information.\n"},
                {NULL, NULL, 0, NULL}
        };

#if PY_MAJOR_VERSION >= 3
/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModPyDem =
{
    PyModuleDef_HEAD_INIT,
    "pastml", "PASTML extension for Python 3",
    -1,
    PastmlMethods
};

PyMODINIT_FUNC
PyInit_pastml(void)
{
    PyObject *m = PyModule_Create(&cModPyDem);
    PyModule_AddStringMacro(m, MARGINAL_APPROXIMATION);
    PyModule_AddStringMacro(m, MARGINAL);
    PyModule_AddStringMacro(m, MAX_POSTERIORI);
    PyModule_AddStringMacro(m, DOWNPASS);
    PyModule_AddStringMacro(m, DELTRAN);
    PyModule_AddStringMacro(m, ACCTRAN);
    PyModule_AddStringMacro(m, JOINT);
    PyModule_AddStringMacro(m, JC);
    PyModule_AddStringMacro(m, F81);
    return m;
}

#else

/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC
initpastml(void) {
    (void) Py_InitModule("pastml", PastmlMethods);
}

#endif