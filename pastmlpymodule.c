//
// Created by azhukova on 1/24/18.
//
#include <Python.h>
#include "runpastml.h"

/*  wrapped pastml function */
static PyObject *infer_ancestral_states(PyObject *self, PyObject *args) {
    char *annotation_name;
    char *tree_name;
    char *out_annotation_name;
    char *out_tree_name;
    char *model;
    int sts;
    double *frequency;

    if (!PyArg_ParseTuple(args, "sssss", &annotation_name, &tree_name, &out_annotation_name, &out_tree_name, &model)) {
        return NULL;
    }

    sts = runpastml(annotation_name, tree_name, out_annotation_name, out_tree_name, model, frequency, "T", "T");
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
                    "   :param out_annotation_file: str, path where the csv file with the inferred annotations will be stored.\n"
                    "   :param out_tree_file: str, path where the output tree (with named internal nodes) in newick format will be stored.\n"
                    "   :param model: str, the model of state evolution, must be either JC or F81.\n"},
    {NULL,        NULL,      0,            NULL}
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
    return PyModule_Create(&cModPyDem);
}

#else

/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC
initpastml(void) {
    (void) Py_InitModule("pastml", PastmlMethods);
}

#endif