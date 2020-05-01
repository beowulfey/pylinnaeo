#include <Python.h>
#include <clustal-omega.h>

static PyObject *
clustalo_clustalo(PyObject *self, PyObject *args, PyObject *keywds)
{
    mseq_t *prMSeq = NULL;

    LogDefaultSetup(&rLog);
    LogMuteAll(&rLog);

    // Initialize sequence and alignment options.
    NewMSeq(&prMSeq);
    // Assuming input sequences are not aligned.
    prMSeq->aligned = false;

    opts_t rAlnOpts;
    SetDefaultAlnOpts(&rAlnOpts);

    // Required
    PyObject *inputDict;
    // Optional
    int seqtype = SEQTYPE_PROTEIN;
    PyObject *mbedGuideTree = NULL;
    PyObject *mbedIteration = NULL;
    int numCombinedIterations = rAlnOpts.iNumIterations;
    int maxGuidetreeIterations = rAlnOpts.iMaxGuidetreeIterations;
    int maxHMMIterations = rAlnOpts.iMaxHMMIterations;
    int numThreads = 1;
    static char *kwlist[] = {
        "seqs",
        "seqtype",
        "mbed_guide_tree",
        "mbed_iteration",
        "num_combined_iterations",
        "max_guidetree_iterations",
        "max_hmm_iterations",
        "num_threads",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O!|iOOiiii", kwlist,
            &PyDict_Type, &inputDict,
            &seqtype,
            &mbedGuideTree,
            &mbedIteration,
            &numCombinedIterations,
            &maxGuidetreeIterations,
            &maxHMMIterations,
            &numThreads))
        return NULL;

    if (PyObject_Not(inputDict))
        return PyDict_New();

    InitClustalOmega(numThreads);

    switch (seqtype) {
        case SEQTYPE_DNA:
        case SEQTYPE_RNA:
        case SEQTYPE_PROTEIN:
            prMSeq->seqtype = seqtype;
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "Bad seqtype, must be one of SEQTYPE_DNA, SEQTYPE_RNA, SEQTYPE_PROTEIN, or SEQTYPE_UNKNOWN");
            return NULL;
    }
    if (mbedGuideTree != NULL)
        rAlnOpts.bUseMbed = PyObject_IsTrue(mbedGuideTree);
    if (mbedIteration != NULL)
        rAlnOpts.bUseMbedForIteration = PyObject_IsTrue(mbedIteration);
    rAlnOpts.iNumIterations = numCombinedIterations;
    rAlnOpts.iMaxGuidetreeIterations = maxGuidetreeIterations;
    rAlnOpts.iMaxHMMIterations = maxHMMIterations;

    // Read in sequences from input.
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(inputDict, &pos, &key, &value)) {
        const char *seq = PyUnicode_AsUTF8(value);
	// Sanitize sequence.
	// THIS IS SUPER DANGEROUS AND WILL OVERWRITE YOUR DATA IN PYTHON!
	// FUCKING FIX THIS ASAP
	int seqPos;
        for (seqPos = 0; seqPos < (int)strlen(seq); seqPos++) {
            char *res = &(seq[seqPos]);
            if (isgap(*res))
                continue;
            if (prMSeq->seqtype == SEQTYPE_PROTEIN && strchr(AMINO_ALPHABET, toupper(*res)) == NULL) {
                *res = AMINOACID_ANY;
            } else if (prMSeq->seqtype==SEQTYPE_DNA && strchr(DNA_ALPHABET, toupper(*res)) == NULL) {
                *res = NUCLEOTIDE_ANY;
            } else if (prMSeq->seqtype==SEQTYPE_RNA && strchr(RNA_ALPHABET, toupper(*res)) == NULL) {
                *res = NUCLEOTIDE_ANY;
            }
        }
        AddSeq(&prMSeq, PyUnicode_AsUTF8(key), seq);
    }
    // Can't align with only 1 sequence.
    if (prMSeq->nseqs <= 1) {
        return PyDict_Copy(inputDict);
    }

    // Perform the alignment.
    int rv;
    Py_BEGIN_ALLOW_THREADS
    rv = Align(prMSeq, NULL, &rAlnOpts);
    Py_END_ALLOW_THREADS

    if (rv) {
        PyErr_SetString(PyExc_RuntimeError, "Error while running clustal omega alignment.");
        return NULL;
    }

    // Return the aligned results in a dict.
    PyObject *returnDict = PyDict_New();
    int idx;
    for (idx = 0; idx < prMSeq->nseqs; idx++) {
        const char *key = prMSeq->sqinfo[idx].name;
        PyObject *value = PyUnicode_FromString(prMSeq->seq[idx]);
        PyDict_SetItemString(returnDict, key, value);
    }
    return returnDict;
}

static PyMethodDef ClustaloMethods[] = {
    {"clustalo",  (PyCFunction)clustalo_clustalo, METH_VARARGS | METH_KEYWORDS,
     "Runs clustal omega.\n"
     "Args:\n"
     "  data (dict): dictionary of sequence_name => bases\n"
     "\n"
     "Kwargs:\n"
     "  seqtype (int): clustalo.DNA (1), clustalo.RNA (2) or clustalo.PROTEIN (3)\n"
     "  mbed_guide_tree (bool): whether mBed-like clustering guide tree should be used\n"
     "  mbed_iteration (bool): whether mBed-like clustering iteration should be used\n"
     "  num_combined_iterations (int): number of (combined guide-tree/HMM) iterations\n"
     "  max_guidetree_iterations (int): max guide tree iterations within combined iterations\n"
     "  max_hmm_iterations (int): max HMM iterations within combined iterations\n"
     "  num_threads (int): number of threads to use (requires libclustalo compiled with OpenMP)\n"
     "\n"
     "Returns dict of sequence_named => aligned_bases ('_' for gaps)\n"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef clustalo =
{
	PyModuleDef_HEAD_INIT,
	"clustalo",
	"",
	-1,
	ClustaloMethods
};

PyMODINIT_FUNC PyInit_clustalo(void)
{
    PyObject *module = PyModule_Create(&clustalo);
    PyModule_AddIntConstant(module, "DNA", SEQTYPE_DNA);
    PyModule_AddIntConstant(module, "RNA", SEQTYPE_RNA);
    PyModule_AddIntConstant(module, "PROTEIN", SEQTYPE_PROTEIN);
    return module;
}
