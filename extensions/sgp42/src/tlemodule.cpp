#define PY_SSIZE_CLEAN_T
#include <Python.h>

#include <SGP4.h>
#include <stdio.h>  // sscanf

#define SAT_NAME_LENGTH 24

typedef struct {
	PyObject_HEAD
	elsetrec satrec;
	char name[SAT_NAME_LENGTH];
} tle_t;


static int split_tle_string(const char* tle_string, char name[SAT_NAME_LENGTH], char line1[130], char line2[130]) 
{
    try {
        if (tle_string[0] == '1' && tle_string[1] == ' ') {
            // no name present
            sscanf(tle_string, "[^\n] [^\n]", line1, line2);
        }
        else {
            // name present
            sscanf(tle_string, "[^\n] [^\n] [^\n]", name, line1, line2);
        }
    }
    catch (...) {
        return -1;
    }


    return 0;
}

static int tle_init(tle_t* tle, PyObject* args, PyObject* kwargs) {
    const char* tle_string;

    if (PyArg_ParseTuple(args, "s", tle_string) < 0) {
        PyErr_SetString(PyExc_TypeError, "argument must be str type");
        return -1;
    }

    char name[SAT_NAME_LENGTH];
    char line1[130], line2[130];
    int result = split_tle_string(tle_string, name, line1, line2);
}

static PyTypeObject tle_type = {
	PyVarObject_HEAD_INIT(NULL, 0)
    "tle.TwoLineElement",   /* tp_name */
    sizeof(tle_t),            /* tp_basicsize */
    0,                      /* tp_itemsize */
    0,                      /* tp_dealloc */
    0,                      /* tp_vectorcall_offset */
    0,                      /* tp_getattr */
    0,                      /* tp_setattr */
    0,                      /* tp_as_async */
    0,                      /* tp_repr */
    0,                      /* tp_as_number */
    0,                      /* tp_as_sequence */
    0,                      /* tp_as_mapping */
    0,                      /* tp_hash */
    0,                      /* tp_call */
    0,                      /* tp_str */
    0,                      /* tp_getattro */
    0,                      /* tp_setattro */
    0,                      /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,     /* tp_flags */
    PyDoc_STR(""),          /* tp_doc */
    0,                      /* tp_traverse */
    0,                      /* tp_clear */
    0,                      /* tp_richcompare */
    0,                      /* tp_weaklistoffset */
    0,                      /* tp_iter */
    0,                      /* tp_iternext */
    0,                      /* tp_methods */
    0,                      /* tp_members */
    0,                      /* tp_getset */
    0,                      /* tp_base */
    0,                      /* tp_dict */
    0,                      /* tp_descr_get */
    0,                      /* tp_descr_set */
    0,                      /* tp_dictoffset */
    tle_init,                      /* tp_init */
    0,                      /* tp_alloc */
    0,                      /* tp_new */
    0,                      /* tp_free */
};

PyDoc_STRVAR(tle_doc, "Module supporting the SGP4 satellite propagation model.");
static PyModuleDef tle_module = {
	PyModuleDef_HEAD_INIT,	/* m_base */
	"tle",		/* m_name */
	tle_doc,	/* m_doc */
	-1			/* m_size */

	/*PyModuleDef_Base m_base;
  const char* m_name;
  const char* m_doc;
  Py_ssize_t m_size;
  PyMethodDef* m_methods;
  struct PyModuleDef_Slot* m_slots;
  traverseproc m_traverse;
  inquiry m_clear;
  freefunc m_free;*/
};

PyMODINIT_FUNC
PyInit_tle(void)
{
	PyObject* module = NULL;
	
	return module;
}
