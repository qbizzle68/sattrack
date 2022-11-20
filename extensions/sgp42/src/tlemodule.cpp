#define PY_SSIZE_CLEAN_T
#include <Python.h>

#include <SGP4.h>
#include <stdio.h>  // sscanf

#define SAT_NAME_LENGTH 24

//#define _DEBUG

/* forward declaration */
//static PyTypeObject tle_type;

#define TLE_Check(o)    PyObject_TypeCheck(o, &tle_type)

typedef struct {
	PyObject_HEAD
	elsetrec satrec;
	char name[SAT_NAME_LENGTH];
} tle_t;


static int split_tle_string(const char* tle_string, char name[], char line1[], char line2[])
{
    int result;
    try {
        if (tle_string[0] == '1' && tle_string[1] == ' ') {
            // no name present
            result = sscanf(tle_string, "%[^\n] %[^\n]", line1, line2);
            return ((result == 2) - 1);
        }
        else {
            // name present
            result = sscanf(tle_string, "%[^\n] %[^\n] %[^\n]", name, line1, line2);
            return ((result == 3) - 1);
        }
    }
    catch (...) {
        return -1;
    }

    return 0;
}

void print_satrec(const elsetrec& rec) {
    std::cout << "satnum: " << rec.satnum
        << "\nepoch year: " << rec.epochyr
        << "\nepoch num rev: " << rec.epochtynumrev
        << "\nerror: " << rec.error
        << "\nop mode: " << rec.operationmode
        << "\ninit: " << rec.init
        << "\nmethod: " << rec.method
        << "\nisimp: " << rec.isimp
        << "\na: " << rec.a << '\n';
}

static int tle_init(PyObject* self, PyObject* args, PyObject* kwargs) {
    tle_t* tle = (tle_t*)self;
    const char* tle_string = "";
    gravconsttype grav = (gravconsttype)(-1);

    if (PyArg_ParseTuple(args, "s|i", &tle_string, (int*)&grav) < 0) {
        if (!tle_string)
            PyErr_SetString(PyExc_TypeError, "first argument must be str type");
        else if (grav == (gravconsttype)(-1))
            PyErr_SetString(PyExc_TypeError, "second argument must be in type");
        else if (grav != gravconsttype::wgs72old && grav != gravconsttype::wgs72 && grav != gravconsttype::wgs84)
            PyErr_SetString(PyExc_ValueError, "second argument must be wgs72old, wgs72 or wgs84");
        return -1;
    }

    char name[SAT_NAME_LENGTH];
    char line1[130], line2[130];
    int result = split_tle_string(tle_string, name, line1, line2);
    if (result < 0) {
        PyErr_SetString(PyExc_ValueError, "error parsing argument during tle initialization");
        return -1;
    }

#ifdef _DEBUG
    printf("tle string name is :%s:\n", name);
    printf("tle string line1 is :%s:\n", line1);
    printf("tle string line2 is :%s:\n", line2);
#endif

    /* call to init elsetrec values */
    SGP4Funcs::twoline2rv(line1, line2, 'i', grav, tle->satrec);
#ifdef _DEBUG
    printf("error code: %i\n", tle->satrec.error);
#endif
    if (tle->satrec.error) {
        PyObject* exc = PyErr_NewException("tle.sgp4Error", NULL, NULL);
        PyErr_Format(exc, "error in sgp4 library (error = %i)", tle->satrec.error);
        Py_DECREF(exc);
    }

    return 0;
}

static PyObject* get_state(PyObject* self, PyObject* args) {
    tle_t* tle = (tle_t*)self;
    double r[3], v[3];
    PyObject* position_tuple = NULL, * velocity_tuple = NULL, * position = NULL, * velocity = NULL;

//    double jdepoch, jdepochF;
  //  if (PyArgs_ParseTuple(args, "O"))

    SGP4Funcs::sgp4(tle->satrec, 0, r, v);

    if (tle->satrec.error) {
        PyObject* exc = PyErr_NewException("tle.sgp4Error", NULL, NULL);
        PyErr_Format(exc, "error in sgp4 library (error = %i)", tle->satrec.error);
        Py_DECREF(exc);
    }

    PyObject* pyevspace = PyImport_ImportModule("pyevspace");
    if (!pyevspace)
        return NULL;
    PyObject* EVector_type = PyObject_GetAttrString(pyevspace, "EVector");
    Py_DECREF(pyevspace);
    if (!EVector_type)
        return NULL;

    /* build values for EVector constructor */
    position_tuple = Py_BuildValue("(ddd)", r[0], r[1], r[2]);
    if (!position_tuple)
        goto cleanup;
    velocity_tuple = Py_BuildValue("(ddd)", v[0], v[1], v[2]);
    if (!velocity_tuple)
        goto cleanup;

    position = PyObject_CallOneArg(EVector_type, position_tuple);
    if (!position)
        goto cleanup;
    velocity = PyObject_CallOneArg(EVector_type, velocity_tuple);
    if (!velocity)
        goto cleanup;
    Py_DECREF(position_tuple);
    Py_DECREF(velocity_tuple);

    /* add block since goto bypasses initialization here otherwise */
    {
        PyObject* rtn = Py_BuildValue("(OO)", position, velocity);
        Py_DECREF(position);
        Py_DECREF(velocity);
        return rtn;
    }

cleanup:
    Py_XDECREF(position_tuple);
    Py_XDECREF(velocity_tuple);
    Py_XDECREF(position);
    Py_XDECREF(velocity);
    return NULL;
}

static PyMethodDef tle_methods[] = {
    {"getState", get_state, METH_NOARGS, "gets the state at time 0"},
    {NULL, NULL}
};

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
    tle_methods,                      /* tp_methods */
    0,                      /* tp_members */
    0,                      /* tp_getset */
    0,                      /* tp_base */
    0,                      /* tp_dict */
    0,                      /* tp_descr_get */
    0,                      /* tp_descr_set */
    0,                      /* tp_dictoffset */
    (initproc)tle_init,     /* tp_init */
    0,                      /* tp_alloc */
    PyType_GenericNew,                      /* tp_new */
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

    if (PyType_Ready(&tle_type) < 0)
        return NULL;

    module = PyModule_Create(&tle_module);
    if (!module)
        return NULL;

    Py_INCREF(&tle_type);
    if (PyModule_AddObject(module, "TwoLineElement", (PyObject*)&tle_type) < 0) {
        Py_DECREF(module);
        Py_DECREF(&tle_type);
        return NULL;
    }
	
	return module;
}
