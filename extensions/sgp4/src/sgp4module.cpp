#define PY_SSIZE_CLEAN_T
#include <Python.h>

#include <SGP4.h>
#include <stdio.h>  // sscanf

/**
 * tle.py
 *  - replace TwoLineElement
 * _orbit.py
 *  - replace _compute_eccentric_vector(position, velocity, MU)
 *  - replace _elements_from_state(position, velocity, MU)
 * _topocentric.py
 *  - replace _elements_from_tle(tle, jd)   (compute state at jd and call _elements_from_state)
 *  - Satellite.anomalyAt(jd, anomalyType)  to use _elements_from_tle
 */

#define TLE_Check(o)    PyObject_TypeCheck(o, &tle_type)
#define TLE_JDNUMBER(o)     o->satrec.jdsatepoch
#define TLE_JDFRACTION(o)   o->satrec.jdsatepochF
#define TLE_PERIOD(o)       o->satrec.period_sec

#define SAT_NAME_LENGTH 25  // 24 chars for TLE format plus null-char
#define EARTH_MU        3.986004418e5 // km^3s^-2

/* struct for TwoLineElement type */
typedef struct {
	PyObject_HEAD
	elsetrec satrec;
	char name[SAT_NAME_LENGTH];
} tle_t;

/* struct for module state */
typedef struct {
    PyObject* EVector;
    PyObject* JulianDate;
    PyObject* TwoLineElement;
} sgp4_state;

static void 
trim_name(char name[]) 
{
    char* end = name + 23;
    while (end > name && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
}

static int 
split_tle_string(const char* tle_string, char name[], char line1[], char line2[])
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
            trim_name(name);
            return ((result == 3) - 1);
        }
    }
    catch (...) {
        return -1;
    }

    return 0;
}

static int 
tle_init(tle_t* self, PyObject* args, PyObject* kwargs) 
{
    const char* tle_string = "";
    gravconsttype grav = gravconsttype::wgs72;

    // get arguments, default grav to wgs72
    if (PyArg_ParseTuple(args, "s|i", &tle_string, (int*)&grav) < 0)
        return -1;
    if (grav != gravconsttype::wgs72old && grav != gravconsttype::wgs72 && grav != gravconsttype::wgs84) {
        PyErr_SetString(PyExc_ValueError, "second argument must be wfs72old, wgs72 or wgs84");
        return -1;
    }

    // parse tle string argument into lines
    char name[SAT_NAME_LENGTH];
    char line1[130], line2[130]; // twoline2rv expects arrays of 130 length
    int result = split_tle_string(tle_string, name, line1, line2);
    if (result < 0) {
        PyErr_SetString(PyExc_ValueError, "error parsing argument during tle initialization");
        return -1;
    }
    strncpy(self->name, name, SAT_NAME_LENGTH);

    // call to init elsetrec values
    SGP4Funcs::twoline2rv(line1, line2, 'i', grav, self->satrec);
    if (self->satrec.error) {
        // do we save this somewhere else?
        PyObject* exc = PyErr_NewException("_sgp4.sgp4Error", NULL, NULL);
        PyErr_Format(exc, "error in sgp4 library (error = %i)", self->satrec.error);
        Py_DECREF(exc);
        return -1;
    }

    return 0;
}

static int
jd_AsDouble(PyObject* const jd, PyTypeObject* jdType, double& rtn, const char* paramStr)
{
    if (!Py_IS_TYPE(jd, jdType)) {
        PyErr_Format(PyExc_TypeError, "%s parameter must be JulianDate type", paramStr);
        return -1;
    }

    PyObject* value = PyObject_GetAttrString(jd, "value");
    if (!value)
        return -1;

    double tmp = PyFloat_AsDouble(value);
    if (tmp == -1.0 && PyErr_Occurred())
        return -1;

    rtn = tmp;
    return 0;
}

static PyObject*
build_vector(PyObject* vectorType, double x, double y, double z)
{
    PyObject* tuple = Py_BuildValue("(ddd)", x, y, z);
    if (!tuple)
        return NULL;

    PyObject* vector = PyObject_CallOneArg(vectorType, tuple);
    Py_DECREF(tuple);
    return vector;
}

static PyObject* 
get_state(PyObject* self, PyObject* const* args, Py_ssize_t size) 
{
    if (size != 2) {
        PyErr_Format(PyExc_TypeError, "_sgp4._getState() takes exactly two arguments (%i given)", size);
        return NULL;
    }

    //* get the module state
    sgp4_state* state = (sgp4_state*)PyModule_GetState(self);
    if (!state)
        return NULL;

    // extract arguments
    if (!Py_IS_TYPE(args[0], (PyTypeObject*)state->TwoLineElement)) {
        PyErr_SetString(PyExc_TypeError, "first parameter must be TwoLineElement type");
        return NULL;
    }

    tle_t* tle = (tle_t*)args[0];
    double jd;
    if (jd_AsDouble(args[1], (PyTypeObject*)state->JulianDate, jd, "second") < 0)
        return NULL;

    // run the sgp4 library function
    double r[3], v[3], dt = (jd - TLE_JDNUMBER(tle) - TLE_JDFRACTION(tle)) * 1440.0;
    SGP4Funcs::sgp4(tle->satrec, dt, r, v);

    // check if the error flag was set
    if (tle->satrec.error) {
        PyObject* exc = PyErr_NewException("_sgp4.sgp4Error", NULL, NULL);
        PyErr_SetString(exc, "error in sgp4 library (error = %i)");
        Py_DECREF(exc);
        return NULL;
    }

    // get buffer from new vectors here 

    PyObject* position = build_vector(state->EVector, r[0], r[1], r[2]);
    if (!position)
        return NULL;

    PyObject* velocity = build_vector(state->EVector, v[0], v[1], v[2]);
    if (!velocity) {
        Py_DECREF(position);
        return NULL;
    }

    // build and return state
    PyObject* rtn = Py_BuildValue("(OO)", position, velocity);
    Py_DECREF(position);
    Py_DECREF(velocity);
    return rtn;
}

/* this should be replacable with access to EVector buffer */
static int
get_vec_items(PyObject* const self, double v[])
{
    PyObject* obj = PySequence_Fast(self, NULL);
    if (!obj)
        return -1;
    PyObject** items = PySequence_Fast_ITEMS(obj);
    v[0] = PyFloat_AS_DOUBLE(items[0]);
    v[1] = PyFloat_AS_DOUBLE(items[1]);
    v[2] = PyFloat_AS_DOUBLE(items[2]);
    return 0;
}

static PyObject*
get_elements_from_state(PyObject* self, PyObject* const* args, Py_ssize_t size)
{
    // signiture: _elements_from_state(position: EVector, velocity: EVector, MU: float = EARTH_MU) 
    //                  -> (raan, inc, aop, ecc, a, m, nu)

    double mu = EARTH_MU;
    if (size == 3) {
        mu = PyFloat_AsDouble(args[2]);
        if (mu == -1.0 && PyErr_Occurred())
            return NULL;
    }
    else if (size != 2) {
        PyErr_Format(PyExc_TypeError, "_sgp4._getState() takes two or three arguments (%i given)", size);
        return NULL;
    }

    // get module state
    sgp4_state* state = (sgp4_state*)PyModule_GetState(self);
    if (!state)
        return NULL;

    // extract arguments
    if (!Py_IS_TYPE(args[0], (PyTypeObject*)state->EVector)) {
        PyErr_SetString(PyExc_TypeError, "first argument must be pyevspace.EVector type");
        return NULL;
    }
    if (!Py_IS_TYPE(args[1], (PyTypeObject*)state->EVector)) {
        PyErr_SetString(PyExc_TypeError, "second argument must be pyevspace.EVector type");
        return NULL;
    }

    //  get EVector buffer here

    double r[3], v[3];
    if (get_vec_items(args[0], r) < 0)
        return NULL;
    if (get_vec_items(args[1], v) < 0)
        return NULL;

    double p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;
    SGP4Funcs::rv2coe_SGP4(r, v, mu, p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper);

    // returning mean anomaly right now with this 
    return Py_BuildValue("ddddddd", omega, incl, argp, ecc, a, m, nu);
}

static PyObject*
get_elements_from_tle(PyObject* self, PyObject* const* args, Py_ssize_t size)
{
    // signiture: _elements_from_tle(tle: TwoLineElement, jd: JulianDate, MU: float = EARTH_MU)
    //                  -> (raan, inc, aop, ecc, a, m, nu)

    double mu = EARTH_MU;
    if (size == 3) {
        mu = PyFloat_AsDouble(args[2]);
        if (mu == -1.0 && PyErr_Occurred())
            return NULL;
    }
    else if (size != 2) {
        PyErr_Format(PyExc_TypeError, "_sgp4._getState() takes two or three arguments (%i given)", size);
        return NULL;
    }

    // get module state
    sgp4_state* state = (sgp4_state*)PyModule_GetState(self);
    if (!state)
        return NULL;

    // extract arguments
    if (!Py_IS_TYPE(args[0], (PyTypeObject*)state->TwoLineElement)) {
        PyErr_SetString(PyExc_TypeError, "first argument must be TwoLineElement type");
        return NULL;
    }
    double jd; 
    jd_AsDouble(args[1], (PyTypeObject*)state->JulianDate, jd, "second");

    // get state
    tle_t* tle = (tle_t*)args[0];
    double r[3], v[3], dt = jd - TLE_JDNUMBER(tle) - TLE_JDFRACTION(tle);
    SGP4Funcs::sgp4(tle->satrec, dt, r, v);

    // get elements
    double p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;
    SGP4Funcs::rv2coe_SGP4(r, v, mu, p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper);

    // returning mean anomaly right now with this 
    return Py_BuildValue("ddddddd", omega, incl, argp, ecc, a, m, nu);
}

static PyObject* 
compute_ecc_vector(PyObject* self, PyObject* const* args, Py_ssize_t size)
{
    // signiture: compute_ecc_vector(position: EVector, velocity: EVector, mu: float = EARTH_MU) -> EVector
    double mu = EARTH_MU;
    if (size == 3) {
        mu = PyFloat_AsDouble(args[2]);
        if (mu == -1.0 && PyErr_Occurred())
            return NULL;
    }
    else if (size != 2) {
        PyErr_Format(PyExc_TypeError, "_sgp4._compute_eccentric_vector() takes two or three arguments (%i given)", size);
        return NULL;
    }

    // get module state
    sgp4_state* state = (sgp4_state*)PyModule_GetState(self);
    if (!state)
        return NULL;

    // check arguments
    if (!Py_IS_TYPE(args[0], (PyTypeObject*)state->EVector)) {
        PyErr_SetString(PyExc_TypeError, "first argument must be EVector type");
        return NULL;
    }
    if (!Py_IS_TYPE(args[1], (PyTypeObject*)state->EVector)) {
        PyErr_SetString(PyExc_TypeError, "second argument must be EVector type");
        return NULL;
    }

    // here is a good spot for the ability to get the buffer of the vector args...

    // replicate computing eccentric vector from rv2coe
    double r[3], v[3];
    // replace this with the buffer from EVector
    if (get_vec_items(args[0], r) < 0)
        return NULL;
    if (get_vec_items(args[1], v) < 0)
        return NULL;

    double ebar[3];
    double magr = SGP4Funcs::mag_SGP4(r);
    double magv = SGP4Funcs::mag_SGP4(v);
    double c1 = magv * magv - mu / magr;
    double rdotv = SGP4Funcs::dot_SGP4(r, v);
    for (int i = 0; i < 3; i++)
        ebar[i] = (c1 * r[i] - rdotv * v[i]) / mu;

    // build and return eccentric vector
    PyObject* eccVec = build_vector(state->EVector, ebar[0], ebar[1], ebar[2]);
    return eccVec;
}

static PyMethodDef sgp4_methods[] = {
    {"_getState", (PyCFunction)get_state, METH_FASTCALL, "_getState(tle, jd, mu) -> (position, velocity)"},
    {"_elements_from_state", (PyCFunction)get_elements_from_state, METH_FASTCALL, "_elements_from_state(position, velocity, mu) -> (sma, ecc, inc, raan, aop, m, nu)"},
    {"_elements_from_tle", (PyCFunction)get_elements_from_tle, METH_FASTCALL, "_elements_from_tle(tle, jd, mu) -> (sma, ecc, inc, raan, aop, m, nu)"},
    {"_compute_eccentric_vector", (PyCFunction)compute_ecc_vector, METH_FASTCALL, "_compute_eccentric_vector(position, velocity, mu) -> eccentricVector"},
    {NULL, NULL}
};

static PyObject*
tle_epoch(tle_t* self, void* Py_UNUSED) 
{
    sgp4_state* state = (sgp4_state*)PyType_GetModuleState(Py_TYPE(self));
    if (!state)
        return NULL;

    double jdValue = TLE_JDNUMBER(self) + TLE_JDFRACTION(self);
    PyObject* jd = PyObject_CallMethod(state->JulianDate, "fromNumber", "d", jdValue);

    return jd;
}

static PyObject*
tle_name(tle_t* self, void* Py_UNUSED) 
{
    return PyUnicode_FromString(self->name);
}

static PyObject*
tle_sma(tle_t* self, void* Py_UNUSED)
{
    return PyFloat_FromDouble(self->satrec.a * 6371.0);
}

static PyObject*
tle_periapsis(tle_t* self, void* Py_UNUSED)
{
    return PyFloat_FromDouble((self->satrec.altp + 1.0) * 6371.0);
}

static PyObject*
tle_apoapsis(tle_t* self, void* Py_UNUSED)
{
    return PyFloat_FromDouble((self->satrec.alta + 1.0) * 6371.0);
}

static PyGetSetDef tle_getset[] = {
    {"epoch",       (getter)tle_epoch},
    {"name",        (getter)tle_name},
    {"sma",         (getter)tle_sma},
    {"apoapsis",    (getter)tle_apoapsis},
    {"periapsis",   (getter)tle_periapsis},
    {NULL}
};

static PyType_Slot tle_slots[] = {
    {Py_tp_getset, tle_getset},
    {Py_tp_init, (initproc)tle_init},
    {Py_tp_new, PyType_GenericNew},
    {0, NULL}
};

static PyType_Spec tle_spec = {
    "tle.TwoLineElement",
    sizeof(tle_t),
    0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HEAPTYPE,
    tle_slots
};

static int
sgp4_exec(PyObject* module)
{
    sgp4_state* state = (sgp4_state*)PyModule_GetState(module);
    if (!state)
        return -1;

    PyObject* type = PyType_FromModuleAndSpec(module, &tle_spec, NULL);
    if (!type)
        return -1;
    if (PyModule_AddObjectRef(module, "TwoLineElement", type) < 0) {
        Py_DECREF(type);
        return -1;
    }
    state->TwoLineElement = type;

    PyObject* module_import = PyImport_ImportModule("pyevspace");
    if (!module_import) {
        Py_DECREF(state->TwoLineElement);
        return -1;
    }
    type = PyObject_GetAttrString(module_import, "EVector");
    Py_DECREF(module_import);
    if (!type) {
        Py_DECREF(state->TwoLineElement);
        return -1;
    }
    state->EVector = type;
    type = NULL;

    module_import = PyImport_ImportModule("sattrack.spacetime.juliandate");
    if (!module_import) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->EVector);
        return -1;
    }
    type = PyObject_GetAttrString(module_import, "JulianDate");
    Py_DECREF(module_import);
    if (!type) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->EVector);
        return -1;
    }
    state->JulianDate = type;
    type = NULL;

    if (PyModule_AddIntConstant(module, "WGS72OLD", (int)gravconsttype::wgs72old) < 0) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->EVector);
        Py_DECREF(state->JulianDate);
        return -1;
    }
    if (PyModule_AddIntConstant(module, "WGS72", (int)gravconsttype::wgs72) < 0) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->EVector);
        Py_DECREF(state->JulianDate);
        return -1;
    }
    if (PyModule_AddIntConstant(module, "WGS84", (int)gravconsttype::wgs84) < 0) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->EVector);
        Py_DECREF(state->JulianDate);
        return -1;
    }

    return 0;
}

static void sgp4_free(void* self) 
{
    sgp4_state* state = (sgp4_state*)PyModule_GetState((PyObject*)self);
    if (!state)
        return;

    Py_DECREF(state->EVector);
    Py_DECREF(state->JulianDate);
    Py_DECREF(state->TwoLineElement);
}

static PyModuleDef_Slot sgp4_slots[]{
    {Py_mod_exec, sgp4_exec},
    {0, NULL}
};

PyDoc_STRVAR(tle_doc, "Module supporting the SGP4 satellite propagation model.");
static PyModuleDef tle_module = {
	PyModuleDef_HEAD_INIT,	/* m_base */
	"_sgp4",		        /* m_name */
	tle_doc,	            /* m_doc */
	sizeof(sgp4_state),	    /* m_size */
    sgp4_methods,           /* m_methods */
    sgp4_slots,             /* m_slots */
    NULL,                   /* m_traverse */
    NULL,                   /* m_clear */
    sgp4_free,              /* m_free */
};

PyMODINIT_FUNC
PyInit__sgp4(void)
{
	return PyModuleDef_Init(&tle_module);
}
