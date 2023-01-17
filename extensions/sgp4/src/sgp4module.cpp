#define PY_SSIZE_CLEAN_T
#include <Python.h>
#include <structmember.h> // PyMemberDef

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

#define SAT_NAME_LENGTH     25  // 24 chars for TLE format plus null-char
#define TLE_LINE_LENGTH     70  // TLE line length
#define SGP4_LINE_LENGTH    130 // sgp4 library methods expect lines
                                // with length of 130
#define EARTH_MU            3.986004418e5 // km^3s^-2

/* struct for TwoLineElement type */
typedef struct {
    PyObject_HEAD
    elsetrec satrec;
    //char name[SAT_NAME_LENGTH];
    char* name;
    char* line1;
    char* line2;
} tle_t;

/* struct for module state */
typedef struct {
    PyObject* Vector;
    PyObject* JulianDate;
    PyObject* TwoLineElement;
} sgp4_state;

static void 
trim_name(char *name)
{
    // start from the end of the string
    char* end = name + 23;
    // work back towards the 'front' looking for the first non-space char
    while (end > name && isspace((unsigned char)*end)) {
        end--;
    }
    // make the next char mark the end of the string
    end[1] = '\0';
}

static int 
split_tle_string(const char* tle_string, char *name, char *line1, char *line2)
{
    int result;
    // try block allows us to convert cpp style exceptions to cpython style errors
    try {
        if (tle_string[0] == '1' && tle_string[1] == ' ') {
            // no name present
            result = sscanf(tle_string, "%[^\n] %[^\n]", line1, line2);
            return (result == 2) ? 0 : -1;
        }
        else {
            // name present
            result = sscanf(tle_string, "%[^\n] %[^\n] %[^\n]",
                            name, line1, line2);
            trim_name(name);
            return (result == 3) ? 0 : -1;
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
    if (PyArg_ParseTuple(args, "s|i", &tle_string, (int*)&grav) < 0) {
        return -1;
    }
    if (grav != gravconsttype::wgs72old && grav != gravconsttype::wgs72 &&
        grav != gravconsttype::wgs84) {
        PyErr_SetString(PyExc_ValueError, "second argument "
                        "must be wgs72old, wgs72 or wgs84");
        return -1;
    }

    // parse tle string argument into lines
    self->name = (char*)calloc(SAT_NAME_LENGTH, 1);
    if (!self->name) {
        PyErr_NoMemory();
        return -1;
    }
    char line1[130], line2[130]; // twoline2rv expects arrays of 130 length
    if (split_tle_string(tle_string, self->name, line1, line2) < 0) {
        PyErr_SetString(PyExc_ValueError, "error parsing "
                        "argument during tle initialization");
        return -1;
    }

    self->line1 = (char*)calloc(TLE_LINE_LENGTH, 1);
    if (!self->line1) {
        PyErr_NoMemory();
        return -1;
    }
    strncpy(self->line1, line1, TLE_LINE_LENGTH);
    self->line2 = (char*)calloc(TLE_LINE_LENGTH, 1);
    if (!self->line2) {
        PyErr_NoMemory();
        return -1;
    }
    strncpy(self->line2, line2, TLE_LINE_LENGTH);

    // call to init elsetrec values
    SGP4Funcs::twoline2rv(line1, line2, 'i', 
                          grav, self->satrec);
    if (self->satrec.error) {
        // do we save this somewhere else?
        // todo: make an exception for each type of error or add an attribute
        // with the error value
        PyObject* exc = PyErr_NewException("sgp4.sgp4Error", NULL, NULL);
        PyErr_Format(exc, "error in sgp4 library (error = %i)",
                     self->satrec.error);
        Py_DECREF(exc);
        return -1;
    }

    return 0;
}

static void
tle_free(void* self)
{
    tle_t* tle = (tle_t*)self;
    free(tle->name);
    free(tle->line1);
    free(tle->line2);
}

static int
jd_AsDouble(PyObject* const jd, PyTypeObject* jdType, double& rtn,
            const char* paramStr)
{
    if (!Py_IS_TYPE(jd, jdType)) {
        PyErr_Format(PyExc_TypeError,
                     "%s parameter must be JulianDate type",
                     paramStr);
        return -1;
    }

    PyObject* value = PyObject_GetAttrString(jd, "value");
    if (!value) {
        return -1;
    }

    double tmp = PyFloat_AsDouble(value);
    Py_DECREF(value);
    if (tmp == -1.0 && PyErr_Occurred()) {
        return -1;
    }

    rtn = tmp;
    return 0;
}

static inline PyObject*
build_vector(PyObject* vectorType, double x, double y, double z)
{
    PyObject* vector = PyObject_CallFunction(vectorType, "ddd",
                                             x, y, z);
    return vector;
}

static PyObject* 
get_state(PyObject* self, PyObject* const* args, Py_ssize_t size) 
{
    // _getState(TwoLineElement, JulianDate) -> (Vector, Vector)

    if (size != 2) {
        PyErr_Format(PyExc_TypeError,
                     "_sgp4._getState() takes exactly two arguments "
                     "(%i given)",
                     size);
        return NULL;
    }

    // get the module state
    sgp4_state* state = (sgp4_state*)PyModule_GetState(self);
    if (!state) {
        return NULL;
    }

    // extract arguments
    if (!Py_IS_TYPE(args[0], (PyTypeObject*)state->TwoLineElement)) {
        PyErr_SetString(PyExc_TypeError,
                        "first parameter must be TwoLineElement type");
        return NULL;
    }
    tle_t* tle = (tle_t*)args[0];

    double jd;
    if (jd_AsDouble(args[1], (PyTypeObject*)state->JulianDate,
                    jd, "second") < 0)
    {
        return NULL;
    }

    /* we're not using the Vector buffer here because we would
    still need to instantiate the Vectors, then change the already
    initialized valuse. seems redundant to not just initialize
    with known values and using a tiny more memory in doing so */

    // run the sgp4 library function
    double r[3], v[3];
    double dt = (jd - TLE_JDNUMBER(tle) - TLE_JDFRACTION(tle)) * 1440.0;
    SGP4Funcs::sgp4(tle->satrec, dt, r, v);

    // check if the error flag was set
    if (tle->satrec.error) {
        PyObject* exc = PyErr_NewException("sgp4.sgp4Error",
                                           NULL, NULL);
        PyErr_Format(exc, "error in sgp4 library (error = %i)",
                     tle->satrec.error);
        Py_DECREF(exc);
        return NULL;
    }

    // get buffer from new vectors here 

    PyObject* position = build_vector(state->Vector,
                                      r[0], r[1], r[2]);
    if (!position) {
        return NULL;
    }

    PyObject* velocity = build_vector(state->Vector,
                                      v[0], v[1], v[2]);
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

/* this should be replacable with access to Vector buffer */
/*static int
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
}*/

#define TLE_STR_FORMAT  "%s\n%s\n%s"

static PyObject*
_get_tle_repr(const tle_t* self, const char* format)
{
    const size_t buffer_size = snprintf(NULL, 0, format,
                                        self->name, self->line1, self->line2);

    char* buffer = (char*)malloc(buffer_size + 1);
    if (!buffer) {
        return PyErr_NoMemory();
    }

    sprintf(buffer, format, self->name, self->line1, self->line2);

    PyObject* rtn = PyUnicode_FromString(buffer);
    free(buffer);

    return rtn;
}

static PyObject*
tle_str(const tle_t* self)
{
    return _get_tle_repr(self, TLE_STR_FORMAT);
}

static PyObject*
tle_repr(const tle_t* self)
{
    const char* format = "TwoLineElement('"
                         TLE_STR_FORMAT"')";
    return _get_tle_repr(self, format);
}

static int
check_buffer(Py_buffer* view, const char* argNum)
{
    if (view->ndim != 1) {
        PyErr_Format(PyExc_TypeError,
                     "expected %s argument to export a buffer with "
                     "1-dimension (%i found)",
                     argNum, view->ndim);
        return -1;
    }
    else if (*view->shape != 3) {
        PyErr_Format(PyExc_TypeError,
                     "expected %s argument to export a buffer with "
                     "shape = (3,) ((%i,) found)",
                     argNum, *view->shape);
        return -1;
    }
    else if (strcmp(view->format, "d")) {
        PyErr_Format(PyExc_TypeError,
                     "expected %s argument to export a buffer with "
                     "format type 'd' (%'s' found)",
                     argNum, view->format);
        return -1;
    }

    return 0;
}

static PyObject*
elements_from_state(PyObject* self, PyObject* args)
{
    // signiture: _elements_from_state(position: Vector, velocity: Vector, MU: float = EARTH_MU)
    //                  -> (raan, inc, aop, ecc, a, m, nu)

    double mu = EARTH_MU;
    Py_buffer* positionView = (Py_buffer*)malloc(sizeof(Py_buffer));
    if (!positionView) {
        return PyErr_NoMemory();
    }
    Py_buffer* velocityView = (Py_buffer*)malloc(sizeof(Py_buffer));
    if (!velocityView) {
        return PyErr_NoMemory();
    }

    if (!PyArg_ParseTuple(args, "y*y*|d", positionView, velocityView, &mu)) {
        return NULL;
    }

    if (check_buffer(positionView, "first") < 0) {
        return NULL;
    }
    if (check_buffer(velocityView, "second") < 0) {
        return NULL;
    }

    double p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;
    SGP4Funcs::rv2coe_SGP4((double*)positionView->buf,
                           (double*)velocityView->buf,
                           mu, p, a, ecc, incl, omega,
                           argp, nu, m, arglat, truelon,
                           lonper);

    PyBuffer_Release(positionView);
    PyBuffer_Release(velocityView);

    // returning mean anomaly right now with this 
    return Py_BuildValue("ddddddd", omega, incl, argp, ecc, a, m, nu);
}

static PyObject*
elements_from_tle(PyObject* self, PyObject* const* args, Py_ssize_t size)
{
    // signiture: _elements_from_tle(tle: TwoLineElement, jd: JulianDate,
    //              MU: float = EARTH_MU) -> (raan, inc, aop, ecc, a, m, nu)

    double mu = EARTH_MU;
    if (size == 3) {
        mu = PyFloat_AsDouble(args[2]);
        if (mu == -1.0 && PyErr_Occurred()) {
            return NULL;
        }
    }
    else if (size != 2) {
        PyErr_Format(PyExc_TypeError, 
                     "_sgp4._getState() takes two or three arguments"
                     "(%i given)", size);
        return NULL;
    }

    sgp4_state* state = (sgp4_state*)PyModule_GetState(self);
    if (!state) {
        return NULL;
    }

    if (!Py_IS_TYPE(args[0], (PyTypeObject*)state->TwoLineElement)) {
        PyErr_SetString(PyExc_TypeError,
                        "first argument must be TwoLineElement type");
        return NULL;
    }

    double jd;
    // also checks type
    jd_AsDouble(args[1], (PyTypeObject*)state->JulianDate,
                jd, "second");

    // get state
    tle_t* tle = (tle_t*)args[0];
    double r[3], v[3];
    double dt = (jd - TLE_JDNUMBER(tle) - TLE_JDFRACTION(tle)) * 1440;
    SGP4Funcs::sgp4(tle->satrec, dt, r, v);

    // get elements
    double p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;
    SGP4Funcs::rv2coe_SGP4(r, v, mu, p, a, ecc, incl, omega,
                           argp, nu, m, arglat, truelon, lonper);

    // returning mean anomaly right now with this 
    return Py_BuildValue("ddddddd", omega, incl, argp, ecc, a, m, nu);
}

static PyObject* 
compute_ecc_vector(PyObject* self, PyObject* args)
{
    // signiture: compute_ecc_vector(position: Vector, velocity: Vector,
    //                               mu: float = EARTH_MU) -> Vector

    double mu = EARTH_MU;
    Py_buffer* positionView = (Py_buffer*)malloc(sizeof(Py_buffer));
    if (!positionView) {
        return PyErr_NoMemory();
    }
    Py_buffer* velocityView = (Py_buffer*)malloc(sizeof(Py_buffer));
    if (!velocityView) {
        return PyErr_NoMemory();
    }

    if (!PyArg_ParseTuple(args, "y*y*|d", positionView, 
                          velocityView, &mu))
    {
        return NULL;
    }

    if (check_buffer(positionView, "first") < 0) {
        return NULL;
    }
    if (check_buffer(velocityView, "second") < 0) {
        return NULL;
    }

    // get module state
    sgp4_state* state = (sgp4_state*)PyModule_GetState(self);
    if (!state) {
        return NULL;
    }

    double ebar[3];
    double* r = (double*)positionView->buf;
    double* v = (double*)velocityView->buf;
    double magr = SGP4Funcs::mag_SGP4(r);
    double magv = SGP4Funcs::mag_SGP4(v);
    double c1 = magv * magv - mu / magr;
    double rdotv = SGP4Funcs::dot_SGP4(r, v);
    for (int i = 0; i < 3; i++) {
        ebar[i] = (c1 * r[i] - rdotv * v[i]) / mu;
    }

    PyBuffer_Release(positionView);
    PyBuffer_Release(velocityView);

    // build and return eccentric vector
    PyObject* eccVec = build_vector(state->Vector, ebar[0],
                                    ebar[1], ebar[2]);
    return eccVec;
}

static PyMethodDef sgp4_methods[] = {
    {"getState", (PyCFunction)get_state, METH_FASTCALL,
     PyDoc_STR("_getState(tle, jd) -> (position, velocity)\n\nReturns a "
     "tuple of state vectors at the given time.")},

    {"elementsFromState", (PyCFunction)elements_from_state, METH_VARARGS,
     PyDoc_STR("elementsFromState(position, velocity, mu) -> (sma, ecc, "
     "inc, raan, aop, m, nu)\n\nReturns a tuple of orbital elements from a "
     "known state. Tuple returned contains the values (raan, inc, aop, ecc, "
     "sma, meanAnomaly, trueAnomaly).")},

    {"elementsFromTle", (PyCFunction)elements_from_tle, METH_FASTCALL,
     PyDoc_STR("elementsFromTle(tle, jd, mu) -> (sma, ecc, inc, raan, "
     "aop, m, nu)\n\nReturns a tuple of elements from a TLE at a given time "
     "(state is implicitly computed). If you have a known state, use "
     "sgp4.elementsFromState() instead. Tuple returned contains the values "
     "(raan, inc, aop, ecc, sma, meanAnomaly, trueAnomaly).")},

    {"computeEccentricVector", (PyCFunction)compute_ecc_vector,
     METH_VARARGS, PyDoc_STR("computeEccentricVector(position, velocity, mu) "
     "-> eccentricVector\n\nReturns the eccentric vector of a satellite from "
     "a known state.")},

    {NULL, NULL}
};

static PyObject*
tle_epoch(tle_t* self, void* Py_UNUSED(_))
{
    sgp4_state* state = (sgp4_state*)PyType_GetModuleState(Py_TYPE(self));
    if (!state) {
        return NULL;
    }

    double jdValue = TLE_JDNUMBER(self) + TLE_JDFRACTION(self);
    return PyObject_CallMethod(state->JulianDate, "fromNumber",
                               "d", jdValue);
}

static PyObject*
tle_sma(tle_t* self, void* Py_UNUSED(_))
{
    return PyFloat_FromDouble(self->satrec.a * 6371.0);
}

static PyObject*
tle_periapsis(tle_t* self, void* Py_UNUSED(_))
{
    return PyFloat_FromDouble((self->satrec.altp + 1.0) * 6371.0);
}

static PyObject*
tle_apoapsis(tle_t* self, void* Py_UNUSED(_))
{
    return PyFloat_FromDouble((self->satrec.alta + 1.0) * 6371.0);
}

static PyObject*
tle_satnum(tle_t* self, void* Py_UNUSED(_))
{
    return PyUnicode_FromString(self->satrec.satnum);
}

static PyGetSetDef tle_getset[] = {
    {"epoch", (getter)tle_epoch, NULL, 
     PyDoc_STR("The epoch of the TLE as a JulianDate type."), NULL},

    {"sma", (getter)tle_sma, NULL,
     PyDoc_STR("The semi-major axis of the orbit, in kilometers."), NULL},

    {"apogee", (getter)tle_apoapsis, NULL,
     PyDoc_STR("The apogee of the orbit, converted from Earth radii to "
               "kilometers using Earth's equitorial radius."), NULL},

    {"perigee", (getter)tle_periapsis, NULL,
     PyDoc_STR("The perigee of the orbit, converted from Earth radii to "
               "kilometers using Earth's equitorial radius."), NULL},

    {"satNumber", (getter)tle_satnum, NULL,
     PyDoc_STR("The satellite catalog number."), NULL},

    {NULL}
};

static PyMemberDef tle_member[] = {
    {"name", T_STRING, offsetof(tle_t, name), READONLY,
     PyDoc_STR("Name of the satllite if provided when creating instance.")},

    {"line1", T_STRING, offsetof(tle_t, line1), READONLY,
     PyDoc_STR("First line of the TLE.")},

    {"line2", T_STRING, offsetof(tle_t, line2), READONLY,
     PyDoc_STR("Second line of the TLE.")},

    {"error", T_INT, offsetof(tle_t, satrec.error), READONLY,
     PyDoc_STR("Error code from SGP4 library methods.")},

    {"nddot", T_DOUBLE, offsetof(tle_t, satrec.nddot), READONLY,
     PyDoc_STR("Second derivative of mean motion (n double dot).")},

    {"ndot", T_DOUBLE, offsetof(tle_t, satrec.ndot), READONLY,
     PyDoc_STR("Derivative of mean motion (n dot).")},

    {"bstar", T_DOUBLE, offsetof(tle_t, satrec.bstar), READONLY,
     PyDoc_STR("BStar value of the satellite.")},

    {"inc", T_DOUBLE, offsetof(tle_t, satrec.inclo), READONLY,
     PyDoc_STR("Inclination term of the TLE.")},

    {"raan", T_DOUBLE, offsetof(tle_t, satrec.nodeo), READONLY,
     PyDoc_STR("Right-ascension term of the TLE.")},

    {"ecc", T_DOUBLE, offsetof(tle_t, satrec.ecco), READONLY,
     PyDoc_STR("Eccentricity term of the TLE.")},

    {"aop", T_DOUBLE, offsetof(tle_t, satrec.argpo), READONLY,
     PyDoc_STR("Argument of perigee of the TLE.")},

    {"meanAnomaly", T_DOUBLE, offsetof(tle_t, satrec.mo), READONLY,
     PyDoc_STR("Mean anomaly term of the TLE.")},

    {"meanMotion", T_DOUBLE, offsetof(tle_t, satrec.no_kozai), READONLY,
     PyDoc_STR("Kozai'd mean motion term of the TLE.")},

    {"meanMotionUnkozai", T_DOUBLE, offsetof(tle_t, satrec.no_unkozai),
     READONLY, PyDoc_STR("Unkozai'd mean motion term of the TLE.")},

    {"classification", T_CHAR, offsetof(tle_t, satrec.classification),
     READONLY, PyDoc_STR("Classification of the satellite. U=unclassified, "
                         "C=classified, S=secret.")},

    {"ephemeris", T_INT, offsetof(tle_t, satrec.ephtype), READONLY,
     PyDoc_STR("Ephemeris type of the TLE.")},

    {"setNumber", T_LONG, offsetof(tle_t, satrec.elnum), READONLY,
     PyDoc_STR("Element set number of the TLE.")},

    {"revNumber", T_LONG, offsetof(tle_t, satrec.revnum), READONLY,
     PyDoc_STR("Revolution number of the satellite.")},

    // testing what these are
    {"dia_mm", T_LONG, offsetof(tle_t, satrec.dia_mm), READONLY,
     PyDoc_STR("RSO dia in mm.")},

    {"period", T_DOUBLE, offsetof(tle_t, satrec.period_sec), READONLY,
     PyDoc_STR("Period in seconds.")},

    {"active", T_UBYTE, offsetof(tle_t, satrec.active), READONLY,
     PyDoc_STR("Active S/C flag (0=n, 1=y")},

    {"not_orbital", T_UBYTE, offsetof(tle_t, satrec.not_orbital), READONLY,
     PyDoc_STR("Orbiting S/C flag (0=n, 1=y)")},

    {"rcs_m2", T_DOUBLE, offsetof(tle_t, satrec.rcs_m2), READONLY,
     PyDoc_STR("'RCS (m^2)' storage")},

    {NULL},
};

static PyType_Slot tle_slots[] = {
    {Py_tp_getset,  tle_getset},
    {Py_tp_init,    (initproc)tle_init},
    {Py_tp_new,     PyType_GenericNew},
    {Py_tp_members, tle_member},
    {Py_tp_free,    (freefunc)tle_free},
    {Py_tp_str,     (reprfunc)tle_str},
    {Py_tp_repr,    (reprfunc)tle_repr},
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
    if (!state) {
        return -1;
    }

    PyObject* type = PyType_FromModuleAndSpec(module, &tle_spec, NULL);
    if (!type) {
        return -1;
    }
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
    type = PyObject_GetAttrString(module_import, "Vector");
    Py_DECREF(module_import);
    if (!type) {
        Py_DECREF(state->TwoLineElement);
        return -1;
    }
    state->Vector = type;
    type = NULL;

    module_import = PyImport_ImportModule("sattrack.spacetime.juliandate");
    if (!module_import) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->Vector);
        return -1;
    }
    type = PyObject_GetAttrString(module_import, "JulianDate");
    Py_DECREF(module_import);
    if (!type) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->Vector);
        return -1;
    }
    state->JulianDate = type;
    type = NULL;

    if (PyModule_AddIntConstant(module, "WGS72OLD", (int)gravconsttype::wgs72old) < 0) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->Vector);
        Py_DECREF(state->JulianDate);
        return -1;
    }
    if (PyModule_AddIntConstant(module, "WGS72", (int)gravconsttype::wgs72) < 0) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->Vector);
        Py_DECREF(state->JulianDate);
        return -1;
    }
    if (PyModule_AddIntConstant(module, "WGS84", (int)gravconsttype::wgs84) < 0) {
        Py_DECREF(state->TwoLineElement);
        Py_DECREF(state->Vector);
        Py_DECREF(state->JulianDate);
        return -1;
    }

    return 0;
}

static void sgp4_free(void* self) 
{
    sgp4_state* state = (sgp4_state*)PyModule_GetState((PyObject*)self);
    if (!state) {
        return;
    }

    Py_DECREF(state->Vector);
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
    "sgp4",		            /* m_name */
    tle_doc,	            /* m_doc */
    sizeof(sgp4_state),	    /* m_size */
    sgp4_methods,           /* m_methods */
    sgp4_slots,             /* m_slots */
    NULL,                   /* m_traverse */
    NULL,                   /* m_clear */
    sgp4_free,              /* m_free */
};

PyMODINIT_FUNC
PyInit_sgp4(void)
{
    return PyModuleDef_Init(&tle_module);
}
