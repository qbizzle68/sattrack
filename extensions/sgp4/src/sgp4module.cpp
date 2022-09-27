#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "SGP4.h"
#define pi 3.14159265358979323846

typedef struct
{
	PyObject_HEAD
	elsetrec satrec;
} _SGP4_Propagator;

static void _elsetrec_init(_SGP4_Propagator* prop, PyObject* tle) {
	const double deg2rad = pi / 180.0;         //   0.0174532925199433
	const double xpdotp = 1440.0 / (2.0 * pi);  // 229.1831180523293
	PyObject* next = NULL;

	PyObject* tle_epoch = PyObject_GetAttrString(tle, "_epoch");
	next = PyObject_GetAttrString(tle_epoch, "_day_number");
	prop->satrec.jdsatepoch = PyLong_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle_epoch, "_day_fraction");
	prop->satrec.jdsatepochF = PyFloat_AsDouble(next);
	Py_DECREF(next);
	Py_DECREF(tle_epoch);
	tle_epoch = NULL;
	next = PyObject_GetAttrString(tle, "_meanMotion");
	prop->satrec.no_kozai = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_ecc");
	prop->satrec.ecco = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_inc");
	prop->satrec.inclo = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_raan");
	prop->satrec.nodeo = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_aop");
	prop->satrec.argpo = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_meanAnom");
	prop->satrec.mo = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_nDDot");
	prop->satrec.nddot = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_bStar");
	prop->satrec.bstar = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_nDot");
	prop->satrec.ndot = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_setNumber");
	prop->satrec.elnum = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = PyObject_GetAttrString(tle, "_revNum");
	prop->satrec.revnum = PyFloat_AsDouble(next);
	Py_DECREF(next);
	next = NULL;
	prop->satrec.classification = 'U'; // this shouldn't matter
	// not setting satrec.intldesg, shouldnt matter
	prop->satrec.ephtype = 0; // hard-coded

	// convert units and initialize
	prop->satrec.no_kozai = prop->satrec.no_kozai / xpdotp; //* rad/min
	prop->satrec.ndot = prop->satrec.ndot / (xpdotp * 1440.0);  //* ? * minperday
	prop->satrec.nddot = prop->satrec.nddot / (xpdotp * 1440.0 * 1440);
	prop->satrec.inclo = prop->satrec.inclo * deg2rad;
	prop->satrec.nodeo = prop->satrec.nodeo * deg2rad;
	prop->satrec.argpo = prop->satrec.argpo * deg2rad;
	prop->satrec.mo = prop->satrec.mo * deg2rad;

	SGP4Funcs::sgp4init(wgs72, 'i', prop->satrec.satnum, prop->satrec.jdsatepoch + prop->satrec.jdsatepochF - 2433281.5, prop->satrec.bstar,
		prop->satrec.ndot, prop->satrec.nddot, prop->satrec.ecco, prop->satrec.argpo, prop->satrec.inclo, prop->satrec.mo, prop->satrec.no_kozai,
		prop->satrec.nodeo, prop->satrec);
}

//static PyObject* update(_SGP4_Propagator prop)

static PyObject* getState(_SGP4_Propagator* prop, PyObject* const* args, Py_ssize_t size)
{
    if (size != 2) {
		PyErr_SetString(PyExc_TypeError, "getState takes exactly two arguments");
		return NULL;
	}
	PyObject* tle = args[0];
	PyObject* jd = args[1];

	if (!tle && !jd) {
		PyErr_SetString(PyExc_ValueError, "Arguments cannon be NULL.");
		return NULL;
	}
	if (static_cast<std::string>(tle->ob_type->tp_name) != "TwoLineElement") {
		PyErr_SetString(PyExc_TypeError, "First argument must be TwoLineElement type.");
		return NULL;
	}
	if (static_cast<std::string>(jd->ob_type->tp_name) != "JulianDate") {
		PyErr_SetString(PyExc_TypeError, "Second argument must be JulianDate type.");
		return NULL;
	}

	PyObject* tle_epoch = PyObject_GetAttrString(tle, "_epoch");
	PyObject* tle_jd_num = PyObject_GetAttrString(tle_epoch, "_day_number");
	PyObject* tle_jd_frac = PyObject_GetAttrString(tle_epoch, "_day_fraction");
	PyObject* jd_num = PyObject_GetAttrString(jd, "_day_number");
	PyObject* jd_frac = PyObject_GetAttrString(jd, "_day_fraction");

	double tsince = (PyFloat_AsDouble(jd_num)
		- PyFloat_AsDouble(tle_jd_num)
		+ PyFloat_AsDouble(jd_frac)
		- PyFloat_AsDouble(tle_jd_frac)) * 1440.0;
	Py_DECREF(tle_epoch);
	Py_DECREF(tle_jd_num);
	Py_DECREF(tle_jd_frac);
	Py_DECREF(jd_num);
	Py_DECREF(jd_frac);
	double r[3], v[3];
	SGP4Funcs::sgp4(prop->satrec, tsince, r, v);

	if (prop->satrec.error > 0) {
		PyErr_SetString(PyExc_Exception, "Error occured in sgp4().");
		return NULL;
	}

	PyObject* rtuple = Py_BuildValue("(ddd)", r[0], r[1], r[2]);
	if (!rtuple) {
		PyErr_SetString(PyExc_ValueError, "Error building position tuple.");
		return NULL;
	}
	PyObject* vtuple = Py_BuildValue("(ddd)", v[0], v[1], v[2]);
	if (!rtuple) {
		PyErr_SetString(PyExc_ValueError, "Error building position tuple.");
		Py_DECREF(rtuple);
		return NULL;
	}

	PyObject* rtn = PyTuple_Pack(2, rtuple, vtuple);
	if (!rtn) {
		PyErr_SetString(PyExc_ValueError, "Error building return tuple.");
		Py_DECREF(rtuple);
		Py_DECREF(vtuple);
		return NULL;
	}

	return rtn;
}

static int _sgp4_init(_SGP4_Propagator* self, PyObject* args, PyObject* Py_UNUSED) {
	PyObject* tle = NULL;

	if (!PyArg_ParseTuple(args, "O", &tle))
		return -1;

	if (static_cast<std::string>(tle->ob_type->tp_name) != "TwoLineElement") {
		PyErr_SetString(PyExc_TypeError, "Argument must be TwoLineElement type.");
		Py_DECREF(tle);
		return -1;
	}

	_elsetrec_init(self, tle);
	if (self->satrec.error > 0) {
		PyErr_SetString(PyExc_Exception, "Error in sgp4init().");
		Py_DECREF(tle);
		return -1;
	}

	return 0;
}

static PyMethodDef sgp4_methods[] = {
	{"getState", (PyCFunction)getState, METH_FASTCALL, "Returns the state of the satellite at a given time."},
	{NULL}
};

static PyTypeObject _sgp4_propType = {
	PyVarObject_HEAD_INIT(NULL, 0)
};

static PyModuleDef sgp4module = {
	PyModuleDef_HEAD_INIT,
	"sgp4",
	"sgp4 base methods",
	-1,
};

PyMODINIT_FUNC
PyInit_sgp4(void)
{
	_sgp4_propType.tp_name		= "sgp4.SGP4_Propagator";
	_sgp4_propType.tp_basicsize	= sizeof(_SGP4_Propagator);
	_sgp4_propType.tp_itemsize	= 0;
	_sgp4_propType.tp_flags		= Py_TPFLAGS_DEFAULT;
	_sgp4_propType.tp_doc		= "Base class which support SGP4 propagation.";
	_sgp4_propType.tp_methods	= sgp4_methods;
	_sgp4_propType.tp_init		= (initproc)_sgp4_init;
	_sgp4_propType.tp_new		= PyType_GenericNew;

	PyObject* m = NULL;

	if (PyType_Ready(&_sgp4_propType) < 0)
		return NULL;

	m = PyModule_Create(&sgp4module);
	if (!m)
		return NULL;

	if (PyModule_AddObject(m, "SGP4_Propagator", (PyObject*)&_sgp4_propType) < 0) {
		Py_DECREF(m);
		return NULL;
	}

	Py_INCREF(&_sgp4_propType);
	return m;
		
}
