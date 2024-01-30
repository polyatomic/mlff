#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "chemutil.h"
#include "util.h"

typedef struct {
   PyObject_HEAD
   int size;
   Iterator<Mol> mols;
} MolIteratorObject;

static PyObject* MolIterator_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
   MolIteratorObject *self;
   self = (MolIteratorObject*)type->tp_alloc(type, 0);
   if (self != NULL) {
      self->size = 0;
   }
   return (PyObject*)self;
}

static PyObject* MolIterator_Size(MolIteratorObject *self, PyObject *Py_UNUSED(ignored)) {
   return PyLong_FromLong(self->size);
}

static PyObject* MolIterator_ReadXYZFile(MolIteratorObject *self, PyObject *args) {
   char *str;
   if (!PyArg_ParseTuple(args, "s", &str)) {
      return NULL;
   }
   self->mols = ReadXYZFile(str);
   self->size = self->mols.Size();
   if (self->size == 0) {
      PyErr_Format(PyExc_Exception, "Error opening %s", str);
      return NULL;
   }
   Py_RETURN_NONE;
}

static PyMethodDef MolIterator_methods[] = {
   { "Size", (PyCFunction)MolIterator_Size, METH_NOARGS, "Return the number of geometries" },
   { "ReadXYZFile", (PyCFunction)MolIterator_ReadXYZFile, METH_VARARGS, "Read geometry file" },
   { NULL }
};

static PyTypeObject MolIteratorType = {
   PyVarObject_HEAD_INIT(NULL, 0)            // ob_base
   "chemutil.MolIterator",                   // tp_name
   sizeof(MolIteratorObject),                // tp_basicsize
   0,                                        // tp_itemsize
   nullptr,                                  // tp_dealloc
   0,                                        // tp_vectorcall_offset
   nullptr,                                  // tp_getattr
   nullptr,                                  // tp_setattr
   nullptr,                                  // tp_as_async
   nullptr,                                  // tp_repr
   nullptr,                                  // tp_as_number
   nullptr,                                  // tp_as_sequence
   nullptr,                                  // tp_as_mapping
   nullptr,                                  // tp_hash
   nullptr,                                  // tp_call
   nullptr,                                  // tp_str
   nullptr,                                  // tp_getattro
   nullptr,                                  // tp_setattro
   nullptr,                                  // tp_as_buffer
   Py_TPFLAGS_DEFAULT,                       // tp_flags
   PyDoc_STR("MolIterator objects"),         // tp_doc
   nullptr,                                  // tp_traverse
   nullptr,                                  // tp_clear
   nullptr,                                  // tp_richcompare
   0,                                        // tp_weaklistoffset
   nullptr,                                  // tp_iter
   nullptr,                                  // tp_iternext
   MolIterator_methods,                      // tp_methods
   nullptr,                                  // tp_members
   nullptr,                                  // tp_getset
   nullptr,                                  // tp_base
   nullptr,                                  // tp_dict
   nullptr,                                  // tp_descr_get
   nullptr,                                  // tp_descr_set
   0,                                        // tp_dictoffset
   nullptr,                                  // tp_init
   nullptr,                                  // tp_alloc
   MolIterator_new,                          // tp_new
   nullptr,                                  // tp_free
   nullptr,                                  // tp_is_gc
   nullptr,                                  // tp_bases
   nullptr,                                  // tp_mro
   nullptr,                                  // tp_cache
   nullptr,                                  // tp_subclasses
   nullptr,                                  // tp_weaklist
   nullptr,                                  // tp_del
   0,                                        // tp_version_tag
   nullptr,                                  // tp_finalize
   nullptr                                   // tp_vectorcall
};

static PyModuleDef chemutilmodule = {
   PyModuleDef_HEAD_INIT,      // m_base
   "chemutil",                 // m_name
   "Chemutil wrapper module.", // m_doc
   -1,                         // m_size
   nullptr,                    // m_methods
   nullptr,                    // m_slots
   nullptr,                    // m_traverse
   nullptr,                    // m_clear
   nullptr                     // m_free
};

PyMODINIT_FUNC PyInit_chemutil() {
   PyObject *m;
   if (PyType_Ready(&MolIteratorType) < 0)
      return NULL;

   m = PyModule_Create(&chemutilmodule);
   if (m == NULL)
      return NULL;

   Py_INCREF(&MolIteratorType);
   if (PyModule_AddObject(m, "MolIterator", (PyObject*)&MolIteratorType) < 0) {
      Py_DECREF(&MolIteratorType);
      Py_DECREF(m);
      return NULL;
   }

   return m;
}
