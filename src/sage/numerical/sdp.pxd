cdef extern from *:
    ctypedef double* const_double_ptr "const double*"
    cdef int BINARY = 1
    cdef int REAL = -1
    cdef int INTEGER = 0

from sage.structure.sage_object cimport SageObject
from sage.numerical.backends.generic_sdp_backend cimport GenericSDPBackend
cdef class SDPVariable


cdef class SemidefiniteProgram(SageObject):
    cdef GenericSDPBackend _backend
    cdef list _sdpvariables
    cdef SDPVariable _default_sdpvariable
    cdef dict _variables
    cdef object _linear_functions_parent
    cdef object _linear_constraints_parent
    cdef int _check_redundant
    cdef list _constraints
    cpdef sum(self, L)

cdef class SDPVariable(SageObject):
    cdef SemidefiniteProgram _p
    cdef dict _dict
    cdef char * _name
    cdef bint _hasname

