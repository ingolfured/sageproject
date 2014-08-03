from sage.structure.parent cimport Parent
from sage.structure.element cimport ModuleElement, RingElement, Element


cdef class LinearSDPFunctionsParent_class(Parent):
    cpdef _element_constructor_(self, x)
    cpdef _coerce_map_from_(self, R)
    cdef object _multiplication_symbol
    cpdef object _get_multiplication_symbol(self)

cdef class LinearSDPFunction(ModuleElement):
    cdef dict _f
    cdef int _dim
    cpdef iteritems(self)
    cpdef ModuleElement _add_(self, ModuleElement b)
    cpdef ModuleElement _sub_(self, ModuleElement b)
    cpdef ModuleElement _lmul_(self, RingElement b)
    cpdef ModuleElement _rmul_(self, RingElement b)
    cpdef ModuleElement _neg_(self)
    cpdef _acted_upon_(self, x, bint self_on_left)
    cdef _richcmp(left, right, int op)
    cpdef is_zero(self)
    cpdef equals(LinearSDPFunction left, LinearSDPFunction right)
    cdef int _cmp_c_impl(left, Element right) except -2

cdef class LinearSDPConstraintsParent_class(Parent):
    cdef LinearSDPFunctionsParent_class _LF
    cpdef _element_constructor_(self, left, right=?, equality=?)
    cpdef _coerce_map_from_(self, R)

cdef class LinearSDPConstraint(Element):
    cdef bint equality
    cdef list constraints
    cdef _richcmp(left, right, int op)
    cdef LinearSDPConstraint _chained_comparator_hack_part1(LinearSDPConstraint left, LinearSDPConstraint right)
    cdef _chained_comparator_hack_part2(self)
    cpdef equals(LinearSDPConstraint left, LinearSDPConstraint right)

cdef LinearSDPConstraint _chained_comparator_hack_search
cdef LinearSDPConstraint _chained_comparator_hack_replace
