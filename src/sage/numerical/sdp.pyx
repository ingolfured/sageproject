r"""
AUTHORS:

-Ingolfur Edvardsson  (2014/06):

Index of functions and methods
------------------------------

Below are listed the methods of :class:`SemidefiniteProgram`. This module
also implements the :class:`SDPSolverException` exception, as well as the
:class:`SDPVariable` class.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~SemidefiniteProgram.add_constraint`                | Adds a constraint to the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.base_ring`                     | Return the base ring
    :meth:`~SemidefiniteProgram.constraints`                   | Returns a list of constraints, as 3-tuples
    :meth:`~SemidefiniteProgram.get_backend`                   | Returns the backend instance used
    :meth:`~SemidefiniteProgram.get_values`                    | Return values found by the previous call to ``solve()``
    :meth:`~SemidefiniteProgram.linear_sdp_constraints_parent` | Return the parent for all linear constraints
    :meth:`~SemidefiniteProgram.linear_sdp_function`           | Construct a new linear function
    :meth:`~SemidefiniteProgram.linear_sdp_functions_parent`   | Return the parent for all linear functions
    :meth:`~SemidefiniteProgram.new_variable`                  | Returns an instance of ``SDPVariable`` associated
    :meth:`~SemidefiniteProgram.remove_constraint`             | Removes a constraint from self
    :meth:`~SemidefiniteProgram.remove_constraints`            | Remove several constraints
    :meth:`~SemidefiniteProgram.set_objective`                 | Sets the objective of the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.set_problem_name`              | Sets the name of the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.show`                          | Displays the ``SemidefiniteProgram`` in a human-readable
    :meth:`~SemidefiniteProgram.solve`                         | Solves the ``SemidefiniteProgram``

Classes and methods
-------------------
"""

#*****************************************************************************
#       Copyright (C) 2014 Ingolfur Edvardsson <ingolfured@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"

from sage.structure.sage_object cimport SageObject
from sage.misc.cachefunc import cached_method
from sage.numerical.linear_sdp_functions import is_LinearSDPFunction, is_LinearSDPConstraint
from sage.misc.superseded import deprecated_function_alias, deprecation
from sage.misc.superseded import deprecated_function_alias

cdef class SemidefiniteProgram(SageObject):

    def __init__(self, solver=None, maximization=True,
                 constraint_generation=False, check_redundant=False):

        from sage.numerical.backends.generic_sdp_backend import get_solver
        self._backend = get_solver(solver=solver,
                                   constraint_generation=constraint_generation)
        if not maximization:
            self._backend.set_sense(-1)

        # Associates an index to the variables
        self._variables = {}

        # Check for redundant constraints
        self._check_redundant = check_redundant
        if check_redundant:
            self._constraints = list()

    def __call__(self, x):
        parent = self.linear_sdp_functions_parent()
        return parent(x)

    linear_sdp_function = __call__

    def linear_sdp_functions_parent(self):
        if self._linear_sdp_functions_parent is None:
            base_ring = self._backend.base_ring()
            from sage.numerical.linear_sdp_functions import LinearSDPFunctionsParent
            self._linear_sdp_functions_parent = LinearSDPFunctionsParent(base_ring)
        return self._linear_sdp_functions_parent

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        A ring. The coefficients that the chosen solver supports.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='cvxopt')
            sage: p.base_ring()
            Rational Field
        """
        return self._backend.base_ring()

    def linear_sdp_constraints_parent(self):
        """
        Return the parent for all linear constraints

        See :mod:`~sage.numerical.linear_functions` for more
        details.

        EXAMPLES::

             sage: p = MixedIntegerLinearProgram()
             sage: p.linear_sdp_constraints_parent()
             Linear constraints over Real Double Field
        """
        if self._linear_sdp_constraints_parent is None:
            from sage.numerical.linear_sdp_functions import LinearSDPConstraintsParent
            LF = self.linear_sdp_functions_parent()
            self._linear_sdp_constraints_parent = LinearSDPConstraintsParent(LF)
        return self._linear_sdp_constraints_parent

    def get_values(self, *lists):
        val = []
        for l in lists:
            if isinstance(l, SDPVariable):
                if l.depth() == 1:
                    c = {}
                    for (k,v) in l.items():
                        #c[k] = self._values[v] if self._values.has_key(v) else None
                        c[k] = self._backend.get_variable_value(self._variables[v])
                    val.append(c)
                else:
                    c = {}
                    for (k,v) in l.items():
                        c[k] = self.get_values(v)
                    val.append(c)
            elif isinstance(l, list):
                if len(l) == 1:
                    val.append([self.get_values(l[0])])
                else:
                    c = []
                    [c.append(self.get_values(ll)) for ll in l]
                    val.append(c)
            elif l in self._variables:
                #val.append(self._values[l])
                val.append(self._backend.get_variable_value(self._variables[l]))

        if len(lists) == 1:
            return val[0]
        else:
            return val

    def set_objective(self,obj):
        cdef list values = []

        #TODO change this comment!!

        # If the objective is None, or a constant, we want to remember
        # that the objective function has been defined ( the user did not
        # forget it ). In some LP problems, you just want a feasible solution
        # and do not care about any function being optimal.

        cdef int i

        if obj is not None:
            f = obj.dict()
        else:
            f = {-1 : 0}

        d = f.pop(-1,self._backend.zero())

        for i in range(self._backend.ncols()):
            values.append(f.get(i,self._backend.zero()))

        self._backend.set_objective(values)

    def add_constraint(self, linear_sdp_function,  name=None):
        if linear_function is 0:
            return

        if is_LinearFunction(linear_function):
            f = linear_function.dict()
            constant_coefficient = f.get(-1,0)

            indices = []
            values = []

            if self._check_redundant:
              b = self._backend
              from __builtin__ import min as min_function
              i = min_function([v for (v,coeff) in f.iteritems() if coeff != 0])
              c = f[i]
              C = [(v,coeff/c) for (v,coeff) in f.iteritems() if v != -1]
              if c > 0:
                min = min/c if min is not None else None
                max = max/c if max is not None else None
              else:
                tempmin = max/c if max is not None else None
                tempmax = min/c if min is not None else None
                min, max = tempmin, tempmax
              if (tuple(C),min,max) in self._constraints:
                return None
              else:
                self._constraints.append((tuple(C),min,max))
            else:
              C = [(v,coeff) for (v,coeff) in f.iteritems() if v != -1]

            if min is None and max is None:
                raise ValueError("Both max and min are set to None ? Weird!")

            self._backend.add_linear_constraint(C, name)

        elif is_LinearConstraint(linear_function):
            constraint = linear_function
            for lhs, rhs in constraint.equations():
                self.add_constraint(lhs-rhs, min=0, max=0, name=name)
            for lhs, rhs in constraint.inequalities():
                self.add_constraint(lhs-rhs, max=0, name=name)
        else:
            raise ValueError('argument must be a linear function or constraint, got '+str(linear_function))

    def solve(self, objective_only=False):
        self._backend.solve()
        return self._backend.get_objective_value()

    def solver_parameter(self, name, value = None):
        if value is None:
            return self._backend.solver_parameter(name)
        else:
            self._backend.solver_parameter(name, value)

    cpdef sum(self, L):
        d = {}
        for v in L:
            for id,coeff  in v.iteritems():
                d[id] = coeff + d.get(id,0)
        return self.linear_sdp_functions_parent()(d)

    def get_backend(self):
        return self._backend

    def new_variable(self, name=""):
        v=SDPVariable(self, name=name)
        return v

class SDPSolverException(RuntimeError):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

cdef class SDPVariable(SageObject):

    def __cinit__(self, p, name=""):
        self._dict = {}
        self._p = p

        self._hasname = (len(name) >0)


        # create a temporary char *
        cdef char *name_c = name
        # and copy it over
        self._name = <char*>sage_malloc(len(name)+1)
        strcpy(self._name, name_c)

    def __dealloc__(self):
        if self._name:
            sage_free(self._name)

    def __getitem__(self, i):
        cdef SDPVariable s = self

        cdef int j

        if i in self._dict:
            return self._dict[i]
        else:
            zero = self._p._backend.zero()
            j = self._p._backend.add_variable(zero , None,  zero,
                                              (str(self._name) + "[" + str(i) + "]")
                                               if self._hasname else None)

            v = self._p.linear_sdp_function({j : 1 })
            self._p._variables[v] = j
            self._dict[i] = v

            return v

        return "SDPVariable of dimension " + str(self._dim) + "."

    def keys(self):
        return self._dict.keys()

    def items(self):
        return self._dict.items()

    def depth(self):
        return self._dim

    def values(self):
        return self._dict.values()

def Sum(x):
    from sage.misc.superseded import deprecation
    deprecation(13646, 'use SemidefiniteProgram.sum() instead')
    if not x:
        return None
    parent = x[0].parent()
    return parent.sum(x)
