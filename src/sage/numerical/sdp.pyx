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

    :meth:`~SemidefiniteProgram.add_constraint`            | Adds a constraint to the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.base_ring`                 | Return the base ring
    :meth:`~SemidefiniteProgram.constraints`               | Returns a list of constraints, as 3-tuples
    :meth:`~SemidefiniteProgram.get_backend`               | Returns the backend instance used
    :meth:`~SemidefiniteProgram.get_max`                   | Returns the maximum value of a variable
    :meth:`~SemidefiniteProgram.get_min`                   | Returns the minimum value of a variable
    :meth:`~SemidefiniteProgram.get_values`                | Return values found by the previous call to ``solve()``
    :meth:`~SemidefiniteProgram.linear_constraints_parent` | Return the parent for all linear constraints
    :meth:`~SemidefiniteProgram.linear_function`           | Construct a new linear function
    :meth:`~SemidefiniteProgram.linear_functions_parent`   | Return the parent for all linear functions
    :meth:`~SemidefiniteProgram.new_variable`              | Returns an instance of ``SDPVariable`` associated
    :meth:`~SemidefiniteProgram.remove_constraint`         | Removes a constraint from self
    :meth:`~SemidefiniteProgram.remove_constraints`        | Remove several constraints
    :meth:`~SemidefiniteProgram.set_max`                   | Sets the maximum value of a variable
    :meth:`~SemidefiniteProgram.set_min`                   | Sets the minimum value of a variable
    :meth:`~SemidefiniteProgram.set_objective`             | Sets the objective of the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.set_problem_name`          | Sets the name of the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.show`                      | Displays the ``SemidefiniteProgram`` in a human-readable
    :meth:`~SemidefiniteProgram.solve`                     | Solves the ``SemidefiniteProgram``

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
from sage.numerical.linear_functions import is_LinearFunction, is_LinearConstraint
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

        self._backend.set_objective(values, d)

    def add_constraint(self, linear_function, max=None, min=None, name=None):
        if linear_function is 0:
            return

        # Raising an exception when min/max are not as expected
        from sage.rings.all import RR
        if ((min is not None and min not in RR)
            or (max is not None and max not in RR)):
            raise ValueError("min and max arguments are required to be numerical")

        if is_LinearFunction(linear_function):
            f = linear_function.dict()
            constant_coefficient = f.get(-1,0)

            # We do not want to ignore the constant coefficient
            max = (max - constant_coefficient) if max is not None else None
            min = (min - constant_coefficient) if min is not None else None

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

            self._backend.add_linear_constraint(C, min, max, name)

        elif is_LinearConstraint(linear_function):
            constraint = linear_function
            for lhs, rhs in constraint.equations():
                self.add_constraint(lhs-rhs, min=0, max=0, name=name)
            for lhs, rhs in constraint.inequalities():
                self.add_constraint(lhs-rhs, max=0, name=name)
        else:
            raise ValueError('argument must be a linear function or constraint, got '+str(linear_function))

    def solve(self, log=None, objective_only=False):
        if log is not None: self._backend.set_verbosity(log)
        self._backend.solve()
        return self._backend.get_objective_value()

    def set_min(self, v, min):
        self._backend.variable_lower_bound(self._variables[v], min)

    def set_max(self, v, max):
        self._backend.variable_upper_bound(self._variables[v], max)

    def get_min(self, v):
        return self._backend.variable_lower_bound(self._variables[v])

    def get_max(self, v):
        return self._backend.variable_upper_bound(self._variables[v])

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
        return self.linear_functions_parent()(d)

    def get_backend(self):
        return self._backend

class SDPSolverException(RuntimeError):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

cdef class SDPVariable(SageObject):

    def __cinit__(self, p, vtype, dim=1, name=""):
        self._dim = dim
        self._dict = {}
        self._p = p
        self._vtype = vtype

        self._hasname = (len(name) >0)

        if dim > 1:
            from sage.misc.superseded import deprecation
            deprecation(15489, "The 'dim' argument will soon disappear. "+
                        "Fortunately variable[1,2] is easier to use than "+
                        "variable[1][2]")

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
        elif self._dim == 1:
            zero = self._p._backend.zero()
            j = self._p._backend.add_variable(zero , None, False, True, False, zero,
                                              (str(self._name) + "[" + str(i) + "]")
                                               if self._hasname else None)

            v = self._p.linear_function({j : 1})
            self._p._variables[v] = j
            self._p._backend.set_variable_type(j,self._vtype)
            self._dict[i] = v

            return v

        else:
            self._dict[i] = SDPVariable(
                self._p,
                self._vtype,
                dim=self._dim-1,
                name = ("" if not self._hasname
                        else (str(self._name) + "[" + str(i) + "]")))

            return self._dict[i]

    def _repr_(self):
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
