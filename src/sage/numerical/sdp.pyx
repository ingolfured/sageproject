r"""
SemiDefinite Programming

A linear program (`LP <http://en.wikipedia.org/wiki/Linear_programming>`_)
is an `optimization problem <http://en.wikipedia.org/wiki/Optimization_%28mathematics%29>`_
in the following form

.. MATH::
     \max \{ c^T x \;|\; A x \leq b, x \geq 0 \}

with given `A \in \mathbb{R}^{m,n}`, `b \in \mathbb{R}^m`,
`c \in \mathbb{R}^n` and unknown `x \in \mathbb{R}^{n}`.
If some or all variables in the vector `x` are restricted over
the integers `\mathbb{Z}`, the problem is called mixed integer
linear program (`MILP <http://en.wikipedia.org/wiki/Mixed_integer_linear_programming>`_).
A wide variety of problems in optimization
can be formulated in this standard form. Then, solvers are
able to calculate a solution.

Imagine you want to solve the following linear system of three equations:

 - `w_0 + w_1 + w_2 - 14 w_3 = 0`
 - `w_1 + 2 w_2 - 8 w_3 = 0`
 - `2 w_2 - 3 w_3 = 0`

and this additional inequality:

 - `w_0 - w_1 - w_2 \geq 0`

where all `w_i \in \mathbb{Z}^+`. You know that the trivial solution is `w_i=0`,
but what is the first non-trivial one with `w_3 \geq 1`?

A mixed integer linear program can give you an answer:

  #. You have to create an instance of :class:`SemidefiniteProgram` and
     -- in our case -- specify that it is a minimization.
  #. Create an dictionary ``w`` of integer variables ``w`` via ``w =
     p.new_variable(integer=True)`` (note that **by default all variables are
     non-negative**, cf :meth:`~SemidefiniteProgram.new_variable`).
  #. Add those three equations as equality constraints via
     :meth:`add_constraint <sage.numerical.sdp.SemidefiniteProgram.add_constraint>`.
  #. Also add the inequality constraint.
  #. Add an inequality constraint `w_3 \geq 1` to exclude the trivial solution.
  #. By default, all variables are non-negative. We remove that constraint
     via ``p.set_min(variable, None)``, see :meth:`set_min <sage.numerical.sdp.SemidefiniteProgram.set_min>`.
  #. Specify the objective function via :meth:`set_objective <sage.numerical.sdp.SemidefiniteProgram.set_objective>`.
     In our case that is just `w_3`. If it
     is a pure constraint satisfaction problem, specify it as ``None``.
  #. To check if everything is set up correctly, you can print the problem via
     :meth:`show <sage.numerical.sdp.SemidefiniteProgram.show>`.
  #. :meth:`Solve <sage.numerical.sdp.SemidefiniteProgram.solve>` it and print the solution.

The following example shows all these steps::

    sage: p = SemidefiniteProgram(maximization=False, solver = "cvxopt")
    sage: w = p.new_variable(integer=True, nonnegative=True)
    sage: p.add_constraint(w[0] + w[1] + w[2] - 14*w[3] == 0)
    sage: p.add_constraint(w[1] + 2*w[2] - 8*w[3] == 0)
    sage: p.add_constraint(2*w[2] - 3*w[3] == 0)
    sage: p.add_constraint(w[0] - w[1] - w[2] >= 0)
    sage: p.add_constraint(w[3] >= 1)
    sage: _ = [ p.set_min(w[i], None) for i in range(1,4) ]
    sage: p.set_objective(w[3])
    sage: p.show()
    Minimization:
       x_3
    Constraints:
      0.0 <= x_0 + x_1 + x_2 - 14.0 x_3 <= 0.0
      0.0 <= x_1 + 2.0 x_2 - 8.0 x_3 <= 0.0
      0.0 <= 2.0 x_2 - 3.0 x_3 <= 0.0
      - x_0 + x_1 + x_2 <= 0.0
      - x_3 <= -1.0
    Variables:
      x_0 is an integer variable (min=0.0, max=+oo)
      x_1 is an integer variable (min=-oo, max=+oo)
      x_2 is an integer variable (min=-oo, max=+oo)
      x_3 is an integer variable (min=-oo, max=+oo)
    sage: print 'Objective Value:', p.solve()
    Objective Value: 2.0
    sage: for i, v in p.get_values(w).iteritems():\
              print 'w_%s = %s' % (i, int(round(v)))
    w_0 = 15
    w_1 = 10
    w_2 = 3
    w_3 = 2

Different backends compute with different base fields, for example::

    sage: p = SemidefiniteProgram(solver='cvxopt')
    sage: p.base_ring()
    Real Double Field
    sage: x = p.new_variable(real=True, nonnegative=True)
    sage: 0.5 + 3/2*x[1]
    0.5 + 1.5*x_0

    sage: p = SemidefiniteProgram(solver='ppl')
    sage: p.base_ring()
    Rational Field
    sage: x = p.new_variable(nonnegative=True)
    sage: 0.5 + 3/2*x[1]
    1/2 + 3/2*x_0


Linear Variables and Expressions
--------------------------------

The underlying linear programming backends always work with matrices
where each column corresponds to a linear variable. These variables
can be accessed using the :meth:`SemidefiniteProgram.gen` method
or by calling with a dictionary variable index to coefficient::

    sage: sdp = SemidefiniteProgram()
    sage: 5 + sdp.gen(0) + 2*sdp.gen(1)
    5 + x_0 + 2*x_1
    sage: sdp({-1:5, 0:1, 1:2})
    5 + x_0 + 2*x_1

However, this alone is often not convenient to construct a linear
program. To make your code more readable, you can construct
:class:`SDPVariable` objects that can be arbitrarily named and
indexed. Internally, this is then translated back to the `x_i`
variables. For example::

    sage: sdp.<a,b> = SemidefiniteProgram()
    sage: a
    SDPVariable of dimension 1.
    sage: 5 + a[1] + 2*b[3]
    5 + x_0 + 2*x_1

Indices can be any object, not necessarily integers. Multi-indices are
also allowed::

    sage: a[4, 'string', QQ]
    x_2
    sage: a[4, 'string', QQ] - 7*b[2]
    x_2 - 7*x_3
    sage: sdp.show()
    Maximization:
    <BLANKLINE>
    Constraints:
    Variables:
      a[1] = x_0 is a continuous variable (min=-oo, max=+oo)
      b[3] = x_1 is a continuous variable (min=-oo, max=+oo)
      a[(4, 'string', Rational Field)] = x_2 is a continuous variable (min=-oo, max=+oo)
      b[2] = x_3 is a continuous variable (min=-oo, max=+oo)


Index of functions and methods
------------------------------

Below are listed the methods of :class:`SemidefiniteProgram`. This module
also implements the :class:`SDPSolverException` exception, as well as the
:class:`SDPVariable` class.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~SemidefiniteProgram.add_constraint`            | Adds a constraint to the ``MixedIntegerLinearProgram``
    :meth:`~SemidefiniteProgram.base_ring`                 | Return the base ring
    :meth:`~SemidefiniteProgram.constraints`               | Returns a list of constraints, as 3-tuples
    :meth:`~SemidefiniteProgram.get_backend`               | Returns the backend instance used
    :meth:`~SemidefiniteProgram.get_values`                | Return values found by the previous call to ``solve()``
    :meth:`~SemidefiniteProgram.linear_constraints_parent` | Return the parent for all linear constraints
    :meth:`~SemidefiniteProgram.linear_function`           | Construct a new linear function
    :meth:`~SemidefiniteProgram.linear_functions_parent`   | Return the parent for all linear functions
    :meth:`~SemidefiniteProgram.new_variable`              | Returns an instance of ``SDPVariable`` associated
    :meth:`~SemidefiniteProgram.number_of_constraints`     | Returns the number of constraints assigned so far
    :meth:`~SemidefiniteProgram.number_of_variables`       | Returns the number of variables used so far
    :meth:`~SemidefiniteProgram.polyhedron`                | Returns the polyhedron defined by the Linear Program
    :meth:`~SemidefiniteProgram.remove_constraint`         | Removes a constraint from self
    :meth:`~SemidefiniteProgram.remove_constraints`        | Remove several constraints
    :meth:`~SemidefiniteProgram.set_objective`             | Sets the objective of the ``MixedIntegerLinearProgram``
    :meth:`~SemidefiniteProgram.set_problem_name`          | Sets the name of the ``MixedIntegerLinearProgram``
    :meth:`~SemidefiniteProgram.show`                      | Displays the ``MixedIntegerLinearProgram`` in a human-readable
    :meth:`~SemidefiniteProgram.solve`                     | Solves the ``MixedIntegerLinearProgram``
    :meth:`~SemidefiniteProgram.solver_parameter`          | Return or define a solver parameter
    :meth:`~SemidefiniteProgram.sum`                       | Efficiently computes the sum of a sequence of LinearFunction elements

AUTHORS:

- Risan (2012/02): added extension for exact computation
"""

#*****************************************************************************
#       Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"

from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.misc.cachefunc import cached_method
from sage.numerical.linear_functions import is_LinearFunction, is_LinearConstraint
from sage.misc.superseded import deprecated_function_alias, deprecation
from sage.matrix.all import Matrix
from sage.matrix.matrix import is_Matrix


cdef class SemidefiniteProgram(SageObject):
    r"""
    The ``SemidefiniteProgram`` class is the link between Sage, linear
    programming (LP) and mixed integer programming (SDP) solvers.

    A Mixed Integer Linear Program (MILP) consists of variables, linear
    constraints on these variables, and an objective function which is to be
    maximised or minimised under these constraints.

    See the :wikipedia:`Linear_programming` for further information on linear
    programming, and the :mod:`MILP module <sage.numerical.sdp>` for its use in
    Sage.

    INPUT:

    - ``solver`` -- selects a solver:

      - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT <http://www.cvxopt.org/>`_
          web site.

      - If ``solver=None`` (default), the default solver is used (see
        :func:`default_sdp_solver`)

    - ``maximization``

      - When set to ``True`` (default), the ``SemidefiniteProgram``
        is defined as a maximization.

      - When set to ``False``, the ``SemidefiniteProgram`` is
        defined as a minimization.

    - ``constraint_generation`` -- Only used when ``solver=None``.

      - When set to ``True``, after solving the ``SemidefiniteProgram``,
        it is possible to add a constraint, and then solve it again.
        The effect is that solvers that do not support this feature will not be
        used.

      - Defaults to ``False``.

    .. WARNING::

        All LP variables are non-negative by default (see :meth:`new_variable`
        and :meth:`set_min`).

    .. SEEALSO::

     - :func:`default_sdp_solver` -- Returns/Sets the default SDP solver.

    EXAMPLES:

    Computation of a maximum stable set in Petersen's graph::

         sage: g = graphs.PetersenGraph()
         sage: p = SemidefiniteProgram(maximization=True)
         sage: b = p.new_variable(binary=True)
         sage: p.set_objective(sum([b[v] for v in g]))
         sage: for (u,v) in g.edges(labels=None):
         ....:     p.add_constraint(b[u] + b[v], max=1)
         sage: p.solve(objective_only=True)
         4.0

    """

    def __init__(self, solver=None, maximization=True,
                 names=tuple()):
        r"""
        Constructor for the ``SemidefiniteProgram`` class.

        INPUT:

        - ``solver`` -- the following solvers should be available through this class:

          - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT <http://www.cvxopt.org/>`_
              web site.

          -If ``solver=None`` (default), the default solver is used (see
           ``default_sdp_solver`` method.

        - ``maximization``

          - When set to ``True`` (default), the ``SemidefiniteProgram``
            is defined as a maximization.
          - When set to ``False``, the ``SemidefiniteProgram`` is
            defined as a minimization.

        - ``constraint_generation`` -- Only used when ``solver=None``.

        - ``names`` -- list/tuple/iterable of string. Default names of
          the SDP variables. Used to enable the ``sdp.<x> =
          SemidefiniteProgram()`` syntax.

        .. SEEALSO::

        - :meth:`default_sdp_solver` -- Returns/Sets the default SDP solver.

        EXAMPLE::

            sage: p = SemidefiniteProgram(maximization=True)

        TESTS:

        Checks that the objects are deallocated without invoking the cyclic garbage
        collector (cf. :trac:`12616`)::

            sage: del p
            sage: def just_create_variables():
            ...       p = SemidefiniteProgram()
            ...       b = p.new_variable(nonnegative=True)
            ...       p.add_constraint(b[3]+b[6] <= 2)
            ...       p.solve()
            sage: C = sage.numerical.sdp.SemidefiniteProgram
            sage: import gc
            sage: _ = gc.collect()  # avoid side effects of other doc tests
            sage: sum([1 for x in gc.get_objects() if isinstance(x,C)])
            0

        We now disable the cyclic garbage collector. Since :trac:`12616` avoids
        a reference cycle, the mixed integer linear program created in
        ``just_create_variables()`` is removed even without the cyclic garbage
        collection::

            sage: gc.disable()
            sage: just_create_variables()
            sage: sum([1 for x in gc.get_objects() if isinstance(x,C)])
            0
            sage: gc.enable()

        Right now the ``nonnegative`` argument is mandatory when a variable is
        created, but later the default will change from nonnegative variables to
        unbounded variables::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(real=True)
            doctest:...: DeprecationWarning: The default value of 'nonnegative' will change, to False instead of True. You should add the explicit 'nonnegative=True'.
            See http://trac.sagemath.org/15521 for details.
            sage: p.get_min(v[0])
            0.0
        """
        self._first_variable_names = list(names)
        from sage.numerical.backends.generic_sdp_backend import get_solver
        self._backend = get_solver(solver=solver)
        if not maximization:
            self._backend.set_sense(-1)

        # Associates an index to the variables
        self._variables = {}


    def linear_functions_parent(self):
        """
        Return the parent for all linear functions

        EXAMPLES::

             sage: p = SemidefiniteProgram()
             sage: p.linear_functions_parent()
             Linear functions over Real Double Field
        """
        if self._linear_functions_parent is None:
            base_ring = self._backend.base_ring()
            from sage.numerical.linear_functions import LinearFunctionsParent
            self._linear_functions_parent = LinearFunctionsParent(base_ring)
        return self._linear_functions_parent

    def linear_constraints_parent(self):
        """
        Return the parent for all linear constraints

        See :mod:`~sage.numerical.linear_functions` for more
        details.

        EXAMPLES::

             sage: p = SemidefiniteProgram()
             sage: p.linear_constraints_parent()
             Linear constraints over Real Double Field
        """
        if self._linear_constraints_parent is None:
            from sage.numerical.linear_functions import LinearConstraintsParent
            LF = self.linear_functions_parent()
            self._linear_constraints_parent = LinearConstraintsParent(LF)
        return self._linear_constraints_parent

    def __call__(self, x):
        """
        Construct a new linear function

        EXAMPLES::

             sage: p = SemidefiniteProgram()
             sage: p.linear_function({1:3, 4:5})
             3*x_1 + 5*x_4

        This is equivalent to::

            sage: p({1:3, 4:5})
            3*x_1 + 5*x_4
        """
        parent = self.linear_functions_parent()
        return parent(x)

    linear_function = __call__

    def _repr_(self):
         r"""
         Returns a short description of the ``SemidefiniteProgram``.

         EXAMPLE::

             sage: p = SemidefiniteProgram()
             sage: v = p.new_variable(nonnegative=True)
             sage: p.add_constraint(v[1] + v[2], max=2)
             sage: print p
             Mixed Integer Program ( maximization, 2 variables, 1 constraints )
         """
         cdef GenericSDPBackend b = self._backend

         return ("Mixed Integer Program "+

                 ( "\"" +self._backend.problem_name()+ "\""
                   if (str(self._backend.problem_name()) != "") else "")+

                 " ( " + ("maximization" if b.is_maximization() else "minimization" ) +

                 ", " + str(b.ncols()) + " variables, " +
                 str(b.nrows()) + " constraints )")

    def __copy__(self):
        r"""
        Returns a copy of the current ``SemidefiniteProgram`` instance.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.add_constraint(v[0] + v[1], max = 10)
            sage: q = copy(p)
            sage: q.number_of_constraints()
            1
        """
        cdef SemidefiniteProgram p = \
            SemidefiniteProgram(solver="cvxopt")
        from copy import copy
        try:
            p._variables = copy(self._variables)
        except AttributeError:
            pass

        try:
            p._default_sdpvariable = self._default_sdpvariable
        except AttributeError:
            pass

        try:
            p._constraints = copy(self._constraints)
        except AttributeError:
            pass

        p._backend = (<GenericSDPBackend> self._backend).copy()
        return p

    def __getitem__(self, v):
        r"""
        Returns the symbolic variable corresponding to the key
        from a default dictionary.

        It returns the element asked, and otherwise creates it.
        If necessary, it also creates the default dictionary.

        This method lets the user define LinearProgram without having to
        define independent dictionaries when it is not necessary for him.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p.set_objective(p['x'] + p['z'])
            sage: p['x']
            x_0
        """

        try:
            return self._default_sdpvariable[v]
        except TypeError:
            self._default_sdpvariable = self.new_variable()
            return self._default_sdpvariable[v]

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        A ring. The coefficients that the chosen solver supports.

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver='cvxopt')
            sage: p.base_ring()
            Real Double Field
            sage: p = SemidefiniteProgram(solver='ppl')
            sage: p.base_ring()
            Rational Field
        """
        return self._backend.base_ring()

    def set_problem_name(self,name):
        r"""
        Sets the name of the ``SemidefiniteProgram``.

        INPUT:

        - ``name`` -- A string representing the name of the
          ``SemidefiniteProgram``.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p.set_problem_name("Test program")
            sage: p
            Semidefinite Program "Test program" ( maximization, 0 variables, 0 constraints )
        """
        self._backend.problem_name(name)

    def new_variable(self, real=False, binary=False, integer=False, nonnegative=None, dim=1, name=""):
        r"""
        Returns an instance of ``SDPVariable`` associated
        to the current instance of ``SemidefiniteProgram``.

        A new variable ``x`` is defined by::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable(nonnegative=True)

        It behaves exactly as an usual dictionary would. It can use any key
        argument you may like, as ``x[5]`` or ``x["b"]``, and has methods
        ``items()`` and ``keys()``.

        INPUT:

        - ``dim`` -- integer. Defines the dimension of the dictionary.
          If ``x`` has dimension `2`, its fields will be of the form
          ``x[key1][key2]``. Deprecated.

        - ``binary, integer, real`` -- boolean. Set one of these
          arguments to ``True`` to ensure that the variable gets the
          corresponding type.

        - ``nonnegative`` -- boolean. Whether the variable should be assumed to
          be nonnegative. Rather useless for the binary type.

        - ``name`` -- string. Associates a name to the variable. This
          is only useful when exporting the linear program to a file
          using ``write_mps`` or ``write_lp``, and has no other
          effect.

        .. SEEALSO::

            - :meth:`set_min`,:meth:`get_min` -- set/get the lower bound of a
              variable. Note that by default, all variables are non-negative.

            - :meth:`set_max`,:meth:`get_max` -- set/get the upper bound of a
              variable.

        EXAMPLE::

            sage: p = SemidefiniteProgram()

         To define two dictionaries of variables, the first being
         of real type, and the second of integer type ::

            sage: x = p.new_variable(real=True, nonnegative=True)
            sage: y = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(x[2] + y[3,5], max=2)
            sage: p.is_integer(x[2])
            False
            sage: p.is_integer(y[3,5])
            True

        An exception is raised when two types are supplied ::

            sage: z = p.new_variable(real = True, integer = True)
            Traceback (most recent call last):
            ...
            ValueError: Exactly one of the available types has to be True

        Unbounded variables::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable(real=True, nonnegative=False)
            sage: y = p.new_variable(integer=True, nonnegative=False)
            sage: p.add_constraint(x[0]+x[3] <= 8)
            sage: p.add_constraint(y[0] >= y[1])
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              x_0 + x_1 <= 8.0
              - x_2 + x_3 <= 0.0
            Variables:
              x_0 is a continuous variable (min=-oo, max=+oo)
              x_1 is a continuous variable (min=-oo, max=+oo)
              x_2 is an integer variable (min=-oo, max=+oo)
              x_3 is an integer variable (min=-oo, max=+oo)

        TESTS:

        Default behaviour (:trac:`15521`)::

            sage: x = p.new_variable(nonnegative=True)
            sage: p.get_min(x[0])
            0.0
        """


        if not name and self._first_variable_names:
            name = self._first_variable_names.pop(0)

        return sdp_variable_parent(self,
                      dim=dim,
                      name=name)

    def _first_ngens(self, n):
        """
        Construct the first `n` SDPVariables.

        This method is used for the generater syntax (see below). You
        probably shouldn't use it for anything else.

        INPUT:

        - ``n`` -- integer. The number of variables to construct.

        OUTPUT:

        A tuple of not necessarily positive :class:`SDPVariable`
        instances.

        EXAMPLES::

            sage: sdp.<a,b> = SemidefiniteProgram()
            sage: a[0] + b[2]
            x_0 + x_1
            sage: sdp.show()
            Maximization:
            <BLANKLINE>
            Constraints:
            Variables:
              a[0] = x_0 is a continuous variable (min=-oo, max=+oo)
              b[2] = x_1 is a continuous variable (min=-oo, max=+oo)
        """
        return tuple(self.new_variable(nonnegative=False) for i in range(n))

    def gen(self, i):
        """
        Return the linear variable `x_i`.

        OUTPUT:

            sage: sdp = SemidefiniteProgram()
            sage: sdp.gen(0)
            x_0
            sage: [sdp.gen(i) for i in range(10)]
            [x_0, x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9]
        """
        return self.linear_functions_parent().gen(i)

    cpdef int number_of_constraints(self):
      r"""
      Returns the number of constraints assigned so far.

      EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)
            sage: p.number_of_constraints()
            2
      """
      return self._backend.nrows()

    cpdef int number_of_variables(self):
      r"""
      Returns the number of variables used so far.

      Note that this is backend-dependent, i.e. we count solver's
      variables rather than user's variables. An example of the latter
      can be seen below: Gurobi converts double inequalities,
      i.e. inequalities like `m <= c^T x <= M`, with `m<M`, into
      equations, by adding extra variables: `c^T x + y = M`, `0 <= y
      <= M-m`.

      EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(p[0] - p[2], max = 4)
            sage: p.number_of_variables()
            2
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)
            sage: p.number_of_variables()
            3
            sage: p = SemidefiniteProgram(solver="glpk")
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.number_of_variables()
            2
            sage: p = SemidefiniteProgram(solver="gurobi")   # optional - Gurobi
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)  # optional - Gurobi
            sage: p.number_of_variables()                          # optional - Gurobi
            3
      """
      return self._backend.ncols()

    def constraints(self, indices = None):
        r"""
        Returns a list of constraints, as 3-tuples.

        INPUT:

        - ``indices`` -- select which constraint(s) to return

            - If ``indices = None``, the method returns the list of all the
              constraints.

            - If ``indices`` is an integer `i`, the method returns constraint
              `i`.

            - If ``indices`` is a list of integers, the method returns the list
              of the corresponding constraints.

        OUTPUT:

        Each constraint is returned as a triple ``lower_bound, (indices,
        coefficients), upper_bound``.  For each of those entries, the
        corresponding linear function is the one associating to variable
        ``indices[i]`` the coefficient ``coefficients[i]``, and `0` to all the
        others.

        ``lower_bound`` and ``upper_bound`` are numerical values.

        EXAMPLE:

        First, let us define a small LP::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)

        To obtain the list of all constraints::

            sage: p.constraints()          # not tested
            [(1.0, ([1, 0], [-1.0, 1.0]), 4.0), (1.0, ([2, 0], [-2.0, 1.0]), None)]

        Or constraint `0` only::

            sage: p.constraints(0)         # not tested
            (1.0, ([1, 0], [-1.0, 1.0]), 4.0)

        A list of constraints containing only `1`::

            sage: p.constraints([1])       # not tested
            [(1.0, ([2, 0], [-2.0, 1.0]), None)]

        TESTS:

        As the ordering of the variables in each constraint depends on the
        solver used, we define a short function reordering it before it is
        printed. The output would look the same without this function applied::

            sage: def reorder_constraint((lb,(ind,coef),ub)):
            ...     d = dict(zip(ind, coef))
            ...     ind.sort()
            ...     return (lb, (ind, [d[i] for i in ind]), ub)

        Running the examples from above, reordering applied::

            sage: p = SemidefiniteProgram(solver = "cvxopt")
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)
            sage: sorted(map(reorder_constraint,p.constraints()))
            [(1.0, ([0, 1], [1.0, -1.0]), 4.0), (1.0, ([0, 2], [1.0, -2.0]), None)]
            sage: reorder_constraint(p.constraints(0))
            (1.0, ([0, 1], [1.0, -1.0]), 4.0)
            sage: sorted(map(reorder_constraint,p.constraints([1])))
            [(1.0, ([0, 2], [1.0, -2.0]), None)]

        """
        from sage.rings.integer import Integer as Integer
        cdef int i
        cdef str s
        cdef GenericSDPBackend b = self._backend

        result = list()

        # If indices is None, we actually want to return all constraints
        if indices is None:
          indices = range(b.nrows())

        # Only one constraint
        if isinstance(indices, int) or isinstance(indices, Integer):
            lb, ub = b.row_bounds(indices)
            return (lb, b.row(indices), ub)

        # List of constraints
        elif isinstance(indices, list):
            for i in indices:
                lb, ub = b.row_bounds(i)
                result.append((lb, b.row(i), ub))

            return result

        # Weird Input
        else:
          raise ValueError, "constraints() requires a list of integers, though it will accommodate None or an integer."

    def polyhedron(self, **kwds):
        r"""
        Returns the polyhedron defined by the Linear Program.

        INPUT:

        All arguments given to this method are forwarded to the constructor of
        the :func:`Polyhedron` class.

        OUTPUT:

        A :func:`Polyhedron` object whose `i`-th variable represents the `i`-th
        variable of ``self``.

        .. warning::

            The polyhedron is built from the variables stored by the LP solver
            (i.e. the output of :meth:`show`). While they usually match the ones
            created explicitely when defining the LP, a solver like Gurobi has
            been known to introduce additional variables to store constraints of
            the type ``lower_bound <= linear_function <= upper bound``. You
            should be fine if you did not install Gurobi or if you do not use it
            as a solver, but keep an eye on the number of variables in the
            polyhedron, or on the output of :meth:`show`. Just in case.

        EXAMPLES:

        A LP on two variables::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] <= 1)
            sage: p.add_constraint(3*p['y'] + p['x'] <= 2)
            sage: P = p.polyhedron(); P
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        3-D Polyhedron::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] + 3*p['z'] <= 1)
            sage: p.add_constraint(2*p['y'] + p['z'] + 3*p['x'] <= 1)
            sage: p.add_constraint(2*p['z'] + p['x'] + 3*p['y'] <= 1)
            sage: P = p.polyhedron(); P
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices

        An empty polyhedron::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] + 3*p['z'] <= 1)
            sage: p.add_constraint(2*p['y'] + p['z'] + 3*p['x'] <= 1)
            sage: p.add_constraint(2*p['z'] + p['x'] + 3*p['y'] >= 2)
            sage: P = p.polyhedron(); P
            The empty polyhedron in QQ^3

        An unbounded polyhedron::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] - p['z'] <= 1)
            sage: P = p.polyhedron(); P
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 3 rays

        A square (see :trac:`14395`) ::

            sage: p = SemidefiniteProgram()
            sage: x,y = p['x'], p['y']
            sage: p.set_min(x,None)
            sage: p.set_min(y,None)
            sage: p.add_constraint( x <= 1 )
            sage: p.add_constraint( x >= -1 )
            sage: p.add_constraint( y <= 1 )
            sage: p.add_constraint( y >= -1 )
            sage: p.polyhedron()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        from copy import copy
        cdef GenericSDPBackend b = self._backend
        cdef int i

        # Constraints
        inequalities = []
        equalities = []
        nvar = self.number_of_variables()
        for lb, (indices, values), ub in self.constraints():
            coeffs = dict(zip(indices, values))
            # Equalities
            if (not lb is None) and lb == ub:
                linear_function = []
                linear_function = [coeffs.get(i,0) for i in range(nvar)]
                linear_function.insert(0,-lb)
                equalities.append(linear_function)
                continue
            # Lower Bound
            if not lb is None:
                linear_function = []
                linear_function = [coeffs.get(i,0) for i in range(nvar)]
                linear_function.insert(0,-lb)
                inequalities.append(linear_function)
            # Upper Bound
            if not ub is None:
                linear_function = []
                linear_function = [-coeffs.get(i,0) for i in range(nvar)]
                linear_function.insert(0,ub)
                inequalities.append(linear_function)

        # Variable bounds
        zero = [0] * nvar
        for 0<= i < nvar:
            lb, ub = b.col_bounds(i)
            # Fixed variable
            if (not lb is None) and lb == ub:
                linear_function = copy(zero)
                linear_function[i] = 1
                linear_function.insert(0,-lb)
                equalities.append(linear_function)
                continue
            # Lower bound
            if not lb is None:
                linear_function = copy(zero)
                linear_function[i] = 1
                linear_function.insert(0,-lb)
                inequalities.append(linear_function)
            # Upper bound
            if not ub is None:
                linear_function = copy(zero)
                linear_function[i] = -1
                linear_function.insert(0,ub)
                inequalities.append(linear_function)
        return Polyhedron(ieqs = inequalities, eqns = equalities)

    def show(self):
        r"""
        Displays the ``SemidefiniteProgram`` in a human-readable
        way.

        EXAMPLES:

        When constraints and variables have names ::

            sage: p = SemidefiniteProgram(solver="cvxopt")
            sage: x = p.new_variable(name="Hey")
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2, name="Constraint_1")
            sage: p.show()
            Maximization:
              Hey[1] + Hey[2]
            Constraints:
              Constraint_1: -3.0 Hey[1] + 2.0 Hey[2] <= 2.0
            Variables:
              Hey[1] = x_0 is a continuous variable (min=0.0, max=+oo)
              Hey[2] = x_1 is a continuous variable (min=0.0, max=+oo)

        Without any names ::

            sage: p = SemidefiniteProgram(solver="cvxopt")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2)
            sage: p.show()
            Maximization:
              x_0 + x_1
            Constraints:
              -3.0 x_0 + 2.0 x_1 <= 2.0
            Variables:
              x_0 is a continuous variable (min=0.0, max=+oo)
              x_1 is a continuous variable (min=0.0, max=+oo)

        """
        cdef int i, j
        cdef GenericSDPBackend b = self._backend

        # inv_variables associates a SDPVariable object to an id
        inv_variables = {}
        for (v, id) in self._variables.iteritems():
            inv_variables[id]=v

        # varid_name associates variables id to names
        varid_name = {}
        for 0<= i < b.ncols():
            s = b.col_name(i)
            varid_name[i] = s if s else 'x_'+str(i)

        ##### Sense and objective function
        print ("Maximization:" if b.is_maximization() else "Minimization:")
        print " ",
        first = True
        for 0<= i< b.ncols():
            c = b.objective_coefficient(i)
            if c == 0:
                continue
            print (("+ " if (not first and c>0) else "") +
                   ("" if c == 1 else ("- " if c == -1 else str(c)+" "))+varid_name[i]
                   ),
            first = False
        if b.obj_constant_term > self._backend.zero(): print "+", b.obj_constant_term
        elif b.obj_constant_term < self._backend.zero(): print "-", -b.obj_constant_term
        print

        ##### Constraints
        print "Constraints:"
        for 0<= i < b.nrows():
            indices, values = b.row(i)
            lb, ub = b.row_bounds(i)
            print " ",
            # Constraint's name
            if b.row_name(i):
                print b.row_name(i)+":",
            # Lower bound
            if lb is not None:
                print str(lb)+" <=",
            first = True
            for j, c in sorted(zip(indices, values)):
                if c == 0:
                    continue
                print (("+ " if (not first and c>0) else "") +
                       ("" if c == 1 else ("- " if c == -1 else (str(c) + " " if first and c < 0 else ("- " + str(abs(c)) + " " if c < 0 else str(c) + " "))))+varid_name[j]
                       ),
                first = False
            # Upper bound
            print ("<= "+str(ub) if ub is not None else "")

        ##### Variables
        print "Variables:"
        for 0<= i < b.ncols():
            var_type = 'a continuous'
            if varid_name[i] == str(self.gen(i)):
                name = varid_name[i]
            else:
                name = '{0} = {1}'.format(varid_name[i], self.gen(i))
            lb, ub = b.col_bounds(i)
            print('  {0} is {1} variable (min={2}, max={3})'.format(
                name, var_type,
                lb if lb is not None else "-oo",
                ub if ub is not None else "+oo"))


    def get_values(self, *lists):
        r"""
        Return values found by the previous call to ``solve()``.

        INPUT:

        - Any instance of ``SDPVariable`` (or one of its elements),
          or lists of them.

        OUTPUT:

        - Each instance of ``SDPVariable`` is replaced by a dictionary
          containing the numerical values found for each
          corresponding variable in the instance.
        - Each element of an instance of a ``SDPVariable`` is replaced
          by its corresponding numerical value.

        .. NOTE::

            While a variable may be declared as binary or integer, its value as
            returned by the solver is of type ``float``.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable(nonnegative=True)
            sage: y = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[3] + 3*y[2,9] + x[5])
            sage: p.add_constraint(x[3] + y[2,9] + 2*x[5], max=2)
            sage: p.solve()
            6.0

        To return  the optimal value of ``y[2,9]``::

            sage: p.get_values(y[2,9])
            2.0

        To get a dictionary identical to ``x`` containing optimal
        values for the corresponding variables ::

            sage: x_sol = p.get_values(x)
            sage: x_sol.keys()
            [3, 5]

        Obviously, it also works with variables of higher dimension::

            sage: y_sol = p.get_values(y)

        We could also have tried ::

            sage: [x_sol, y_sol] = p.get_values(x, y)

        Or::

            sage: [x_sol, y_sol] = p.get_values([x, y])

        TESTS:

        When 'dim' will be removed, also remove from this function the code that
        uses it::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable(dim=2)
            doctest:...: DeprecationWarning: The 'dim' argument will soon disappear. Fortunately variable[1,2] is easier to use than variable[1][2]
            See http://trac.sagemath.org/15489 for details.
            sage: p.add_constraint(b[1][2] +  b[2][3] == 0)
            sage: _ = p.solve()
            sage: _ = p.get_values(b)
        """
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
        r"""
        Sets the objective of the ``SemidefiniteProgram``.

        INPUT:

        - ``obj`` -- A linear function to be optimized.
          ( can also be set to ``None`` or ``0`` when just
          looking for a feasible solution )

        EXAMPLE:

        Let's solve the following linear program::

            Maximize:
              x + 5 * y
            Constraints:
              x + 0.2 y       <= 4
              1.5 * x + 3 * y <= 4
            Variables:
              x is Real (min = 0, max = None)
              y is Real (min = 0, max = None)

        This linear program can be solved as follows::

            sage: p = SemidefiniteProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 2/10*x[2], max=4)
            sage: p.add_constraint(1.5*x[1]+3*x[2], max=4)
            sage: round(p.solve(),5)
            6.66667
            sage: p.set_objective(None)
            sage: _ = p.solve()
        """
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

    def add_constraint(self, linear_function, name=None):
        r"""
        Adds a constraint to the ``SemidefiniteProgram``.

        INPUT:

        - ``linear_function`` -- Two different types of arguments are possible:
            - A linear function. In this case, arguments ``min`` or ``max``
              have to be specified.
            - A linear constraint of the form ``A <= B``, ``A >= B``,
              ``A <= B <= C``, ``A >= B >= C`` or ``A == B``. In this
              case, arguments ``min`` and ``max`` will be ignored.
        - ``name`` -- A name for the constraint.

        EXAMPLE:

        Consider the following linear program::

            Maximize:
              x + 5 * y
            Constraints:
              x + 0.2 y       <= 4
              1.5 * x + 3 * y <= 4
            Variables:
              x is Real (min = 0, max = None)
              y is Real (min = 0, max = None)

        It can be solved as follows::

            sage: p = SemidefiniteProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2], max=4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2], max=4)
            sage: round(p.solve(),6)
            6.666667

        To add a constrain we give ``add_constraint`` this
        very expression::

            sage: p.add_constraint( x[5] + 3*x[7] <= x[6] + 3 )

        One can also define double-bounds or equality using symbols
        ``<=``, ``>=`` and ``==``::

            sage: p.add_constraint( x[5] + 3*x[7] == x[6] + 3 )
            sage: p.add_constraint( x[5] + 3*x[7] <= x[6] + 3 <= x[8] + 27 )

        The previous program can be rewritten::

            sage: p = SemidefiniteProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2] <= 4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2] <= 4)
            sage: round(p.solve(), 5)
            6.66667

        TESTS:

        Complex constraints::

            sage: p = SemidefiniteProgram(solver = "cvxopt")
            sage: b = p.new_variable(nonnegative=True)
            sage: p.add_constraint( b[8] - b[15] <= 3*b[8] + 9)
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              -2.0 x_0 - x_1 <= 9.0
            Variables:
              x_0 is a continuous variable (min=0.0, max=+oo)
              x_1 is a continuous variable (min=0.0, max=+oo)

        Empty constraint::

            sage: p=SemidefiniteProgram()
            sage: p.add_constraint(sum([]),min=2)


        """
        if linear_function is 0:
            return

        from sage.numerical.linear_tensor_constraints import is_LinearTensorConstraint
        from sage.numerical.linear_tensor import is_LinearTensor

        if is_LinearTensorConstraint(linear_function):
            c = linear_function
            if c.is_equation():
                self.add_constraint(c.lhs()-c.rhs(), name=name)
                self.add_constraint(-c.lhs()+c.rhs(), name=name)
            else:
                self.add_constraint(c.lhs()-c.rhs(), name=name)

        elif is_LinearTensor(linear_function):
            l = linear_function.dict().items()
            l.sort()
            self._backend.add_linear_constraint(l, name)

        else:
            raise ValueError('argument must be a linear function or constraint, got '+str(linear_function))

    def get_matrix(self):
        return self._backend.get_matrix()


    def solve(self, log=None, objective_only=False):
        r"""
        Solves the ``SemidefiniteProgram``.

        INPUT:

        - ``log`` -- integer (default: ``None``) The verbosity level. Indicates
          whether progress should be printed during computation. The solver is
          initialized to report no progress.

        - ``objective_only`` -- Boolean variable.

          - When set to ``True``, only the objective function is returned.
          - When set to ``False`` (default), the optimal numerical values
            are stored (takes computational time).

        OUTPUT:

        The optimal value taken by the objective function.

        .. WARNING::

            By default, all variables of a LP are assumed to be
            non-negative. See :meth:`set_min` to change it.

        EXAMPLES:

        Consider the following linear program::

            Maximize:
              x + 5 * y
            Constraints:
              x + 0.2 y       <= 4
              1.5 * x + 3 * y <= 4
            Variables:
              x is Real (min = 0, max = None)
              y is Real (min = 0, max = None)

        This linear program can be solved as follows::

            sage: p = SemidefiniteProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2], max=4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2], max=4)
            sage: round(p.solve(),6)
            6.666667
            sage: x = p.get_values(x)
            sage: round(x[1],6)
            0.0
            sage: round(x[2],6)
            1.333333

         Computation of a maximum stable set in Petersen's graph::

            sage: g = graphs.PetersenGraph()
            sage: p = SemidefiniteProgram(maximization=True)
            sage: b = p.new_variable(nonnegative=True)
            sage: p.set_objective(sum([b[v] for v in g]))
            sage: for (u,v) in g.edges(labels=None):
            ...       p.add_constraint(b[u] + b[v], max=1)
            sage: p.set_binary(b)
            sage: p.solve(objective_only=True)
            4.0

        Constraints in the objective function are respected::

            sage: p = SemidefiniteProgram()
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve()
            9.0
        """
        if log is not None: self._backend.set_verbosity(log)
        self._backend.solve()
        return self._backend.get_objective_value()


    def solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        The solver parameters are by essence solver-specific, which
        means their meaning heavily depends on the solver used.

        (If you do not know which solver you are using, then you use
        use cvxopt).

        Aliases:

        Very common parameters have aliases making them
        solver-independent. For example, the following::

            sage: p = SemidefiniteProgram(solver = "cvxopt")
            sage: p.solver_parameter("timelimit", 60)

        Sets the solver to stop its computations after 60 seconds, and
        works with GLPK, CPLEX and Gurobi.

            - ``"timelimit"`` -- defines the maximum time spent on a
              computation. Measured in seconds.

        Solver-specific parameters:

            - GLPK : We have implemented very close to comprehensive coverage of
              the cvxopt solver parameters for the simplex and integer
              optimization methods. For details, see the documentation of
              :meth:`cvxoptBackend.solver_parameter
              <sage.numerical.backends.glpk_backend.cvxoptBackend.solver_parameter>`.


        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        EXAMPLE::

            sage: p = SemidefiniteProgram(solver = "cvxopt")
            sage: p.solver_parameter("timelimit", 60)
            sage: p.solver_parameter("timelimit")
            60.0
        """
        if value is None:
            return self._backend.solver_parameter(name)
        else:
            self._backend.solver_parameter(name, value)

    cpdef sum(self, L):
        r"""
        Efficiently computes the sum of a sequence of
        :class:`~sage.numerical.linear_functions.LinearFunction` elements

        INPUT:

        - ``sdp`` -- the :class:`SemidefiniteProgram` parent.

        - ``L`` -- list of
          :class:`~sage.numerical.linear_functions.LinearFunction` instances.

        .. NOTE::

            The use of the regular ``sum`` function is not recommended
            as it is much less efficient than this one

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(nonnegative=True)

        The following command::

            sage: s = p.sum([v[i] for i in xrange(90)])

        is much more efficient than::

            sage: s = sum([v[i] for i in xrange(90)])
        """
        d = {}
        for v in L:
            for id,coeff  in v.iteritems():
                d[id] = coeff + d.get(id,0)
        return self.linear_functions_parent()(d)

    def get_backend(self):
        r"""
        Returns the backend instance used.

        This might be useful when acces to additional functions provided by
        the backend is needed.

        EXAMPLE:

        This example uses the simplex algorthm and prints information::

            sage: p = SemidefiniteProgram(solver="cvxopt")
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: b = p.get_backend()
            sage: b.solver_parameter("simplex_or_intopt", "simplex_only")
            sage: b.solver_parameter("verbosity_simplex", "GLP_MSG_ALL")
            sage: p.solve()  # tol 0.00001
            cvxopt Simplex Optimizer, v4.44
            2 rows, 2 columns, 4 non-zeros
            *     0: obj =   7.000000000e+00  infeas =  0.000e+00 (0)
            *     2: obj =   9.400000000e+00  infeas =  0.000e+00 (0)
            OPTIMAL SOLUTION FOUND
            9.4
        """
        return self._backend


class SDPSolverException(RuntimeError):
    r"""
    Exception raised when the solver fails.
    """

    def __init__(self, value):
        r"""
        Constructor for ``SDPSolverException``.

        ``SDPSolverException`` is the exception raised when the solver fails.

        EXAMPLE::

            sage: from sage.numerical.sdp import SDPSolverException
            sage: SDPSolverException("Error")
            SDPSolverException()

        TESTS:

        No continuous solution::

            sage: p=SemidefiniteProgram(solver="cvxopt")
            sage: v=p.new_variable(nonnegative=True)
            sage: p.add_constraint(v[0],max=5.5)
            sage: p.add_constraint(v[0],min=7.6)
            sage: p.set_objective(v[0])


        No integer solution::

            sage: p=SemidefiniteProgram(solver="cvxopt")
            sage: v=p.new_variable(nonnegative=True)
            sage: p.add_constraint(v[0],max=5.6)
            sage: p.add_constraint(v[0],min=5.2)
            sage: p.set_objective(v[0])
            sage: p.set_integer(v)

        Tests of cvxopt's Exceptions::

            sage: p.solve()
            Traceback (most recent call last):
            ...
            SDPSolverException: 'cvxopt : Solution is undefined'
        """
        self.value = value

    def __str__(self):
        r"""
        Returns the value of the instance of ``SDPSolverException``.

        EXAMPLE::

            sage: from sage.numerical.sdp import SDPSolverException
            sage: e = SDPSolverException("Error")
            sage: print e
            'Error'
        """
        return repr(self.value)


cdef class SDPVariable(Element):
    r"""
    ``SDPVariable`` is a variable used by the class
    ``SemidefiniteProgram``.

    .. warning::

        You should not instantiate this class directly. Instead, use
        :meth:`SemidefiniteProgram.new_variable`.
    """

    def __init__(self, parent, sdp, dim, name):
        r"""
        Constructor for ``SDPVariable``.

        INPUT:

        - ``parent`` -- :class:`SDPVariableParent`. The parent of the
          SDP variable.

        - ``sdp`` -- :class:`SemidefiniteProgram`. The
          underlying linear program.

        - ``dim`` -- the integer defining the definition of the variable.

        - ``name`` -- A name for the ``SDPVariable``.

        - ``lower_bound``, ``upper_bound`` -- lower bound and upper
          bound on the variable. Set to ``None`` to indicate that the
          variable is unbounded.

        For more informations, see the method
        ``SemidefiniteProgram.new_variable``.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p.new_variable(nonnegative=True)
            SDPVariable of dimension 1.
        """
        super(SDPVariable, self).__init__(parent)
        self._dim = dim
        self._dict = {}
        self._p = sdp
        self._name = name

        if dim > 1:
            from sage.misc.superseded import deprecation
            deprecation(15489, "The 'dim' argument will soon disappear. "+
                        "Fortunately variable[1,2] is easier to use than "+
                        "variable[1][2]")

    def __getitem__(self, i):
        r"""
        Returns the symbolic variable corresponding to the key.

        Returns the element asked, otherwise creates it.
        (When depth>1, recursively creates the variables).

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v[0]
            x_0

        TESTS:

        This function contains "dim" code that will have to be removed::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable(binary=True, dim=2)
        """
        cdef int j
        if i in self._dict:
            return self._dict[i]
        elif self._dim == 1:
            zero = self._p._backend.zero()
            name = self._name + "[" + str(i) + "]" if self._name else None
            j = self._p._backend.add_variable(
                obj=zero,
                name=name,
            )
            v = self._p.linear_function({j : 1})
            self._p._variables[v] = j
            self._dict[i] = v
            return v
        else:
            name = self._name + "[" + str(i) + "]" if self._name else ''
            self._dict[i] = self.parent().element_class(
                self.parent(),
                self._p,
                dim=self._dim - 1,
                name=name,
            )
            return self._dict[i]


    def _repr_(self):
        r"""
        Returns a representation of self.

        EXAMPLE::

            sage: p=SemidefiniteProgram()
            sage: v=p.new_variable(dim=3)
            sage: v
            SDPVariable of dimension 3.
            sage: v[2][5][9]
            x_0
            sage: v
            SDPVariable of dimension 3.
        """
        return "SDPVariable"

    def keys(self):
        r"""
        Returns the keys already defined in the dictionary.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v.keys()
            [0, 1]
        """
        return self._dict.keys()

    def items(self):
        r"""
        Returns the pairs (keys,value) contained in the dictionary.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v.items()
            [(0, x_0), (1, x_1)]
        """
        return self._dict.items()

    def depth(self):
        r"""
        Returns the current variable's depth.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v.depth()
            1
        """
        return self._dim

    def values(self):
        r"""
        Returns the symbolic variables associated to the current dictionary.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v.values()
            [x_0, x_1]
        """
        return self._dict.values()

    cdef _matrix_rmul_impl(self, m):
        """
        Implement the action of a matrix multiplying from the right.
        """
        result = dict()
        for i, row in enumerate(m.rows()):
            x = self[i]
            assert len(x.dict()) == 1
            x_index = x.dict().keys()[0]
            result[x_index] = row
        from sage.modules.free_module import FreeModule
        V = FreeModule(self._p.base_ring(), m.ncols())
        T = self._p.linear_functions_parent().tensor(V)
        return T(result)

    cdef _matrix_lmul_impl(self, m):
        """
        Implement the action of a matrix multiplying from the left.
        """
        result = dict()
        for i, col in enumerate(m.columns()):
            x = self[i]
            assert len(x.dict()) == 1
            x_index = x.dict().keys()[0]
            result[x_index] = col
        from sage.modules.free_module import FreeModule
        V = FreeModule(self._p.base_ring(), m.nrows())
        T = self._p.linear_functions_parent().tensor(V)
        return T(result)

    cpdef _acted_upon_(self, mat, bint self_on_left):
        """
        Act with matrices on SDPVariables.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()
            sage: m = matrix([[1,2], [3,4]])
            sage: v * m
            (1.0, 2.0)*x_0 + (3.0, 4.0)*x_1
            sage: m * v
            (1.0, 3.0)*x_0 + (2.0, 4.0)*x_1
        """
        from sage.matrix.matrix import is_Matrix
        if is_Matrix(mat):
            return self._matrix_rmul_impl(mat) if self_on_left else self._matrix_lmul_impl(mat)


cdef class SDPVariableParent(Parent):
    """
    Parent for :class:`SDPVariable`.

    .. warning::

        This class is for internal use. You should not instantiate it
        yourself. Use :meth:`SemidefiniteProgram.new_variable`
        to generate sdp variables.
    """

    Element = SDPVariable

    def _repr_(self):
        r"""
        Return representation of self.

        OUTPUT:

        String.

        EXAMPLES::

            sage: sdp.<v> = SemidefiniteProgram()
            sage: v.parent()
            Parent of SDPVariables
        """
        return 'Parent of SDPVariables'

    def _an_element_(self):
        """
        Construct a SDP variable.

        OUTPUT:

        This is required for the coercion framework. We raise a
        ``TypeError`` to abort search for any coercion to another
        parent for binary operations. The only interesting operations
        involving :class:`SDPVariable` elements are actions by
        matrices.

        EXAMPLES::

            sage: sdp.<x> = SemidefiniteProgram()
            sage: parent = x.parent()
            sage: parent.an_element()    # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: disallow coercion
        """
        raise TypeError('disallow coercion')

    def _element_constructor_(self, sdp, dim=1, name=""):
        """
        The Element constructor

        INPUT/OUTPUT:

        See :meth:`SDPVariable.__init__`.

        EXAMPLES::

            sage: sdp = SemidefiniteProgram()
            sage: sdp.new_variable()    # indirect doctest
            SDPVariable of dimension 1.
        """
        return self.element_class(self, sdp, dim, name)


sdp_variable_parent = SDPVariableParent()


def Sum(x):
    """
    Only for legacy support, use :meth:`SemidefiniteProgram.sum` instead.

    EXAMPLES::

        sage: from sage.numerical.sdp import Sum
        sage: Sum([])
        doctest:...: DeprecationWarning: use SemidefiniteProgram.sum() instead
        See http://trac.sagemath.org/13646 for details.

        sage: p = SemidefiniteProgram()
        sage: x = p.new_variable(nonnegative=True)
        sage: Sum([ x[0]+x[1], x[1]+x[2], x[2]+x[3] ])   # deprecation is only shown once
        x_0 + 2*x_1 + 2*x_2 + x_3
    """
    from sage.misc.superseded import deprecation
    deprecation(13646, 'use SemidefiniteProgram.sum() instead')
    if not x:
        return None
    parent = x[0].parent()
    return parent.sum(x)
