r"""
CVXOPT Backend


AUTHORS:

- Ingolfur Edvardsson (2014-05)        : initial implementation

"""

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################



cdef class CVXOPTBackend(GenericBackend):
    cdef list objective_function #c_matrix
    cdef list G_matrix
    cdef list h_matrix
    cdef str prob_name
    cdef int is_maximize

    cdef list row_lower_bound
    cdef list row_upper_bound
    cdef list col_lower_bound
    cdef list col_upper_bound

    cdef list row_name_var
    cdef list col_name_var
    cdef dict answer

    def __cinit__(self, maximization = True):
        """
        Cython constructor

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")

        """

        self.objective_function = [] #c_matrix in the example for cvxopt
        self.G_matrix = []
        self.h_matrix = []
        self.prob_name = None
        self.obj_constant_term = 0
        self.is_maximize = 0

        self.row_lower_bound = []
        self.row_upper_bound = []
        self.col_lower_bound = []
        self.col_upper_bound = []

        self.row_name_var = []
        self.col_name_var = []


        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

        #self.obj_constant_term = 0.0

    def __dealloc__(self):
        r"""
        Destructor function
        """
        #needed?
        #del self.solver

    cpdef int add_variable(self, lower_bound=None, upper_bound=None, binary=False, continuous=True, integer=False, obj=None, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        INPUT:

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(binary=True)
            1
            sage: p.add_variable(lower_bound=-2.0, integer=True)
            2
            sage: p.add_variable(continuous=True, integer=True)       # optional - CVXOPT
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1.0)                    # optional - CVXOPT
            3
            sage: p.col_name(3)                                       # optional - CVXOPT
            'x'
            sage: p.objective_coefficient(3)                          # optional - CVXOPT
            1.0
        """
        if obj == None:
            obj = 0.0
        self.G_matrix.append([0 for i in range(self.nrows())])
        self.col_lower_bound.append(lower_bound)
        self.col_upper_bound.append(upper_bound)
        self.objective_function.append(obj)
        self.col_name_var.append(name)
        return len(self.objective_function) - 1


    cpdef int add_variables(self, int n, lower_bound=None, upper_bound=None, binary=False, continuous=True, integer=False, obj=None, names=None) except -1:
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        INPUT:

        - ``n`` - the number of new variables (must be > 0)

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of all variables in the objective function (default: 0.0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, names=['a','b'])
            6
        """
        for i in range(n):
            self.add_variable()
        return len(self.objective_function) - 1;


    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            *  -1  Continuous

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")   # optional - CVXOPT
            sage: p.ncols()                                        # optional - CVXOPT
            0
            sage: p.add_variable()                                  # optional - CVXOPT
            1
            sage: p.set_variable_type(0,1)                          # optional - CVXOPT
            sage: p.is_variable_integer(0)                          # optional - CVXOPT
            True
        """
        pass

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        if sense == 1:
            self.is_maximize = 1
        else:
            self.is_maximize = 0

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
        """
        if coeff is not None:
            self.objective_function[variable] = float(coeff);
        else:
            return self.objective_function[variable]

    cpdef set_objective(self, list coeff, d = 0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: map(lambda x :p.objective_coefficient(x), range(5)) #optional - CVXOPT
            [1.0, 1.0, 2.0, 1.0, 3.0]

        Constants in the objective function are respected::

            sage: p = MixedIntegerLinearProgram(solver='CVXOPT')
            sage: x=p.new_variable(nonnegative=True)[0]
            sage: y=p.new_variable(nonnegative=True)[0]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve() #optional - CVXOPT
            9.0
        """
        for i in range(len(coeff)):
            self.objective_function[i] = coeff[i];
        obj_constant_term = d;

    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.set_verbosity(2)
        """
        pass

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint.

        INPUT::

        - ``i`` -- index of the constraint to remove.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_constraint(p[0] + p[1], max = 10)           # optional - CVXOPT
            sage: p.remove_constraint(0)                            # optional - CVXOPT
        """
        raise NotImplementedError()

    cpdef remove_constraints(self, constraints):
        r"""
        Remove several constraints.

        INPUT:

        - ``constraints`` -- an iterable containing the indices of the rows to remove.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_constraint(p[0] + p[1], max = 10)           # optional - CVXOPT
            sage: p.remove_constraints([0])                         # optional - CVXOPT
        """
        raise NotImplementedError()

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``name`` - an optional name for this row (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                                       # optional - CVXOPT
            sage: p.add_variables(5)                                                      # optional - CVXOPT
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0) # optional - CVXOPT
            sage: p.row(0)                                                                # optional - CVXOPT
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])                                          # optional - CVXOPT
            sage: p.row_bounds(0)                                                         # optional - CVXOPT
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo') # optional - CVXOPT
            sage: p.row_name(-1)                                                          # optional - CVXOPT
            "foo"
        """
        for column in self.G_matrix:
            column.append(0)
        for a in coefficients:
            self.G_matrix[a[0]][-1] = a[1]

        self.row_lower_bound.append(lower_bound)
        self.row_upper_bound.append(upper_bound)
        self.row_name_var.append(name)

    cpdef add_col(self, list indices, list coeffs):
        """
        Add a column.

        INPUT:

        - ``indices`` (list of integers) -- this list constains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` (list of real values) -- associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the ith entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the ith entry in ``indices``.

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
        """
        column = []
        for i in range(len(indices)):
            column.append(0)

        for i in range(len(indices)):
            column[indices[i]] = coeffs[i]

        self.G_matrix.append(column)

        self.col_lower_bound.append(None)
        self.col_upper_bound.append(None)
        self.objective_function.append(0)
        self.col_name_var.append(None)

    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")   # optional - CVXOPT
            sage: p.add_variables(5)                                # optional - CVXOPT
            5
            sage: p.add_linear_constraints(5, None, 2)          # optional - CVXOPT
            sage: p.row(4)                                      # optional - CVXOPT
            ([], [])
            sage: p.row_bounds(4)                               # optional - CVXOPT
            (None, 2.0)
        """
        for i in range(number):
            for j in range(len(self.objective_function)):
                self.G_matrix[j].append(0)
            self.row_lower_bound.append(lower_bound)
            self.row_upper_bound.append(upper_bound)
            if names is not None:
                self.row_name_var.append(names)
            else:
                self.row_name_var.append(None)


    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT") # optional - CVXOPT
            sage: p.add_linear_constraints(5, 0, None)             # optional - CVXOPT
            sage: p.add_col(range(5), range(5))                    # optional - CVXOPT
            sage: p.solve()                                        # optional - CVXOPT
            0
            sage: p.objective_coefficient(0,1)                     # optional - CVXOPT
            sage: p.solve()                                        # optional - CVXOPT
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """



        """
            sage: p = MixedIntegerLinearProgram(solver = "cvxopt")
            sage: x=p.new_variable(nonnegative=True)[0]
            sage: y=p.new_variable(nonnegative=True)[0]
            sage: z=p.new_variable(nonnegative=True)[0]
            sage: p.set_objective(x + y + 3*z)
            sage: p.add_constraint(x + 2*y <= 4)
            sage: p.add_constraint(5*z - y <= 8)
            sage: round(p.solve(), 2)
            8.8

        """
        from cvxopt import matrix, solvers
        #multiply by -1 if necessary
        #print str("G_matrix is: " ) + str(self.G_matrix)
        #print str("lower bound is: " ) + str(self.row_lower_bound)
        #print str("upper bound is: " ) + str(self.row_upper_bound)
        for eq_index in range(self.nrows()):
            if self.row_lower_bound[eq_index] != None:
                #switch between upper and lower bounds
                self.row_upper_bound[eq_index] = -1 * self.row_lower_bound[eq_index]
                self.row_lower_bound[eq_index] = None
                #-1 to all the elements in G
                for j in range(self.ncols() ):
                    self.G_matrix[j] = -1 * self.G_matrix[j]

        G = []
        for col in self.G_matrix:
            tempcol = []
            for i in range(len(col)):
                tempcol.append( float(col[i]))
            G.append(tempcol)
        G = matrix(G)
        c = matrix([float(e) for e in self.objective_function])
        h = matrix([float(e) for e in self.row_upper_bound])

        #to hide the output from cvxopt
        import sys, StringIO
        actualstdout = sys.stdout
        sys.stdout = StringIO.StringIO()
        #solvers comes from the cvxopt library
        self.answer = solvers.lp(c,G,h)
        sys.stdout = actualstdout
        sys.stdout.flush()
        return 0


    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT") # optional - CVXOPT
            sage: p.add_variables(2)                               # optional - CVXOPT
            2
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3) # optional - CVXOPT
            sage: p.set_objective([2, 5])                          # optional - CVXOPT
            sage: p.solve()                                        # optional - CVXOPT
            0
            sage: p.get_objective_value()                          # optional - CVXOPT
            7.5
            sage: p.get_variable_value(0)                          # optional - CVXOPT
            0.0
            sage: p.get_variable_value(1)                          # optional - CVXOPT
            1.5
        """
        sum = self.obj_constant_term
        i = 0
        for v in self.objective_function:
            sum += v * float(self.answer['x'][i])
            i+=1
        return sum

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT") # optional - CVXOPT
            sage: p.add_variables(2)                              # optional - CVXOPT
            2
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3) # optional - CVXOPT
            sage: p.set_objective([2, 5])                         # optional - CVXOPT
            sage: p.solve()                                       # optional - CVXOPT
            0
            sage: p.get_objective_value()                         # optional - CVXOPT
            7.5
            sage: p.get_variable_value(0)                         # optional - CVXOPT
            0.0
            sage: p.get_variable_value(1)                         # optional - CVXOPT
            1.5
        """
        return self.answer['x'][variable]


    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.ncols()                                       # optional - CVXOPT
            0
            sage: p.add_variables(2)                               # optional - CVXOPT
            2
            sage: p.ncols()                                       # optional - CVXOPT
            2
        """

        return len(self.objective_function)

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT") # optional - CVXOPT
            sage: p.nrows()                                        # optional - CVXOPT
            0
            sage: p.add_linear_constraints(2, 2.0, None)         # optional - CVXOPT
            sage: p.nrows()                                      # optional - CVXOPT
            2
        """
        if (len(self.G_matrix) == 0): return 0
        return len(self.G_matrix[0])


    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT") # optional - CVXOPT
            sage: p.is_maximization()                             # optional - CVXOPT
            True
            sage: p.set_sense(-1)                             # optional - CVXOPT
            sage: p.is_maximization()                             # optional - CVXOPT
            False
        """
        if self.is_maximize == 1:
            return 1
        else:
            return 0

    cpdef problem_name(self, char * name = NULL):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")   # optional - CVXOPT
            sage: p.problem_name("There once was a french fry") # optional - CVXOPT
            sage: print p.get_problem_name()                        # optional - CVXOPT
            There once was a french fry
        """
        if name == NULL:
            return self.name
        self.name = str(<bytes>name)

    cpdef write_lp(self, char * name):
        """
        Write the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_variables(2)                               # optional - CVXOPT
            2
            sage: p.add_linear_constraint([(0, 1], (1, 2)], None, 3) # optional - CVXOPT
            sage: p.set_objective([2, 5])                          # optional - CVXOPT
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))            # optional - CVXOPT
        """
        raise NotImplementedError()

    cpdef write_mps(self, char * name, int modern):
        """
        Write the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_variables(2)                               # optional - CVXOPT
            2
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3) # optional - CVXOPT
            sage: p.set_objective([2, 5])                          # optional - CVXOPT
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))            # optional - CVXOPT
        """
        raise NotImplementedError()

    cpdef row(self, int i):
        """
        Return a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)                              # optional - CVXOPT
            (2.0, 2.0)
        """
        coeff = []
        idx = []
        for j in range(len(self.G_matrix)):
            if self.G_matrix[i] != 0:
                idx.append(j)
                coeff.append(self.G_matrix[i])
        coeff.append(self.h_matrix[i])
        return (idx, coeff)



    cpdef row_bounds(self, int index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_variables(5)                               # optional - CVXOPT
            5
            sage: p.add_linear_constraint(range(5), range(5), 2, 2) # optional - CVXOPT
            sage: p.row(0)                                     # optional - CVXOPT
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)                              # optional - CVXOPT
            (2.0, 2.0)
        """
        raise NotImplementedError()

    cpdef col_bounds(self, int index):
        """
        Return the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_variable()                                 # optional - CVXOPT
            1
            sage: p.col_bounds(0)                              # optional - CVXOPT
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                 # optional - CVXOPT
            sage: p.col_bounds(0)                              # optional - CVXOPT
            (0.0, 5.0)
        """
        raise NotImplementedError()

    cpdef bint is_variable_binary(self, int index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.ncols()                                       # optional - CVXOPT
            0
            sage: p.add_variable()                                 # optional - CVXOPT
            1
            sage: p.set_variable_type(0,0)                         # optional - CVXOPT
            sage: p.is_variable_binary(0)                          # optional - CVXOPT
            True

        """
        raise NotImplementedError()

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.ncols()                                       # optional - CVXOPT
            0
            sage: p.add_variable()                                 # optional - CVXOPT
            1
            sage: p.set_variable_type(0,1)                         # optional - CVXOPT
            sage: p.is_variable_integer(0)                         # optional - CVXOPT
            True
        """
        raise NotImplementedError()

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.ncols()                                       # optional - CVXOPT
            0
            sage: p.add_variable()                                 # optional - CVXOPT
            1
            sage: p.is_variable_continuous(0)                      # optional - CVXOPT
            True
            sage: p.set_variable_type(0,1)                         # optional - CVXOPT
            sage: p.is_variable_continuous(0)                      # optional - CVXOPT
            False

        """
        raise NotImplementedError()

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_linear_constraints(1, 2, None, name="Empty constraint 1")  # optional - CVXOPT
            sage: p.row_name(0)                                     # optional - CVXOPT
            'Empty constraint 1'

        """
        raise NotImplementedError()

    cpdef col_name(self, int index):
        """
        Return the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_variable(name="I am a variable")            # optional - CVXOPT
            1
            sage: p.col_name(0)                                     # optional - CVXOPT
            'I am a variable'
        """
        raise NotImplementedError()

    cpdef variable_upper_bound(self, int index, value = None):
        """
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_variable()                                 # optional - CVXOPT
            1
            sage: p.col_bounds(0)                              # optional - CVXOPT
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                 # optional - CVXOPT
            sage: p.col_bounds(0)                              # optional - CVXOPT
            (0.0, 5.0)
        """
        raise NotImplementedError()

    cpdef variable_lower_bound(self, int index, value = None):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.add_variable()                                 # optional - CVXOPT
            1
            sage: p.col_bounds(0)                              # optional - CVXOPT
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)                 # optional - CVXOPT
            sage: p.col_bounds(0)                              # optional - CVXOPT
            (5.0, None)
        """
        raise NotImplementedError()

    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        .. NOTE::

           The list of available parameters is available at
           :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solver_parameter`.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")  # optional - CVXOPT
            sage: p.solver_parameter("timelimit")                   # optional - CVXOPT
            sage: p.solver_parameter("timelimit", 60)               # optional - CVXOPT
            sage: p.solver_parameter("timelimit")                   # optional - CVXOPT
        """
        raise NotImplementedError()
