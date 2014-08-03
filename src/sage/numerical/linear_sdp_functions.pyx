"""
Linear SDP Functions and Constraints

This module implements linear functions (see :class:`LinearSDPFunction`)
in formal variables and chained (in)equalities between them (see
:class:`LinearSDPConstraint`). By convention, these are always written as
either equalities or less-or-equal. For example::

    sage: p = SemidefiniteProgram()
    sage: x = p.new_variable()
    sage: f = 1 + x[1] + 2*x[2];  f     #  a linear function
    1 + x_0 + 2*x_1
    sage: type(f)
    <type 'sage.numerical.linear_sdp_functions.LinearSDPFunction'>

    sage: c = (0 <= f);  c    # a constraint
    0 <= 1 + x_0 + 2*x_1
    sage: type(c)
    <type 'sage.numerical.linear_sdp_functions.LinearSDPConstraint'>

Note that you can use this module without any reference to linear
programming, it only implements linear functions over a base ring and
constraints. However, for ease of demonstration we will always
construct them out of linear programs (see
:mod:`~sage.numerical.mip`).

Constraints can be equations or (non-strict) inequalities. They can be
chained::

    sage: p = SemidefiniteProgram()
    sage: x = p.new_variable()
    sage: x[0] == x[1] == x[2] == x[3]
    x_0 == x_1 == x_2 == x_3

    sage: ieq_01234 = x[0] <= x[1] <= x[2] <= x[3] <= x[4]
    sage: ieq_01234
    x_0 <= x_1 <= x_2 <= x_3 <= x_4

If necessary, the direction of inequality is flipped to always write
inqualities as less or equal::

    sage: x[5] >= ieq_01234
    x_0 <= x_1 <= x_2 <= x_3 <= x_4 <= x_5

    sage: (x[5]<=x[6]) >= ieq_01234
    x_0 <= x_1 <= x_2 <= x_3 <= x_4 <= x_5 <= x_6
    sage: (x[5]<=x[6]) <= ieq_01234
    x_5 <= x_6 <= x_0 <= x_1 <= x_2 <= x_3 <= x_4

TESTS:

See :trac:`12091` ::

    sage: p = SemidefiniteProgram()
    sage: b = p.new_variable()
    sage: b[0] <= b[1] <= 2
    x_0 <= x_1 <= 2
    sage: list(b[0] <= b[1] <= 2)
    [x_0, x_1, 2]
    sage: 1 >= b[1] >= 2*b[0]
    2*x_0 <= x_1 <= 1
    sage: b[2] >= b[1] >= 2*b[0]
    2*x_0 <= x_1 <= x_2
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
from cpython.object cimport *

cdef extern from "limits.h":
    long LONG_MAX

from sage.structure.parent cimport Parent
from sage.structure.element cimport ModuleElement, Element
from sage.misc.cachefunc import cached_function
from sage.matrix.matrix import Matrix

#*****************************************************************************
#
# Utility functions to test that something is a linear function / constraint
#
#*****************************************************************************

def is_LinearSDPFunction(x):
    """
    Test whether ``x`` is a linear function

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = SemidefiniteProgram()
        sage: x = p.new_variable()
        sage: from sage.numerical.linear_sdp_functions import is_LinearSDPFunction
        sage: is_LinearSDPFunction(x[0] - 2*x[2])
        True
        sage: is_LinearSDPFunction('a string')
        False
    """
    return isinstance(x, LinearSDPFunction)

def is_LinearSDPConstraint(x):
    """
    Test whether ``x`` is a linear constraint

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = SemidefiniteProgram()
        sage: x = p.new_variable()
        sage: ieq = (x[0] <= x[1])
        sage: from sage.numerical.linear_sdp_functions import is_LinearSDPConstraint
        sage: is_LinearSDPConstraint(ieq)
        True
        sage: is_LinearSDPConstraint('a string')
        False
    """
    return isinstance(x, LinearSDPConstraint)

#*****************************************************************************
#
# Factory functions for the parents to ensure uniqueness
#
#*****************************************************************************

@cached_function
def LinearSDPFunctionsParent(base_ring):
    """
    Return the parent for linear functions over ``base_ring``.

    The output is cached, so only a single parent is ever constructed
    for a given base ring.

    INPUT:

    - ``base_ring`` -- a ring. The coefficient ring for the linear
      funcitons.

    OUTPUT:

    The parent of the linear functions over ``base_ring``.

    EXAMPLES::

        sage: from sage.numerical.linear_sdp_functions import LinearSDPFunctionsParent
        sage: LinearSDPFunctionsParent(QQ)
        Linear SDP functions over Rational Field
    """
    return LinearSDPFunctionsParent_class(base_ring)

@cached_function
def LinearSDPConstraintsParent(linear_sdp_functions_parent):
   """
   Return the parent for linear functions over ``base_ring``.

   The output is cached, so only a single parent is ever constructed
   for a given base ring.

    INPUT:

    - ``linear_sdp_functions_parent`` -- a
      :class:`LinearSDPFunctionsParent_class`. The type of linear
      functions that the constraints are made out of.

    OUTPUT:

    The parent of the linear constraints with the given linear functions.

    EXAMPLES::

        sage: from sage.numerical.linear_sdp_functions import \
        ...       LinearSDPFunctionsParent, LinearSDPConstraintsParent
        sage: LF = LinearSDPFunctionsParent(QQ)
        sage: LinearSDPConstraintsParent(LF)
        Linear SDP constraints over Rational Field
   """
   return LinearSDPConstraintsParent_class(linear_sdp_functions_parent)



#*****************************************************************************
#
# Parent of linear functions
#
#*****************************************************************************

cdef class LinearSDPFunctionsParent_class(Parent):
    r"""
    The parent for all linear functions over a fixed base ring.

    .. warning::

        You should use :func:`LinearSDPFunctionsParent` to construct
        instances of this class.

    INPUT/OUTPUT:

    See :func:`LinearSDPFunctionsParent`

    EXAMPLES::

        sage: from sage.numerical.linear_sdp_functions import LinearSDPFunctionsParent_class
        sage: LinearSDPFunctionsParent_class
        <type 'sage.numerical.linear_sdp_functions.LinearSDPFunctionsParent_class'>
    """

    def __init__(self, base_ring):
        """
        The Python constructor

        TESTS::

            sage: from sage.numerical.linear_sdp_functions import LinearSDPFunctionsParent
            sage: LinearSDPFunctionsParent(RDF)
            Linear SDP functions over Real Double Field
        """
        from sage.categories.modules_with_basis import ModulesWithBasis
        Parent.__init__(self, base=base_ring, category=ModulesWithBasis(base_ring))
        self._multiplication_symbol = '*'

    def set_multiplication_symbol(self, symbol='*'):
        """
        Set the multiplication symbol when pretty-printing linear functions.

        INPUT:

        - ``symbol`` -- string, default: ``'*'``. The multiplication
          symbol to be used.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()
            sage: f = -1-2*x[0]-3*x[1]
            sage: LF = f.parent()
            sage: LF._get_multiplication_symbol()
            '*'
            sage: f
            -1 - 2*x_0 - 3*x_1
            sage: LF.set_multiplication_symbol(' ')
            sage: f
            -1 - 2 x_0 - 3 x_1
            sage: LF.set_multiplication_symbol()
            sage: f
            -1 - 2*x_0 - 3*x_1
        """
        self._multiplication_symbol = symbol

    cpdef _get_multiplication_symbol(self):
        """
        Return the multiplication symbol.

        OUTPUT:

        String. By default, this is ``'*'``.

        EXAMPLES::

            sage: LF = SemidefiniteProgram().linear_sdp_functions_parent()
            sage: LF._get_multiplication_symbol()
            '*'
        """
        return self._multiplication_symbol

    def _repr_(self):
        """
        Return as string representation

        EXAMPLES::

            sage: SemidefiniteProgram().linear_sdp_functions_parent()
            Linear SDP functions over Real Double Field
        """
        return 'Linear SDP functions over '+str(self.base_ring())

    cpdef _element_constructor_(self, x):
        """
        Construt a :class:`LinearSDPFunction` from ``x``.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: LF = p.linear_sdp_functions_parent()
            sage: LF._element_constructor_(123)
            123
            sage: p(123)    # indirect doctest
            123
            sage: type(_)
            <type 'sage.numerical.linear_sdp_functions.LinearSDPFunction'>

            sage: p_QQ = SemidefiniteProgram(solver='cvxopt')
            sage: LF_QQ = p_QQ.linear_sdp_functions_parent()
            sage: f = LF({-1:1/2, 2:3/4});  f
            0.5 + 0.75*x_2
            sage: LF(f) is f
            True
            sage: LF_QQ(f)
            0.5 + 0.75*x_2
            sage: LF_QQ(f) is f
            True
        """
        if is_LinearSDPFunction(x):
            if x.parent() is self:
                return x
            else:
                return LinearSDPFunction(self, (<LinearSDPFunction>x)._f)
        return LinearSDPFunction(self, x)

    cpdef _coerce_map_from_(self, R):
        """
        Allow coercion of scalars and matrices into linear functions.

        INPUT:

        - ``R`` -- a ring or a matrix.

        OUTPUT:

        Boolean. Whether there is a coercion map.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: parent = p.linear_sdp_functions_parent()
            sage: parent.coerce(int(2))
            2
            sage: parent._coerce_map_from_(int)
            True
            sage: m = matrix([4])
            sage: parent._coerce_map_from_(m)
            True
            sage: m = 4
            sage: parent._coerce_map_from_(m)
            False
        """
        if R in [int, float, long, complex]:
            return True
        from sage.rings.real_mpfr import mpfr_prec_min
        from sage.rings.all import RealField
        if RealField(mpfr_prec_min()).has_coerce_map_from(R):
            return True
        if isinstance(R,Matrix):
            return True
        return False

    def _an_element_(self):
        """
        Returns an element

        OUTPUT:

        A linear function.

        EXAMPLES::

            sage: p = SemidefiniteProgram().linear_sdp_functions_parent()
            sage: p._an_element_()
            5*x_2 + 7*x_5
            sage: p.an_element()   # indirect doctest
            5*x_2 + 7*x_5
        """
        return self._element_constructor_({2:5, 5:7})


#*****************************************************************************
#
# Elements of linear functions
#
#*****************************************************************************

cdef class LinearSDPFunction(ModuleElement):
    r"""
    An elementary algebra to represent symbolic linear functions.

    .. warning::

        You should never instantiate :class:`LinearSDPFunction`
        manually. Use the element constructor in the parent
        instead. For convenience, you can also call the
        :class:`SemidefiniteProgram` instance directly.

    EXAMPLES:

    For example, do this::

        sage: p = SemidefiniteProgram()
        sage: p({0 : 1, 3 : -8})
        x_0 - 8*x_3

    or this::

        sage: parent = p.linear_sdp_functions_parent()
        sage: parent({0 : 1, 3 : -8})
        x_0 - 8*x_3

    or even this::

        sage: parent = p.linear_sdp_functions_parent()
        sage: parent({0 : matrix([1]), 3 : -8})
        x_0 - 8*x_3

    instead of this::

        sage: from sage.numerical.linear_sdp_functions import LinearSDPFunction
        sage: LinearSDPFunction(p.linear_sdp_functions_parent(), {0 : 1, 3 : -8})
        x_0 - 8*x_3
    """

    def __init__(self, parent, f):
        r"""
        Constructor taking a dictionary or a numerical value as its argument.

        A linear function is represented as a dictionary. The
        values are the coefficient of the variable represented
        by the keys ( which are integers ). The key ``-1``
        corresponds to the constant term.

        EXAMPLES:

        With a dictionary::

            sage: p = SemidefiniteProgram()
            sage: p({0 : 1, 3 : -8})
            x_0 - 8*x_3

        With a dictionary and matrices::

            sage: p = SemidefiniteProgram()
            sage: p({0 : matrix([1]), 3 : -8})
            x_0 - 8*x_3

        Using the constructor with a numerical value::

            sage: p = SemidefiniteProgram()
            sage: p(25)
            25

        The matrices have to have the same dimensions though::

            sage: p = SemidefiniteProgram()
            sage: p({0 : matrix([[1,2],[3,4]]), 3 : -8})
            Traceback (most recent call last):
            ...
            Exception: ...

        The matrices also have to be squares::

            sage: p = SemidefiniteProgram()
            sage: p({0 : matrix([[1,2],[3,4],[5,6]])})
            Traceback (most recent call last):
            ...
            Exception: ...
        """
        ModuleElement.__init__(self, parent)
        R = self.base_ring()
        self._dim = -1

        if isinstance(f, dict):
            self._f = { }
            for key, value in f.iteritems():
                k = int(key)
                if not isinstance(value, Matrix):
                    if self._dim == -1:
                        self._dim = 1
                    if self._dim != 1:
                        raise Exception('The matrices have to have the same dimensions!')
                    self._f[k] = R(value)
                else:
                    if self._dim == -1:
                        self._dim = value.dimensions()[0]
                    if not value.is_square():
                        raise Exception('The matrices have to be squares!')
                    if self._dim != value.dimensions()[0]:
                        raise Exception('The matrices have to have the same dimensions!')
                    self._f[k] = value
        else:
            self._f = {-1: f}

    cpdef iteritems(self):
        """
        Iterate over the index, coefficient pairs

        OUTPUT:

        An iterator over the ``(key, coefficient)`` pairs. The keys
        are integers indexing the variables. The key ``-1``
        corresponds to the constant term.

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver = 'cvxopt')
            sage: x = p.new_variable()
            sage: f = 0.5 + 3/2*x[1] + 0.6*x[3]
            sage: for id, coeff in f.iteritems():
            ...      print 'id =', id, '  coeff =', coeff
            id = 0   coeff = 1.5
            id = 1   coeff = 0.6
            id = -1   coeff = 0.5
        """
        return self._f.iteritems()

    def dict(self):
        r"""
        Returns the dictionary corresponding to the Linear SDP Function.

        OUTPUT:

        The linear function is represented as a dictionary. The value
        are the coefficient of the variable represented by the keys (
        which are integers ). The key ``-1`` corresponds to the
        constant term.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: lf = p({0 : 1, 3 : -8})
            sage: lf.dict()
            {0: 1.0, 3: -8.0}
        """
        return dict(self._f)

    cpdef ModuleElement _add_(self, ModuleElement b):
        r"""
        Defining the + operator

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p({0 : matrix([[1,2],[3,4]])}) + p({2 : matrix([[2,2],[2,2]])})
            [1 2]       [2 2]
            [3 4]*x_0 + [2 2]*x_2
        """
        #YEAH ingo!!! YOU CAN DO IT!!!! laga thessar add adferdir, bua til matrix daemi.
        from sage.matrix.constructor import matrix
        from sage.rings.all import ZZ
        e = dict(self._f)
        P = self.parent()
        for (id,coeff) in b.dict().iteritems():
            print "coeff: " + str(coeff)
            print "f: " + str(self._f)
            if isinstance(coeff, Matrix):
                e[id] = self._f.get(id,matrix.zero(ZZ,coeff.dimensions()[0])) + coeff
                print "e[id]: " +  str(e)
            else:
                e[id] = self._f.get(id,0) + coeff
        return P(e)

    cpdef ModuleElement _neg_(self):
        r"""
        Defining the - operator (opposite).

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: - p({0 : 1, 3 : -8})
            -1*x_0 + 8*x_3
            sage: - p({0 : matrix([1]), 3 : -8})
            [-1]*x_0 + 8*x_3
        """
        P = self.parent()
        return P(dict([(id,-coeff) for (id, coeff) in self._f.iteritems()]))

    cpdef ModuleElement _sub_(self, ModuleElement b):
        r"""
        Defining the - operator (substraction).

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p({2 : 5, 3 : 2}) - 3
            -3 + 5*x_2 + 2*x_3
            sage: p({0 : 1, 3 : -8}) - p({2 : 5, 3 : 2}) - 16
            -16 + x_0 - 5*x_2 - 10*x_3
        """
        e = dict(self._f)
        for (id,coeff) in b.dict().iteritems():
            e[id] = self._f.get(id,0) - coeff
        P = self.parent()
        return P(e)

    cpdef ModuleElement _rmul_(self, RingElement b):
        r"""
        Left multiplication by scalars

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: p({2 : 5, 3 : 2}) * 3
            15*x_2 + 6*x_3
        """
        P = self.parent()
        return P(dict([(id,b*coeff) for (id, coeff) in self._f.iteritems()]))

    cpdef ModuleElement _lmul_(self, RingElement b):
        r"""
        Right multiplication by scalars

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: 3 * p({2 : 5, 3 : 2})
            15*x_2 + 6*x_3
        """
        return self._rmul_(b)

    cpdef _acted_upon_(self, x, bint self_on_left):
       """
       Act with scalars that do not have a natural coercion into
       ``self.base_ring()``

       EXAMPLES:

       Note that there is no coercion from ``RR`` to ``QQ``, but you
       can explicitly convert them::

           sage: 1/2 * 0.6
           0.300000000000000

       If there were a coercion, then 0.6 would have been coerced into
       6/10 and the result would have been rational. This is
       undesirable, which is why there cannot be a coercion on the
       level of the coefficient rings.

       By declaring an action of ``RR`` on the linear functions over
       ``QQ`` we make the following possible::

           sage: p = SemidefiniteProgram(solver='cvxopt')
           sage: p.base_ring()
           Real Double Field
           sage: x = p.new_variable()
           sage: x[0] * 0.6
           0.6*x_0
       """
       R = self.base_ring()
       try:
           x_R = R(x)
       except TypeError:
           return None
       return self._rmul_(x_R)

    def _coeff_formatter(self, coeff, constant_term=False):
        """
        Pretty-print multiplicative coefficient ``x``

        OUTPUT:

        String, including a trailing space if the coefficient is not
        one. Empty string otherwise.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: f = p(1);  type(f)
            <type 'sage.numerical.linear_sdp_functions.LinearSDPFunction'>
            sage: f._coeff_formatter(1)
            ''
            sage: f._coeff_formatter(1, constant_term=True)
            '1'
            sage: f._coeff_formatter(RDF(12.0))
            '12*'
            sage: f._coeff_formatter(RDF(12.3))
            '12.3*'

            sage: q = SemidefiniteProgram(solver='cvxopt')
            sage: f = p(1)
            sage: f._coeff_formatter(13/45)
            '13/45*'
        """
        R = self.base_ring()
        if coeff == R.one() and not constant_term:
            return ''
        try:
            from sage.rings.all import ZZ
            coeff = ZZ(coeff)    # print as integer if possible
        except TypeError:
            pass
        if constant_term:
            return str(coeff)
        else:
            return str(coeff) + self.parent()._get_multiplication_symbol()

    def _repr_(self):
        r"""
        Returns a string version of the linear function.

        EXAMPLE::

            sage: p = SemidefiniteProgram(solver='cvxopt')
            sage: p({-1: -15, 2 : -5.1, 3 : 2/3})
            -15 - 5.1*x_2 + 0.666666666667*x_3
        """
        cdef dict d = dict(self._f)
        cdef bint first = True
        t = ""

        if -1 in d:
            coeff = d.pop(-1)
            if coeff!=0:
                t = self._coeff_formatter(coeff, constant_term=True)
                first = False

        cdef list l = sorted(d.items())
        for id,coeff in l:
            sign = cmp(coeff,0)
            if sign == 0:
                continue
            if not first:
                if sign == -1:
                    t += ' - '
                if sign == +1:
                    t += ' + '
                t += self._coeff_formatter(abs(coeff)) + 'x_' + str(id)
            else:
                t += self._coeff_formatter(coeff) + 'x_' + str(id)
            first = False

        if first:
            return '0'
        else:
            return t

    cpdef is_zero(self):
        """
        Test whether ``self`` is zero.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()
            sage: (x[1] - x[1] + 0*x[2]).is_zero()
            True
        """
        for coeff in self._f.values():
            if not coeff.is_zero():
                return False
        return True

    cpdef equals(LinearSDPFunction left, LinearSDPFunction right):
        """
        Logically compare ``left`` and ``right``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()
            sage: (x[1] + 1).equals(3/3 + 1*x[1] + 0*x[2])
            True
        """
        return (left-right).is_zero()

    def __richcmp__(left, right, int op):
        """
        Override the rich comparison.

        The Sage framework sometimes expects that rich comparison
        results in a boolean value, but we want to return
        :class:`~sage.numerical.linear_sdp_functions.LinearSDPConstraint`
        objects.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()
            sage: x[0].__le__(x[1])    # indirect doctest
            x_0 <= x_1
        """
        return (<LinearSDPFunction>left)._richcmp(right, op)

    cdef _richcmp(left, right, int op):
        """
        Create a inequality or equality object.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: from sage.numerical.linear_sdp_functions import LinearSDPFunction
            sage: p({2 : 5, 3 : 2}) <= p({2 : 3, 9 : 2})
            5*x_2 + 2*x_3 <= 3*x_2 + 2*x_9

            sage: p({2 : 5, 3 : 2}) >= p({2 : 3, 9 : 2})
            3*x_2 + 2*x_9 <= 5*x_2 + 2*x_3

            sage: p({2 : 5, 3 : 2}) == p({2 : 3, 9 : 2})
            5*x_2 + 2*x_3 == 3*x_2 + 2*x_9

            sage: p({2 : 5, 3 : 2}) < p({2 : 3, 9 : 2})
            Traceback (most recent call last):
            ...
            ValueError: strict < is not allowed, use <= instead.

            sage: p({2 : 5, 3 : 2}) > p({2 : 3, 9 : 2})
            Traceback (most recent call last):
            ...
            ValueError: strict > is not allowed, use >= instead.

        TESTS::

            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.explain(10, p(1), operator.le)
            Coercion on left operand via
                Conversion map:
                  From: Integer Ring
                  To:   Linear SDP functions over Real Double Field
            Arithmetic performed after coercions.
            Result lives in Linear SDP functions over Real Double Field
            Linear SDP functions over Real Double Field

            sage: x = p.new_variable()
            sage: operator.le(10, x[0])
            10 <= x_0
            sage: x[0] <= 1
            x_0 <= 1
            sage: x[0] >= 1
            1 <= x_0
            sage: 1 <= x[0]
            1 <= x_0
            sage: 1 >= x[0]
            x_0 <= 1
       """
        LF = left.parent()
        LC = LinearSDPConstraintsParent(LF)
        equality = (op == Py_EQ)
        cdef LinearSDPConstraint  left_constraint = LC(left,  equality=equality)
        cdef LinearSDPConstraint right_constraint = LC(right, equality=equality)
        if op == Py_LT:
            raise ValueError("strict < is not allowed, use <= instead.")
        elif op == Py_EQ:
            return left_constraint._richcmp(right_constraint, op)
        elif op == Py_GT:
            raise ValueError("strict > is not allowed, use >= instead.")
        elif op == Py_LE:
            return left_constraint._richcmp(right_constraint, op)
        elif op == Py_NE:
            raise ValueError("inequality != is not allowed, use one of <=, ==, >=.")
        elif op == Py_GE:
            return left_constraint._richcmp(right_constraint, op)
        else:
            assert(False)   # unreachable

    def __hash__(self):
        r"""
        Return a hash.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: f = p({2 : 5, 3 : 2})
            sage: f.__hash__()   # random output
            103987752
            sage: d = {}
            sage: d[f] = 3
        """
        # see _cmp_c_impl() if you want to change the hash function
        return id(self) % LONG_MAX

    def __cmp__(left, right):
        """
        Part of the comparison framework.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: f = p({2 : 5, 3 : 2})
            sage: cmp(f, f)
            0
        """
        return (<Element>left)._cmp(right)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Implement comparison of two linear functions.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: f = p({2 : 5, 3 : 2})
            sage: cmp(f, f)
            0
            sage: abs(cmp(f, f+0))     # since we are comparing by id()
            1
            sage: abs(cmp(f, f+1))
            1
            sage: len(set([f, f]))
            1
            sage: len(set([f, f+0]))
            2
            sage: len(set([f, f+1]))
            2
        """
        # Note: if you want to implement smarter comparison, you also
        # need to change __hash__(). The comparison function must
        # satisfy cmp(x,y)==0 => hash(x)==hash(y)
        return cmp(id(left), id(right))

#*****************************************************************************
#
# Parent of linear constraints
#
#*****************************************************************************

cdef class LinearSDPConstraintsParent_class(Parent):
    """
    Parent for :class:`LinearSDPConstraint`

    .. warning::

        This class has no reason to be instanciated by the user, and
        is meant to be used by instances of
        :class:`SemidefiniteProgram`. Also, use the
        :func:`LinearSDPConstraintsParent` factory function.

    INPUT/OUTPUT:

        See :func:`LinearSDPFunctionsParent`

    EXAMPLES::

        sage: p = SemidefiniteProgram()
        sage: LC = p.linear_sdp_constraints_parent();  LC
        Linear SDP constraints over Real Double Field
        sage: from sage.numerical.linear_sdp_functions import LinearSDPConstraintsParent
        sage: LinearSDPConstraintsParent(p.linear_sdp_functions_parent()) is LC
        True
    """

    def __init__(self, linear_sdp_functions_parent):
        """
        The Python constructor

        INPUT/OUTPUT:

        See :func:`LinearSDPFunctionsParent`

        TESTS::

            sage: from sage.numerical.linear_sdp_functions import LinearSDPFunctionsParent
            sage: LF = LinearSDPFunctionsParent(RDF)
            sage: from sage.numerical.linear_sdp_functions import LinearSDPConstraintsParent
            sage: LinearSDPConstraintsParent(LF)
            Linear SDP constraints over Real Double Field
        """
        Parent.__init__(self)
        self._LF = linear_sdp_functions_parent

    def linear_sdp_functions_parent(self):
        """
        Return the parent for the linear functions

        EXAMPLES::

            sage: LC = SemidefiniteProgram().linear_sdp_constraints_parent()
            sage: LC.linear_sdp_functions_parent()
            Linear SDP functions over Real Double Field
        """
        return self._LF

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: SemidefiniteProgram().linear_sdp_constraints_parent()
            Linear SDP constraints over Real Double Field
        """
        return 'Linear SDP constraints over '+str(self.linear_sdp_functions_parent().base_ring())

    cpdef _element_constructor_(self, left, right=None, equality=False):
        """
        Construt a :class:`LinearSDPConstraint`.

        INPUT:

        - ``left`` -- a :class:`LinearSDPFunction`, or something that can
          be converted into one, a list/tuple of
          :class:`LinearSDPFunction`, or an existing
          :class:`LinearSDPConstraint`.

        - ``right`` -- a :class:`LinearSDPFunction` or ``None``
          (default).

        - ``equality`` -- boolean (default: ``True``). Whether to
          construct an equation or an inequality.

        OUTPUT:

        The :class:`LinearSDPConstraint` constructed from the input data.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: LC = p.linear_sdp_constraints_parent()
            sage: LC._element_constructor_(1, 2)
            1 <= 2

            sage: x = p.new_variable()
            sage: LC([x[0], x[1], x[2]])
            x_0 <= x_1 <= x_2

            sage: LC([x[0], x[1], x[2]], equality=True)
            x_0 == x_1 == x_2

            sage: type(_)
            <type 'sage.numerical.linear_sdp_functions.LinearSDPConstraint'>

        TESTS::

            sage: inequality = LC([x[0], 1/2*x[1], 3/4*x[2]]); inequality
            x_0 <= 0.5*x_1 <= 0.75*x_2
            sage: LC(inequality) is inequality
            True
            sage: p_QQ = SemidefiniteProgram(solver='cvxopt')
            sage: LC_QQ = p_QQ.linear_sdp_constraints_parent()
            sage: LC_QQ(inequality)
            x_0 <= 0.5*x_1 <= 0.75*x_2
        """
        if right is None and is_LinearSDPConstraint(left):
            if (left.parent() is self) and (left.is_equation() == equality):
                return left
            else:
                return LinearSDPConstraint(self, (<LinearSDPConstraint>left).constraints,
                                        equality=equality)
        if right is None:
            if isinstance(left, (list,tuple)):
                return LinearSDPConstraint(self, left, equality=equality)
            else:
                return LinearSDPConstraint(self, [left], equality=equality)
        else:
            return LinearSDPConstraint(self, [left, right], equality=equality)

    cpdef _coerce_map_from_(self, R):
        """
        Allow coercion of scalars into linear functions.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: parent = p.linear_sdp_constraints_parent()
            sage: parent.coerce(int(2))
            trivial constraint starting with 2
            sage: parent._coerce_map_from_(int)
            True
        """
        return self.linear_sdp_functions_parent().has_coerce_map_from(R)

    def _an_element_(self):
        """
        Returns an element

        EXAMPLES::

            sage: p = SemidefiniteProgram().linear_sdp_functions_parent()
            sage: p._an_element_()
            5*x_2 + 7*x_5
            sage: p.an_element()   # indirect doctest
            5*x_2 + 7*x_5
        """
        LF = self.linear_sdp_functions_parent()
        return self(0) <= LF.an_element()



#*****************************************************************************
#
# Elements of linear constraints
#
#*****************************************************************************

_chained_comparator_hack_search = None
_chained_comparator_hack_replace = None

cdef class LinearSDPConstraint(Element):
    """
    A class to represent formal Linear SDP Constraints.

    A Linear SDP Constraint being an inequality between
    two linear functions, this class lets the user
    write ``LinearSDPFunction1 <= LinearSDPFunction2``
    to define the corresponding constraint, which
    can potentially involve several layers of such
    inequalities (``(A <= B <= C``), or even equalities
    like ``A == B``.

    Trivial constraints (meaning that they have only one term and no
    relation) are also allowed. They are required for the coercion
    system to work.

    .. warning::

        This class has no reason to be instanciated by the user, and
        is meant to be used by instances of
        :class:`SemidefiniteProgram`.

    INPUT:

    - ``parent`` -- the parent, a :class:`LinearSDPConstraintsParent_class`

    - ``terms`` -- a list/tuple/iterable of two or more linear
      functions (or things that can be converted into linear
      functions).

    - ``equality`` -- boolean (default: ``False``). Whether the terms
      are the entries of a chained less-or-equal (``<=``) inequality
      or a chained equality.

    EXAMPLE::

        sage: p = SemidefiniteProgram()
        sage: b = p.new_variable()
        sage: b[2]+2*b[3] <= b[8]-5
        x_0 + 2*x_1 <= -5 + x_2
    """

    def __init__(self, parent, terms, equality=False):
        r"""
        Constructor for ``LinearSDPConstraint``

        INPUT:

        See :class:`LinearSDPConstraint`.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: b[2]+2*b[3] <= b[8]-5
            x_0 + 2*x_1 <= -5 + x_2
        """
        assert len(terms) > 0
        super(LinearSDPConstraint, self).__init__(parent)
        self.equality = equality
        LF = parent.linear_sdp_functions_parent()
        self.constraints = [ LF(term) for term in terms ]

    cpdef equals(LinearSDPConstraint left, LinearSDPConstraint right):
        """
        Compare ``left`` and ``right``.

        OUTPUT:

        Boolean. Whether all terms of ``left`` and ``right`` are
        equal. Note that this is stronger than mathematical
        equivalence of the relations.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()
            sage: (x[1] + 1 >= 2).equals(3/3 + 1*x[1] + 0*x[2] >= 8/4)
            True
            sage: (x[1] + 1 >= 2).equals(x[1] + 1-1 >= 1-1)
            False
        """
        if len(left.constraints) != len(right.constraints):
            return False
        if left.equality != right.equality:
            return False
        cdef LinearSDPFunction l, r
        for i in range(len(left.constraints)):
            l = <LinearSDPFunction>(left.constraints[i])
            r = <LinearSDPFunction>(right.constraints[i])
            if not l.equals(r):
                return False
        return True

    cdef LinearSDPConstraint _chained_comparator_hack_part1(LinearSDPConstraint left, LinearSDPConstraint right):
        """
        Evil hack to allow chained constraints

        Python translates ``x < y < z`` into:

        .. code-block:: python

             temp = x <= y      # calls x.__richcmp__(y)
             if temp:           # calls temp.__nonzero__()
               return y <= z    # calls y.__richcmp__(z)
             else:
               return temp

        but we would like ``x<=y<=z`` as output. The trick to make it
        work is to store ``y`` in the first call to ``__richcmp__()``
        and ``temp`` in the call to ``__nonzero__()``. Then we can
        replace ``y`` by ``x<=y`` in the second call to
        ``__richcmp__``.

        This function implements the first part of this hack, to be
        called from :meth:`__richcmp__`.
        """
        # print '__richcmp__', left, ' compared with', right
        global _chained_comparator_hack_search
        global _chained_comparator_hack_replace
        cdef LinearSDPConstraint search = _chained_comparator_hack_search
        cdef LinearSDPConstraint replace = _chained_comparator_hack_replace
        _chained_comparator_hack_search = right
        if replace is None:
            return left
        assert search is not None
        if search.equals(left):
            _chained_comparator_hack_replace = None
            return replace
        else:
            return left

    cdef _chained_comparator_hack_part2(self):
        """
        This function implements the first part of this hack, to be
        called from :meth:`__nonzero__`.
        """
        # print '__nonzero__', self.constraints
        global _chained_comparator_hack_replace
        _chained_comparator_hack_replace = self

    def is_equation(self):
        """
        Whether the constraint is a chained equation

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: (b[0] == b[1]).is_equation()
            True
            sage: (b[0] <= b[1]).is_equation()
            False
        """
        return self.equality

    def is_less_or_equal(self):
        """
        Whether the constraint is a chained less-or_equal inequality

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: (b[0] == b[1]).is_less_or_equal()
            False
            sage: (b[0] <= b[1]).is_less_or_equal()
            True
        """
        return not self.equality

    def is_trivial(self):
        """
        Test whether the constraint is trivial.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: LC = p.linear_sdp_constraints_parent()
            sage: ieq = LC(1,2);  ieq
            1 <= 2
            sage: ieq.is_trivial()
            False

            sage: ieq = LC(1);  ieq
            trivial constraint starting with 1
            sage: ieq.is_trivial()
            True
        """
        return len(self.constraints) < 2

    def __iter__(self):
        """
        Iterate over the terms of the chained (in)-equality

        OUTPUT:

        A generator yielding the individual terms of the constraint in
        left-to-right order.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: ieq = 1 <= b[0] <= b[2] <= 3 <= b[3];  ieq
            1 <= x_0 <= x_1 <= 3 <= x_2
            sage: list(ieq)
            [1, x_0, x_1, 3, x_2]
            sage: for term in ieq:
            ...       print term
            1
            x_0
            x_1
            3
            x_2
        """
        for term in self.constraints:
            yield term

    def equations(self):
        """
        Iterate over the unchained(!) equations

        OUTPUT:

        An iterator over pairs ``(lhs, rhs)`` such that the individual
        equations are ``lhs == rhs``.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: eqns = 1 == b[0] == b[2] == 3 == b[3];  eqns
            1 == x_0 == x_1 == 3 == x_2

            sage: for lhs, rhs in eqns.equations():
            ...       print str(lhs) + ' == ' + str(rhs)
            1 == x_0
            x_0 == x_1
            x_1 == 3
            3 == x_2
        """
        if not self.is_equation() or self.is_trivial():
            raise StopIteration
        term_iter = iter(self)
        lhs = term_iter.next()
        rhs = term_iter.next()
        while True:
            yield (lhs, rhs)
            lhs = rhs
            rhs = term_iter.next()

    def inequalities(self):
        """
        Iterate over the unchained(!) inequalities

        OUTPUT:

        An iterator over pairs ``(lhs, rhs)`` such that the individual
        equations are ``lhs <= rhs``.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: ieq = 1 <= b[0] <= b[2] <= 3 <= b[3]; ieq
            1 <= x_0 <= x_1 <= 3 <= x_2

            sage: for lhs, rhs in ieq.inequalities():
            ...       print str(lhs) + ' <= ' + str(rhs)
            1 <= x_0
            x_0 <= x_1
            x_1 <= 3
            3 <= x_2
        """
        if not self.is_less_or_equal() or self.is_trivial():
            raise StopIteration
        term_iter = iter(self)
        lhs = term_iter.next()
        rhs = term_iter.next()
        while True:
            yield (lhs, rhs)
            lhs = rhs
            rhs = term_iter.next()

    def _repr_(self):
        r"""
        Returns a string representation of the constraint.

        OUTPUT:

        String.

        EXAMPLE::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: b[3] <= b[8] + 9
            x_0 <= 9 + x_1

            sage: LC = p.linear_sdp_constraints_parent()
            sage: LC(b[3], b[8] + 9)
            x_0 <= 9 + x_1
            sage: LC(b[3])
            trivial constraint starting with x_0
        """
        comparator = ( ' == ' if self.equality else ' <= ' )
        result = comparator.join(map(str, self))
        if self.is_trivial():
            return 'trivial constraint starting with '+result
        return result

    def __nonzero__(self):
        """
        Part of the hack to allow chained (in)equalities

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: ieq = (b[3] <= b[8] + 9)
            sage: ieq <= ieq <= ieq
            x_0 <= 9 + x_1 <= x_0 <= 9 + x_1 <= x_0 <= 9 + x_1
        """
        self._chained_comparator_hack_part2()
        return True

    def __richcmp__(left, right, int op):
        """
        Override the rich comparison.

        The Sage framework sometimes expects that rich comparison
        results in a boolean value, but we want to return
        :class:`~sage.numerical.linear_sdp_functions.LinearSDPConstraint`
        objects.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: b[0] <= b[1] <= b[2] <= b[3]
            x_0 <= x_1 <= x_2 <= x_3
            sage: b[0] <= 1 <= b[1] <= 2 <= b[2] <= 3
            x_0 <= 1 <= x_1 <= 2 <= x_2 <= 3
        """
        return (<LinearSDPConstraint>left)._richcmp(right, op)

    cdef _richcmp(py_left, py_right, int op):
        """
        Chain (in)equalities
        """
        #  print 'richcmp', py_left, ', ', py_right
        LC = py_left.parent()
        if not is_LinearSDPConstraint(py_right):
            py_right = LC(py_right, equality=py_left.is_equation())
        elif py_right.parent() is not LC:
            py_right = LC(py_right.constraints, equality=py_left.is_equation())
        cdef LinearSDPConstraint right = <LinearSDPConstraint>py_right
        cdef LinearSDPConstraint left = py_left._chained_comparator_hack_part1(right)
        if op == Py_LT:
            raise ValueError("strict < is not allowed, use <= instead.")
        elif op == Py_EQ:
            if not (left.is_equation() and right.is_equation()):
                raise ValueError("can only chain together equations")
            return LC(left.constraints + right.constraints, equality=True)
        elif op == Py_GT:
            raise ValueError("strict > is not allowed, use >= instead.")
        elif op == Py_LE:
            if not (left.is_less_or_equal() and right.is_less_or_equal()):
                raise ValueError("can only chain together inequalities")
            return LC(left.constraints + right.constraints)
        elif op == Py_NE:
            raise ValueError("inequality != is not allowed, use one of <=, ==, >=.")
        elif op == Py_GE:
            if not (left.is_less_or_equal() and right.is_less_or_equal()):
                raise ValueError("can only chain together inequalities")
            return LC(right.constraints + left.constraints)
        else:
            assert(False)   # unreachable
