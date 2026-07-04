import random

from dataclasses import dataclass, field


from sage.all import (
    ascii_art,
    block_matrix,
    floor,
    latex,
    lcm,
    matrix,
    Parent,
    QQ,
    RingElement,
    solve_diophantine,
    var,
    ZZ,
)
from sage.categories.rings import Rings
from .normalizer import standardize_sol_matrix
from .space_groups import SpaceGroup_gap
from .meta import Option
from .congruences import solve_qz_integral, solve_qz_rational, Solution_QZ_rational


# TODO:
#  - somehow this should be in the category of Rings, but it doesn't appear there
class _Q_modZ(Parent):

    def __init__(self):
        Parent.__init__(self, Rings)

    def _element_constructor_(self, x):
        return Q_modZ_Element(self, x)

    def _repr_(self):
        return "Q/Z"


# TODO:
#  - change it to Z-module
#  - add sage coercion
#  - make this shit work with matrix(QmodZ, [[...]])
class Q_modZ_Element(RingElement):

    def __init__(self, parent, r) -> None:
        RingElement.__init__(self, parent)
        self.r = QQ(r)
        self._normalize()

    def _normalize(self):
        self.r = self.r - floor(self.r)

    def __add__(self, a):
        if isinstance(a, Q_modZ_Element):
            return Q_modZ_Element(self.parent(), self.r + a.r)
        else:
            return Q_modZ_Element(self.parent(), self.r + QQ(a))

    def __sub__(self, a):
        if isinstance(a, Q_modZ_Element):
            return Q_modZ_Element(self.parent(), self.r - a.r)
        else:
            return Q_modZ_Element(self.parent(), self.r - QQ(a))

    def __radd__(self, a):
        return self.__add__(a)

    def __rsub__(self, a):
        if isinstance(a, Q_modZ_Element):
            return Q_modZ_Element(self.parent(), a.r - self.r)
        else:
            return Q_modZ_Element(self.parent(), QQ(a) - self.r)

    def __eq__(self, value: object, /) -> bool:
        if isinstance(value, Q_modZ_Element):
            return self.r == value.r
        else:
            return self.r == value

    def __mul__(self, a):
        if isinstance(a, Q_modZ_Element):
            return Q_modZ_Element(self.parent(), self.r * a.r)
        else:
            return Q_modZ_Element(self.parent(), self.r * QQ(a))

    def inv(self):
        res = 1 / self
        if res == 0:
            raise ValueError('inverse is zero')
        return res

    def _div_(self, a):
        if isinstance(a, Q_modZ_Element):
            return Q_modZ_Element(self.parent(), self.r / a.r)
        else:
            return Q_modZ_Element(self.parent(), self.r / QQ(a))

    def _repr_(self):
        return repr(self.r)


QmodZ = _Q_modZ()


@dataclass
class SolutionSimpleMatrix:
    except_values: list
    extra_eqs: list = field(default_factory=list)

    def add_equation(self, eq):
        self.extra_eqs.append(eq)
        return self

    def solution_list(self):
        # TODO: add combination with extra_eqs
        return sorted(self.except_values)


def _solve_simple_mat2(A) -> Option[SolutionSimpleMatrix]:
    """Given matrix of integral indeterminates A, find when it forms a simple virtual endomorphism
    of Z^n, using Nekrashevych theorem.


    Notes
    -----
    A matrix A^{-1} forms a simple virtual endomorphism if and only its characteristic polynomial isn't divisible
    by a monic polynomial with integer coefficients. In other words, A^{-1} is simple if and only if A^{-1} has no
    eigenvalue which is algebraic integer.

    For a 2x2 matrix it simplifies to the conditions:

    1. |det(A)| != 1
    2. det(A) + tr(A) != -1
    2. det(A) - tr(A) != -1

    if tr(A) == 0, then we can further simplify condition |det(A)| != 1 by having |bc| != 1 for integral b, c entities
    of the matrix A. From which follow four simpler equations
    [
        b !=  1 && c !=  1,
        b !=  1 && c != -1
        b != -1 && c !=  1
        b != -1 && c != -1,
    ]
    """

    # if A.trace() == 0:
    #     cond1 = (A[0, 1] - 1).simplify_rational(), (A[1, 0] - 1).simplify_rational()
    #     cond2 = (A[0, 1] - 1).simplify_rational(), (A[1, 0] + 1).simplify_rational()
    #     cond3 = (A[0, 1] + 1).simplify_rational(), (A[1, 0] - 1).simplify_rational()
    #     cond4 = (A[0, 1] + 1).simplify_rational(), (A[1, 0] + 1).simplify_rational()

    #     sols = []
    #     for cond in [cond1, cond2, cond3, cond4]:
    #         sol = solve(cond)
    #         if isinstance(sol, dict):
    #             sol = [sol,]

    cond1 = (A.det() + A.trace()).simplify_rational()
    cond2 = (A.det() - A.trace()).simplify_rational()

    if len(A.det().variables()) > 2:
        return Option.error('cant solve for more than 2 variables')

    det_res1 = solve_diophantine(A.det() - 1, A.det().variables(), solution_dict=True)
    det_res2 = solve_diophantine(A.det() + 1, A.det().variables(), solution_dict=True)

    if isinstance(det_res1, dict):
        det_res1 = [
            det_res1,
        ]

    if isinstance(det_res2, dict):
        det_res2 = [
            det_res2,
        ]

    det_res1 = [tuple([val for _, val in sorted(sol.items())]) for sol in det_res1]
    det_res2 = [tuple([val for _, val in sorted(sol.items())]) for sol in det_res2]

    eig_res1 = solve_diophantine(cond1 + 1, cond1.variables(), solution_dict=True)
    eig_res2 = solve_diophantine(cond2 + 1, cond2.variables(), solution_dict=True)

    if isinstance(eig_res1, dict):
        eig_res1 = [
            eig_res1,
        ]

    if isinstance(eig_res2, dict):
        eig_res2 = [
            eig_res2,
        ]

    eig_res1 = [tuple([val for _, val in sorted(sol.items())]) for sol in eig_res1]
    eig_res2 = [tuple([val for _, val in sorted(sol.items())]) for sol in eig_res2]

    return Option.success(
        SolutionSimpleMatrix(list(set(det_res1 + det_res2 + eig_res1 + eig_res2)))
    )


def _subminor(A: matrix, i):
    n = A.dimensions()[0]
    return matrix([[A[j, k] for k in range(n) if k != i] for j in range(n) if j != i])


def factorize(A: matrix):
    n = A.dimensions()[0]
    for i, row in enumerate(A):
        if list(row).count(0) == (n-1) and row[i] != 0:
            return _subminor(A, i), row[i]

    for i, col in enumerate(A.T):
        if list(col).count(0) == (n-1) and col[i] != 0:
            return _subminor(A, i), col[i]
    return None


def tau(A: matrix):
    row1, row2, row3 = A
    return _tau(*row1, *row2, *row3)


def _tau(y0, y1, y2, y3, y4, y5, y6, y7, y8):
    return y0 * y4 + y0 * y8 + y4 * y8 - y1 * y3 - y2 * y6 - y5 * y7


def _solve_simple_mat3(A) -> Option[SolutionSimpleMatrix]:
    subA = factorize(A)
    # case of permutational matrix
    if subA is None:
        if not A.trace() == 0 or not tau(A) == 0:
            raise ValueError("something went wrong, this shouldn't be like that")
        return Option.success(
            SolutionSimpleMatrix([]).add_equation(A.det() == 1).add_equation(A.det() == -1)
        )
    else:
        subA, varx = subA
        res = _solve_simple_mat2(subA)
        if res.failed:
            return res
        else:
            res.result.add_equation(varx == 1)
            return res


def solve_simple_mat(A) -> Option[SolutionSimpleMatrix]:
    assert len(list(A)) == len(list(A[0]))
    if len(list(A)) == 2:
        return _solve_simple_mat2(A)
    elif len(list(A)) == 3:
        return _solve_simple_mat3(A)
    else:
        return Option.error('not implemented for higher dimensions.')


@dataclass
class SC_DegreesSolution:
    matrices: list = field(default_factory=list)
    _unique_matrices: dict = field(default_factory=dict)

    def add_matrix(self, m: matrix):
        h = str(m)
        if h not in self._unique_matrices:
            self._unique_matrices[h] = m
            self.matrices.append(m)

    def unique_determinants(self):
        found = set()
        for M in self.matrices:
            p = M.det()
            p = self._standardize_poly(p)
            if str(p) not in found:
                yield p
                found.add(str(p))

    @staticmethod
    def _standardize_poly(p):
        vs = p.variables()
        tmp_vars = [var(f't{i}') for i in range(len(vs))]
        tmp_p = p.subs({v: tmp for v, tmp in zip(vs, tmp_vars)})
        resvars = [var(f'y{i}') for i in range(len(vs))]
        return tmp_p.subs({tmp: new_v for tmp, new_v in zip(tmp_vars, resvars)})


class SR_Degrees:
    def __init__(self, group_index, method="ascii", verbose=2, dim=2):
        self.congruence_print_thresh = 10
        self.disp = method
        self.group_index = group_index
        self.verbose = verbose
        if method == "latex":
            self._display = latex
            self.title = "\\subsection{Group %d}" % group_index
            self.section = "\\subsubsection{%s}"
            self.neq = "\\neq"
            self.pm = "\\pm"
            self.pref = "$$"
            self.in_sym = "\\in"
            self.notin_sym = "\\notin"
            self.z = "\\mathbb{Z}"
        elif method == "markdown":
            self._display = latex
            self.title = f"## Group {group_index}"
            self.section = "### %s"
            self.pref = "$$"
            self.pm = "\\pm"
            self.neq = "\\neq"
            self.in_sym = "\\in"
            self.notin_sym = "\\notin"
            self.z = "\\mathbb{Z}"
        else:
            self._display = ascii_art
            self.title = f"=========================== Group {group_index} =================================="
            self.section = "\n----------------%s-------------------------\n"
            self.pref = "\n"
            self.pm = '+-'
            self.neq = '!='
            self.in_sym = "  ∈  "
            self.notin_sym = "  /∈  "
            self.z = "Z"

        self.G = SpaceGroup_gap.from_gap_cryst(group_index, dim=dim, change_basis=True)

    def display(self, *args, use_pref=True, sep='') -> str:
        bod = self._display("") if self.disp == "ascii" else ""
        first = True
        for el in args:
            if not first and sep:
                if self.disp in ['latex', 'markdown']:
                    bod += sep
                else:
                    bod += self._display(sep)
            first = False

            if isinstance(el, str) and self.disp in ["latex", "markdown"]:
                bod += el
            else:
                bod += self._display(el)

        if use_pref:
            return self.pref + str(bod) + self.pref
        else:
            return str(bod)

    def print(self, *args, **kwargs):
        if self.verbose > 1:
            print(*args, **kwargs)

    def header(self):
        if self.verbose == 0:
            return
        print(self.title)
        print("Generators of group:")
        print(self.display([matrix(QQ, el) for el in self.G.gap_G.GeneratorsOfGroup()]))  # type: ignore
        print("SNoT")
        print(self.display(self.G.snot))

    def texdoc_header(self):
        if self.verbose == 0:
            return

        if self.disp == "latex":
            print("\\documentclass[12pt]{article}")
            print("\\usepackage{a4wide}")
            print("\\usepackage{amsmath,amssymb,amsthm}")
            print("\\title{Planar Groups}")
            print("\\begin{document}")
            print("\\maketitle")
            print("\\tableofcontents")
            print("\\section{Planar Groups}")

    def texdoc_ending(self):
        if self.verbose == 0:
            return

        if self.disp == "latex":
            print("\\end{document}")

    def construct_congruences_v2(self, A_inv, A):
        """Construct congruences for the subgroup check.

        The congruences have the following form:
            (A, t)(B, alpha(b))(A^{-1}, -A^{-1}t) = (ABA^{-1}, alpha(ABA^{-1})) mod L

        The action on the lattice L coincides with a linear operator A,
        so the index equals to det(A):
            (A, t)(E, e)(A, t)^{-1} = (E, At)

        Calculate the conjugation for every element in the SNoT:
            (A, t)(B, alpha(B))(A^{-1}, -A^{-1}t) =
            = (A, t)(BA^{-1}, -BA^{-1}t + alpha(B)) =

            = (ABA^{-1}, -ABA^{-1}t + A alpha(B) + t)
            = (ABA^{-1}, (E - ABA^{-1})t + A alpha(B)) = (ABA^{-1}, alpha(ABA^{-1})) mod L

            Hence, the result system can be written

            (E - ABA^{-1})t + A alpha(B) = alpha(ABA^{-1})

            or

            A alpha(B) - alpha(ABA^{-1}) = (ABA^{-1} - E)t  mod L

        We can obviously abuse multiple appearence of ABA^{-1}.
        """
        G = self.G
        P = self.G.point_group_elements()
        trans_a = list(var(' '.join(f"a{i}" for i in range(G.dim))))
        x = matrix(trans_a).T
        E = matrix.identity(QQ, G.dim)

        conds = []
        base_variables = list(A_inv.variables()) + list(trans_a)

        if len(P) > 5:
            self.print(f'Note that the point group consists of {len(P)} elements. Omitting most matrices..')
            indexes_to_print = list(range(len(P)))
            random.shuffle(indexes_to_print)
            indexes_to_print = indexes_to_print[:5]
        else:
            indexes_to_print = list(range(len(P)))
        for i, g in enumerate(P):
            # there is no point of checking identity element, since it always satisfies cocycle compatibility
            # a^{-1} (E, 0) a = (E, 0)
            if g.is_one():
                continue
            conj = (A * g * A_inv).simplify_rational()
            assert G.in_alpha(conj), conj

            # condition = A * G.alpha(g) - G.alpha(conj)
            # condition = condition.simplify_rational()

            # sym_res = (conj - E) * x

            alpha_conj = (E - conj) * x + (A * G.alpha(g))
            real_alpha = G.alpha(conj)

            for ((ll,), (rr,)) in zip(alpha_conj, real_alpha):
                conds.append(ll - rr)

            if i in indexes_to_print:
                self.print_congruence(g, conj, alpha_conj, real_alpha)

        return conds, base_variables

    def construct_congruences_v3(self, A_inv, A):
        G = self.G
        P = self.G.point_group_elements()

        M = []
        b_var = []
        b_num = []
        E = matrix.identity(QQ, G.dim)
        for i, g in enumerate(P):
            if g.is_one():
                continue

            conj = (A * g * A_inv).simplify_rational()
            assert G.in_alpha(conj), conj
            M.extend((conj - E).rows())
            b_var.extend((A * G.alpha(g)).rows())
            b_num.extend(G.alpha(conj).rows())

        M = matrix(ZZ, M)

        self.print('matrix M of (sigma(P) - E):')
        self.print(self.display(M))

        b_var = matrix(b_var)
        b_num = matrix(QQ, b_num)

        S, U, _ = M.smith_form()
        r = M.rank()

        ui = U.matrix_from_rows(range(r, U.nrows()))

        res_left = ui * b_var
        res_right = ui * b_num
        return res_left, res_right

    def solve_congruences_v5(self, left, right, base_variables=None):
        if base_variables is None:
            base_variables = left.variables()

        M = matrix(QQ, [[cond.coefficient(v) for v in base_variables] for cond in left.list()])
        return solve_qz_integral(M, right, base_variables)

    def print_system(self, *eqs):
        if self.disp in ["latex", "markdown"]:
            if self.disp == 'markdown':
                self.print('$$')
            self.print(r"\begin{align*}")
            for eq in set(eqs):
                # print(eq)
                self.print(self.display(eq, '&=', 0, r'\quad', r'\mod' r'\mathbb{Z}\\', use_pref=False))
            self.print(r"\end{align*}")
            if self.disp == 'markdown':
                self.print('$$')
            self.print()
        else:
            for eq in set(eqs):
                self.print(self.display(eq == 0))

    def print_congruence(self, g, conj, alpha_conj, real_alpha):
        G = self.G

        def construct_el(x, y):
            return block_matrix([[x, y], [0, 1]])

        for j in range(self.G.dim):
            alpha_conj[j, 0] = alpha_conj[j, 0].simplify_rational()
        if self.disp in ["latex", "markdown"]:
            if self.disp == 'latex':
                self.print(r'\begin{footnotesize}')
            self.print(
                self.pref
                + "(A, a)"
                + self.display(construct_el(g, G.alpha(g)), use_pref=False)
                + "(A^{-1}, -A^{-1}a) = "
            )
            self.print(self.display(construct_el(conj, alpha_conj), use_pref=False) + "=")
            self.print(self.display(construct_el(conj, real_alpha), use_pref=False) + self.pref)
            if self.disp == 'latex':
                self.print(r'\end{footnotesize}')
        else:
            self.print(
                ascii_art("\na\n\n")
                + ascii_art(" ")
                + ascii_art(construct_el(g, G.alpha(g)))
                + ascii_art("\na_inv\n\n")
                + ascii_art(" ")
                + ascii_art("\n=\n\n")
                + ascii_art(" ")
                + ascii_art(construct_el(conj, alpha_conj))
                + ascii_art(" ")
                + ascii_art("\n=\n\n")
                + ascii_art(" ")
                + ascii_art(construct_el(conj, real_alpha))
            )
            self.print()

    def lattice_compat(self, A) -> Option[matrix]:
        left = A
        right = matrix(QQ, [0 for _ in range(len(left.list()))]).T
        res = self.solve_congruences_v5(left, right)
        if res.failed:
            return Option.error(res.error_msg)

        ans = res.result
        A_new = A.subs({var: exp for var, exp in zip(list(ans.base_variables), list(ans.expressions))})
        return Option.success(standardize_sol_matrix(A_new, new_var='z'))

    def cocycle_compat(self, A) -> bool:
        A_inv = A.inverse().simplify_rational()
        eqs, base_vars = self.construct_congruences_v2(A_inv, A)
        ans = self.solve_congruences_v4(eqs, base_vars, list())
        return not ans.failed

    def cocycle_compat_v2(self, A) -> Option[matrix]:
        A_inv = A.inverse().simplify_rational()
        left, right = self.construct_congruences_v3(A_inv, A)
        self.print('congruence system:')
        self.print(self.display(left, right, use_pref=True, sep='='))
        res = self.solve_congruences_v5(left, right, base_variables=A.variables())
        if res.failed:
            return Option.error(res.error_msg)

        ans = res.result
        A_new = A.subs({var: exp for var, exp in zip(list(ans.base_variables), list(ans.expressions))})
        return Option.success(standardize_sol_matrix(A_new, new_var='y'))

    def print_large_mtx(self, mtx):
        rows = mtx.nrows()
        cols = mtx.ncols()
        self.print(r'\begin{pmatrix}')
        for i in range(self.congruence_print_thresh // 2):
            self.print(self.display(*mtx[i], sep=' & ', use_pref=False), r'\\')
        self.print(' & '.join([r'\cdots'] * cols), r'\\')
        for i in range(rows - self.congruence_print_thresh // 2, rows):
            self.print(self.display(*mtx[i], sep=' & ', use_pref=False), r'\\')
        self.print(r'\end{pmatrix}')

    def solve_congruences_v4(self, conds, base_variables, variables) -> Option[Solution_QZ_rational]:
        subs = {str(el): 0 for el in base_variables}

        # ugly sage method resolution
        def is_scalar(cond):
            if hasattr(cond, 'is_constant'):
                return cond.is_constant()
            return cond.substitute(**subs) == cond

        # basic clearing of trivial equations
        cleared_conds = []
        for rr in conds:
            sc = is_scalar(rr)
            # ignore trivial
            if sc and rr == 0:
                continue
            # check on conditions of kind 1/2 = 0  mod Z
            if sc and rr != 0:
                self.print(f'no solutions: {rr} != 0')
                return Option.error(f'no solutions: {rr} != 0')

            # check on conditions of kind 2 = 4  mod Z
            if sc and rr.is_integer():   # hope sage finds method `is_integer`
                continue
            cleared_conds.append(rr)

        # m * f = 0 in G-module Z^n for m=|G| => we have m as a universal common denominator
        m_bound = len(self.G.point_group_elements())

        M = matrix(QQ, [[cond.coefficient(v) for v in base_variables] for cond in cleared_conds])
        # f = matrix(base_variables).T

        # N_conds = len(cleared_conds)
        # N_vars = len(base_variables)
        # f_prox = matrix([var(f"y{i}") for i in range(len(base_variables))])

        # compute -f(0, ... 0) to get equations' free part b in Ax = b mod Z
        b = [-cond({v: 0 for v in base_variables + variables}) for cond in cleared_conds]
        b = matrix(QQ, b).T

        x = matrix(base_variables).T
        # FIXME: works only for latex
        if len(cleared_conds) > self.congruence_print_thresh:
            self.print('$$')
            self.print_large_mtx(M)
            self.print(self.display('*', x, '=', use_pref=False))
            self.print_large_mtx(b)
            self.print(self.display(r'\quad', r'\mod \mathbb{Z}', use_pref=False))
            self.print('$$')
        else:
            self.print(self.display(M, '*', x, '=', b, r'\quad', r'\mod \mathbb{Z}'))

        m = lcm(M.denominator(), b.denominator())
        assert (m_bound % m) == 0
        return solve_qz_rational(M, b)

    def generate_texdoc(self):
        self.texdoc_header()
        self.sc_degrees()
        self.texdoc_ending()

    def sc_degrees(self) -> SC_DegreesSolution:
        self.header()

        G = self.G
        norms = self.G.point_group_normalizer(verbose=False, ignore_trivial=True)

        self.print(self.display(norms))
        self.print(self.section % "Dilation")

        if G.is_symmorphic():
            self.print(
                f"Group {self.group_index} is a semi-direct product, therefore the dilation part"
                "is trivial and only consists of integral vectors."
            )

        if self.disp == "latex" and norms:
            self.print("\\begin{enumerate}")
            pref2 = "\\item"
        elif self.disp == 'markdown':
            pref2 = "- "
        else:
            pref2 = "..."

        sc_degrees = SC_DegreesSolution()

        for A in norms:
            self.print(pref2 + " testing A (should have integral entities):")
            self.print(self.pref + "A = \n" + self.display(A, use_pref=False) + self.pref)

            B = self.lattice_compat(A).result

            # A_inv = A.inverse().simplify_rational()
            if G.is_symmorphic():
                sc_degrees.add_matrix(standardize_sol_matrix(B, 'y'))
            else:
                new_a_opt = self.cocycle_compat_v2(B)
                if not new_a_opt.failed:
                    sc_degrees.add_matrix(new_a_opt.result)
                else:
                    self.print('not compatible with cocycle.')
                    continue

            # self.print("Simplicity")
            # self.print(self.pref + "A^{-1} = \n" + self.display(A_inv.simplify_rational(), use_pref=False) + self.pref)
            # self.print("\neigenvalues:")
            # self.print(self.display([el[0] for el in A_inv.charpoly().roots()]))
            # self.print("charpoly:")
            # chp = A_inv.charpoly()(SR("x"))
            # chp = factor(chp)
            # self.print(self.display(chp))
            # self.print("\nindex of subgroup:")
            # self.print(self.pref + "[G : H] = \n" + self.display(A.det(), use_pref=False) + self.pref)

        if self.disp == "latex" and norms:
            self.print("\\end{enumerate}")

        return sc_degrees

    def sr_degrees(self):
        """Compute conditions on self replicating degrees of a crystallographic group.

        Notes
        -----

        Given a matrix

            A  =   [
                [x0, x1],
                [x2, x3]
            ]

        that generates a surjective virtual endomorphism phi : H -> G by a
        conjugation phi(g) = (A, t)^{-1} g (A, t).

        From this follows that A has integral coefficients and |det(A)| > 1.
        Then, phi is simple if and only if
            1) det(A) - tr(A) != -1
            2) det(A) + tr(A) != -1     (sic!)

        which is equivalent that eigenvalues of A are not +-1.

        A proof comes from the h(1) != 0 and h(-1) != 0, where h(x) is a charpoly of A.
        """
        sc_degrees = self.sc_degrees()

        self.print(self.section % "Self-replicating degrees.")

        if self.disp == "latex" and sc_degrees:
            self.print("\\begin{enumerate}")
            pref2 = "\\item"
        elif self.disp == 'markdown':
            pref2 = '--- '
        else:
            pref2 = "..."

        eig_table = [["matrix", "det", "det(A) + tr(A) != -1", "det(A) - tr(A) != -1"]]
        res_table = [["matrix", "det", "except"]]
        for A in sc_degrees.matrices:
            self.print(pref2)
            self.print(self.display("A = ", A))
            self.print("Determinant:")
            self.print(self.display(A.det(), self.in_sym, self.z))

            self.print("\nSelf-replicating degrees:")
            self.print(self.display(A.det(), self.neq, self.pm, 1))

            # condition that eigenvalues isnt equal to +-1
            cond1 = (A.det() + A.trace()).simplify_rational()
            cond2 = (A.det() - A.trace()).simplify_rational()
            self.print(self.display(cond1, self.neq, -1))
            self.print(self.display(cond2, self.neq, -1))

            res = solve_simple_mat(A).result

            eig_table.append([A, A.det(), cond1 + 1, cond2 + 1])
            res_table.append([A, A.det(), res.solution_list()])
        if self.disp == "latex" and sc_degrees:
            self.print("\\end{enumerate}")
            self.print("\\newpage")

        return eig_table, res_table
