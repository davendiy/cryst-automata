from sage.all import (
    QQ,
    SR,
    ascii_art,
    block_matrix,
    factor,
    floor,
    latex,
    matrix,
    Parent,
    RingElement,
    solve,
    solve_diophantine,
    table,
    var,
)
from sage.categories.rings import Rings
from sage.modules.free_module_integer import IntegerLattice
from .space_groups import SpaceGroup_gap


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
        return 1 / self

    def _div_(self, a):
        if isinstance(a, Q_modZ_Element):
            return Q_modZ_Element(self.parent(), self.r / a.r)
        else:
            return Q_modZ_Element(self.parent(), self.r / QQ(a))

    def _repr_(self):
        return repr(self.r)


QmodZ = _Q_modZ()


# ugly gauss from internet
def gauss_elim(mat):
    N = len(mat)
    M = len(mat[0])
    for k in range(M):
        i_max = max(range(k, N), key=lambda i: abs(mat[i][k].r))
        if mat[i_max][k] == 0:
            continue
        mat[k], mat[i_max] = mat[i_max], mat[k]

        for i in range(k + 1, N):
            f = mat[i][k] / mat[k][k]

            for j in range(k + 1, M):
                mat[i][j] = mat[i][j] - f * mat[k][j]

            mat[i][k] = 0


class SC_DegreesSolution:

    def __init__(self):
        self._matrices = {}
        self._body = {}

    def add_solution(self, mat, inv_mat, conditions, solutions):
        self._matrices[str(mat)] = mat
        self._body[str(mat)] = inv_mat, conditions, solutions

    def __iter__(self):
        for mat in self._matrices:
            yield self._matrices[mat], *self._body[mat]


class SR_DegreesSolution:

    def __init__(self):
        pass

    def enumerate(self):
        pass

    def _latex_(self):
        pass


def _solve_simple_mat2(A):
    """Given matrix of indeterminates A, find when it forms a simple virtual endomorphism
    of Z^n, using Nekrashevych theorem.


    Notes
    -----
    A matrix A forms a simple virtual endomorphism if and only its characteristic polynomial isn't divisible
    by a monic polynomial with integer coefficients. In other words, A is simple if and only if A has no
    eigenvalue which is algebraic integer.

    For a 2x2 matrix it simplifies to the conditions:

    1. |det(A)| != 1
    2. det(A) + tr(A) != -1
    2. det(A) - tr(A) != -1
    """
    cond1 = (A.det() + A.trace()).simplify_rational()
    cond2 = (A.det() - A.trace()).simplify_rational()

    det_res1 = solve_diophantine(A.det() - 1, A.det().variables(), solution_dict=True)
    det_res2 = solve_diophantine(A.det() + 1, A.det().variables(), solution_dict=True)

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

    return list(set(det_res1 + det_res2 + eig_res1 + eig_res2))


def _factorize(A):
    n = len(A)
    for i, row in enumerate(A):
        if row.count(0) == (n-1) and row[i] != 0:
            return A.from_rows_and_columns([j for j in range(n) if j != i], [j for j in range(n) if j != i]), row[i]

    for i, col in enumerate(A.T):
        if col.count(0) == (n-1) and col[i] != 0:
            return A.from_rows_and_columns([j for j in range(n) if j != i], [j for j in range(n) if j != i]), row[i]
    return None


def _solve_simple_mat3(A):
    subA = _factorize(A)
    # case of permutational matrix
    if subA is None:
        raise NotImplementedError
    else:
        subA, varx = subA
        return _solve_simple_mat2(subA) + [varx == 1]


def solve_simple_mat(A):
    assert len(A) == len(A[0])
    if len(A) == 2:
        return _solve_simple_mat2(A)
    elif len(A) == 3:
        return _solve_simple_mat3(A)
    else:
        raise NotImplementedError


class SR_Degrees:
    def __init__(self, group_index, method="ascii", verbose=2, dim=2):
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

    def display(self, *args, use_pref=True) -> str:
        bod = self._display("") if self.disp == "ascii" else ""
        for el in args:
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

    def construct_congruences(self, A_inv, A):
        """Construct congruences for the subgroup check.

        The congruences have the following form:
            (A, t)(B, alpha(b))(A^{-1}, -A^{-1}t) = (ABA^{-1}, alpha(ABA^{-1})) mod L

        The lattice just multiplies on matrix, so the index equals to det(A):
            (A, t)(E, e)(A, t)^{-1} = (E, At)

        Calculate the conjugation:
            (A, t)(B, alpha(B))(A^{-1}, -A^{-1}t) =
            = (A, t)(BA^{-1}, -BA^{-1}t + alpha(B)) =

            = (ABA^{-1}, -ABA^{-1}t + A alpha(B) + t)
            = (ABA^{-1}, (E - ABA^{-1})t + A alpha(B)) = (ABA^{-1}, alpha(ABA^{-1})) mod L

            or

            A alpha(B) - alpha(ABA^{-1}) = (ABA^{-1} - E)t  mod L

        We can obviously abuse multiple appearence of ABA^{-1}.
        """
        G = self.G
        P = self.G.point_group_elements()
        snot = self.G.snot
        a0, a1 = var("a0 a1")
        x = matrix([[a0], [a1]])
        E = matrix(QQ, [[1, 0], [0, 1]])
        conds = []
        variables = []
        base_variables = []
        base_variables.extend(A_inv.variables())
        base_variables.extend([a0, a1])

        for i, g in enumerate(P):
            # hope sage simplifies to one from point group
            conj = (A * g * A_inv).simplify_rational()
            assert G.in_alpha(conj)
            # A alpha(B) - alpha(ABA^{-1})
            condition = A * G.alpha(g) - G.alpha(conj)
            condition = condition.simplify_rational()

            # (ABA^{-1} - E)t
            sym_res = (conj - E) * x

            # entire left part to print
            alpha_conj = (E - conj) * x + (A * G.alpha(g))
            for j in range(self.G.dim):
                alpha_conj[j, 0] = alpha_conj[j, 0].simplify_rational()

            if self.disp in ["latex", "markdown"]:
                self.print(
                    self.pref
                    + "(A, a)"
                    + self.display(block_matrix([[g, G.alpha(g)], [0, 1]]), use_pref=False)
                    + "(A^{-1}, -A^{-1}a) = "
                )
                self.print(self.display(block_matrix([[conj, alpha_conj], [0, 1]]), use_pref=False) + "=")
                self.print(self.display(block_matrix([[conj, G.alpha(conj)], [0, 1]]), use_pref=False) + self.pref)
            else:
                self.print(
                    ascii_art("\na\n\n")
                    + ascii_art(" ")
                    + ascii_art(snot[i])
                    + ascii_art("\na_inv\n\n")
                    + ascii_art(" ")
                    + ascii_art("\n=\n\n")
                    + ascii_art(" ")
                    + ascii_art(block_matrix([[conj, alpha_conj], [0, 1]]))
                    + ascii_art(" ")
                    + ascii_art("\n=\n\n")
                    + ascii_art(" ")
                    + ascii_art(
                        block_matrix(
                            [
                                [snot[i].linear_part(), G.alpha(snot[i].linear_part())],
                                [0, 1],
                            ]
                        )
                    )
                )
                self.print()

            conds.append(sym_res[0][0] - (condition[0][0] + var(f"n{i}")))
            conds.append(sym_res[1][0] - (condition[1][0] + var(f"m{i}")))

            variables.append(var(f"n{i}"))
            variables.append(var(f"m{i}"))
        return conds, base_variables, variables

    def solve_congruences_v3(self, conds, base_vars, variables):
        variables = base_vars
        # print(matrix(QQ, [[cond.coefficient(v) for v in variables] for cond in conds]))
        M = [[QmodZ(cond.coefficient(v)) for v in variables] for cond in conds]
        # print(M)
        gauss_elim(M, verbose=True)
        # print(M)
        # print(variables)
        # print(matrix(QQ, [[el.r if isinstance(el, Q_modZ_Element) else el for el in row] for row in M]))

        def to_qq(el):
            if isinstance(el, Q_modZ_Element):
                return el.r
            else:
                return QQ(el)

        # FIXME: works only if you believe that every system has a solution
        return {variables[i]: to_qq(M[i][i]) for i in range(len(variables))}

    def solve_congruences_v2(self, conds, base_vars, variables):
        # https://ask.sagemath.org/question/62549/solve-equation-of-matrices-over-integers/

        variables = base_vars + variables
        M = matrix(QQ, [[cond.coefficient(v) for v in variables] for cond in conds])

        r = matrix(QQ, [0 for _ in variables])

        # NotImplementedError: only integer lattices supported
        B, U = M.LLL(transformation=True)
        nz = sum(1 for r in B.rows() if r == 0)  # number of zero rows in B
        B = B.delete_rows(range(nz))  # trimming first nz rows of B
        U = U.delete_rows(range(nz))  # trimming first nz rows of U
        assert U * M == B  # the key property of U

        L = IntegerLattice(B, lll_reduce=False)  # our basis is already reduced and should not be altered
        assert r in L  # just in case checking that r belongs to L
        v = L.coordinate_vector(r) * U
        assert v * M == r
        return v

    def solve_congruences(self, conds, base_variables, variables):
        tmp = [con == 0 for con in conds]
        self.print("\nequations: ")
        printable = []
        for i, el in enumerate(tmp):
            if i % 4 == 0:
                printable.append(list())
            printable[-1].append(el)

        self.print(self.display(table(printable)))

        res = solve(tmp, *variables)
        self.print("\nanswer:")
        printable = []
        for i, el in enumerate(res[0]):
            if i % 4 == 0:
                printable.append([])
            printable[-1].append(el)
        self.print(self.display(table(printable)))

        if not res:
            self.print("Couldn't solve:", res)
        for cond in res[0]:
            # ugly check if n_i is rational
            if cond.left() in variables and not cond.right().is_integer() and not cond.right().variables():
                self.print(f"(*) Contradiction: {cond}!")
                return None
        return res

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

            A_inv = A.inverse().simplify_rational()
            if not G.is_symmorphic():
                eqs, base_vars, variables = self.construct_congruences(A_inv, A)
                conds = self.solve_congruences(eqs, base_vars, variables)
                if conds is None:
                    self.print("A doesn't form a virtual endomorphism.")
                    continue
                predicted_solutions = self.solve_congruences_v3(eqs, base_vars, variables)

                sc_degrees.add_solution(A, A_inv, conds, predicted_solutions)
            else:
                sc_degrees.add_solution(A, A_inv, [], {})
            self.print("Simplicity")
            self.print(self.pref + "A^{-1} = \n" + self.display(A_inv.simplify_rational(), use_pref=False) + self.pref)
            self.print("\neigenvalues:")
            self.print(self.display([el[0] for el in A_inv.charpoly().roots()]))
            self.print("charpoly:")
            chp = A_inv.charpoly()(SR("x"))
            chp = factor(chp)
            self.print(self.display(chp))
            self.print("\nindex of subgroup:")
            self.print(self.pref + "[G : H] = \n" + self.display(A.det(), use_pref=False) + self.pref)

        if self.disp == "latex" and norms:
            self.print("\\end{enumerate}")

        return sc_degrees

    def sr_degrees_calc(self) -> SR_DegreesSolution:
        pass

    def sr_degrees(self) -> SR_DegreesSolution:
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
        cond_table = [["matrix", "conds"]]
        pred_sols_table = [['matrix', 'sols']]
        for A, _, conds, pred_sols in sc_degrees:
            self.print(pref2)
            self.print(self.display("A = ", A))
            self.print("Determinant:")
            self.print(self.display(A.det(), self.in_sym, self.z))
            self.print("Conditions on endomorphism (self-coverings):")
            self.print(self.display(""))
            nontriv_conds = []
            for r in conds:
                for cond in r:
                    # skip trivial case int \in Z
                    if cond.right().is_integer():
                        continue

                    self.print(self.display(cond.right(), self.in_sym, self.z))
                    nontriv_conds.append(cond.right())

            self.print("\nPredicted solutions:")
            self.print(pred_sols)
            pred_sols_table.append([A, pred_sols])

            self.print("\nSelf-replicating degrees:")
            self.print(self.display(A.det(), self.neq, self.pm, 1))

            # condition that eigenvalues isnt equal to +-1
            cond1 = (A.det() + A.trace()).simplify_rational()
            cond2 = (A.det() - A.trace()).simplify_rational()
            self.print(self.display(cond1, self.neq, -1))
            self.print(self.display(cond2, self.neq, -1))

            # sympy handles quadratic diophantine equations pretty well:
            # https://www.alpertron.com.ar/METHODS.HTM
            # https://web.archive.org/web/20160323033111/http://www.jpr2718.org/ax2p.pdf
            #
            det_res1 = solve_diophantine(A.det() - 1, A.det().variables(), solution_dict=True)
            det_res2 = solve_diophantine(A.det() + 1, A.det().variables(), solution_dict=True)

            det_res1 = [tuple([val for _, val in sorted(sol.items())]) for sol in det_res1]
            det_res2 = [tuple([val for _, val in sorted(sol.items())]) for sol in det_res2]

            self.print("Solve diophantine determinant det(A) = +- 1:")
            self.print(self.display(det_res1))
            self.print(self.display(det_res2))

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

            self.print("Solve diophantine eigenvalues det(A) +- tr(A) = -1:")
            self.print(self.display(eig_res1))
            self.print(self.display(eig_res2))

            eig_table.append([A, A.det(), cond1 + 1, cond2 + 1])
            res_table.append([A, A.det(), sorted(list(set(det_res1 + det_res2 + eig_res1 + eig_res2)))])  # type: ignore
            cond_table.append([A, nontriv_conds])  # type: ignore
        if self.disp == "latex" and sc_degrees:
            self.print("\\end{enumerate}")
            self.print("\\newpage")

        return eig_table, res_table, cond_table, pred_sols_table
