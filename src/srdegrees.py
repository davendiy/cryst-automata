from sage.all import QQ, SR, ascii_art, block_matrix, factor, latex, matrix, solve, table, var
from sage.modules.free_module_integer import IntegerLattice


from .space_groups import SpaceGroup_gap


class SR_Degrees:
    def __init__(self, group_index, method="ascii", verbose=2):
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

        self.G = SpaceGroup_gap.from_gap_cryst(group_index, dim=2, change_basis=True)

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

    def solve_congruences(self, conds, _, variables):
        tmp = [con == 0 for con in conds]
        self.print("\nequations: ")
        self.print(self.display(tmp))
        self.print("\nanswer:")
        res = solve(tmp, *variables)
        self.print(self.display(*res))

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

    def algorithm(self):
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

        sc_degrees = {}

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
                sc_degrees[str(A)] = A_inv, A, conds
            else:
                sc_degrees[str(A)] = A_inv, A, []
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

    def sc_degrees(self):
        return self.algorithm()

    def sr_degrees(self):
        sc_degrees = self.algorithm()

        self.print(self.section % "Self-replicating degrees.")

        if self.disp == "latex" and sc_degrees:
            self.print("\\begin{enumerate}")
            pref2 = "\\item"
        elif self.disp == 'markdown':
            pref2 = '--- '
        else:
            pref2 = "..."

        res_table = [["matrix", "det", "eigenvalue != 1", "eigenvalue != -1", "conds"]]
        for _, (_, A, conds) in sc_degrees.items():
            self.print(pref2)
            self.print(self.display("A = ", A))
            self.print("Determinant:")
            self.print(self.display(A.det(), self.in_sym, self.z))
            self.print("Conditions on endomorphism (self-coverings):")
            for r in conds:
                for cond in r:
                    if cond.right().is_integer():
                        continue
                    self.print(self.display(cond.right(), self.in_sym, self.z))

            self.print("Self-replicating degrees:")
            self.print(self.display(A.det(), self.neq, self.pm, 1))
            x0, x3 = A[0, 0], A[1, 1]

            # condition that eigenvalues isnt equal to +-1
            cond1 = (A.det() + x0 + x3).simplify_rational()
            cond2 = (A.det() - x0 - x3).simplify_rational()
            self.print(self.display(cond1, self.neq, -1))
            self.print(self.display(cond2, self.neq, -1))
            res_table.append([A, A.det(), cond1, cond2, conds])
        if self.disp == "latex" and sc_degrees:
            self.print("\\end{enumerate}")
            self.print("\\newpage")

        return table(res_table, header_row=True)
