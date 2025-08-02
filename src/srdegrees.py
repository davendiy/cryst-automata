from sage.all import QQ, SR, ascii_art, block_matrix, factor, latex, matrix, solve, var
from sage.modules.free_module_integer import IntegerLattice


from .space_groups import SpaceGroup_gap


class SR_Degrees:
    def __init__(self, group_index, method="ascii"):
        self.disp = method
        self.group_index = group_index
        if method == "latex":
            self._display = latex
            self.title = "\\subsection{Group %d}" % group_index
            self.section = "\\subsubsection{%s}"
            self.pref = "$$"
            self.in_sym = "\\in"
            self.notin_sym = "\\notin"
            self.z = "\\mathbb{Z}"
        else:
            self._display = ascii_art
            self.title = f"=========================== Group {group_index} =================================="
            self.section = "\n----------------%s-------------------------\n"
            self.pref = "\n"
            self.in_sym = "  ∈  "
            self.notin_sym = "  /∈  "
            self.z = "Z"

        self.G = SpaceGroup_gap.from_gap_cryst(group_index, dim=2, change_basis=True)

    def display(self, *args, use_pref=True) -> str:
        bod = "" if self.disp == "latex" else self._display("")
        for el in args:
            if isinstance(el, str) and self.disp == "latex":
                bod += el
            else:
                bod += self._display(el)

        if use_pref:
            return self.pref + str(bod) + self.pref
        else:
            return str(bod)

    def header(self):
        print(self.title)
        print("Generators of group:")
        print(self.display([matrix(QQ, el) for el in self.G.gap_G.GeneratorsOfGroup()]))
        print("SNoT")
        print(self.display(self.G.snot))

    def texdoc_header(self):
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
        if self.disp == "latex":
            print("\\end{document}")

    def construct_congruences(self, A, A_inv, include_Avar=False):
        G = self.G
        P = self.G.point_group_elements()
        snot = self.G.snot
        a0, a1 = var("a0 a1")
        x = matrix([[a0], [a1]])
        E = matrix(QQ, [[1, 0], [0, 1]])
        conds = []
        variables = []
        base_variables = []
        base_variables.extend(A.variables())
        base_variables.extend([a0, a1])

        for i, g in enumerate(P):
            conj = (A * g * A_inv).simplify_rational()
            # conj = g
            condition = A_inv * G.alpha(g) - G.alpha(conj)
            condition = condition.simplify_rational()

            sym_res = (conj - E) * x

            alpha_conj = (E - conj) * x + (A_inv * G.alpha(g))
            for j in range(self.G.dim):
                alpha_conj[j, 0] = alpha_conj[j, 0].simplify_rational()

            if self.disp == "latex":
                print(
                    self.pref
                    + "(A^{-1}, a)"
                    + self.display(
                        block_matrix([[g, G.alpha(g)], [0, 1]]), use_pref=False
                    )
                    + "(A, -Aa) = "
                )
                print(
                    self.display(
                        block_matrix([[conj, alpha_conj], [0, 1]]), use_pref=False
                    )
                    + "="
                )
                print(
                    self.display(
                        block_matrix([[conj, G.alpha(conj)], [0, 1]]), use_pref=False
                    )
                    + self.pref
                )
            else:
                print(
                    ascii_art("\na_inv\n\n")
                    + ascii_art(" ")
                    + ascii_art(snot[i])
                    + ascii_art("\na\n\n")
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
                print()

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

        L = IntegerLattice(
            B, lll_reduce=False
        )  # our basis is already reduced and should not be altered
        assert r in L  # just in case checking that r belongs to L
        v = L.coordinate_vector(r) * U
        assert v * M == r
        return v

    def solve_congruences(self, conds, base_variables, variables):
        tmp = [con == 0 for con in conds]
        print("\nequations: ")
        print(self.display(tmp))
        print("\nanswer:")
        res = solve(tmp, *variables)
        print(self.display(*res))

        if not res:
            print("Couldn't solve:", res)
        for cond in res[0]:
            # ugly check if n_i is rational
            if (
                cond.left() in variables
                and not cond.right().is_integer()
                and not cond.right().variables()
            ):
                print(f"[*] Contradiction: {cond}!")
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

        print(self.display(norms))
        print(self.section % "Dilation")

        if G.is_symmorphic():
            print(
                f"Group {self.group_index} is a semi-direct product, therefore the dilation part is trivial and only consists of integral vectors. "
            )

        if self.disp == "latex" and norms:
            print("\\begin{enumerate}")
            pref2 = "\\item"
        else:
            pref2 = "..."

        sc_degrees = {}

        for A_inv in norms:
            print(pref2 + " testing inverse A (should have integral entities):")
            print(
                self.pref
                + "A^{-1} = \n"
                + self.display(A_inv, use_pref=False)
                + self.pref
            )

            A = A_inv.inverse().simplify_rational()
            if not G.is_symmorphic():
                eqs, base_vars, variables = self.construct_congruences(A, A_inv)
                conds = self.solve_congruences(eqs, base_vars, variables)
                if conds is None:
                    print("A doens't form a virtual endomorphism.")
                    continue
                sc_degrees[str(A)] = A, A_inv, conds
            else:
                sc_degrees[str(A)] = A, A_inv, []
            print("Simplicity")
            print(
                self.pref
                + "A = \n"
                + self.display(A.simplify_rational(), use_pref=False)
                + self.pref
            )
            print("\neigenvalues:")
            print(self.display([el[0] for el in A.charpoly().roots()]))
            print("charpoly:")
            chp = A.charpoly()(SR("x"))
            chp = factor(chp)
            print(self.display(chp))
            print("\nindex of subgroup:")
            print(
                self.pref
                + "[G : H] = \n"
                + self.display(A_inv.det(), use_pref=False)
                + self.pref
            )

        if self.disp == "latex" and norms:
            print("\\end{enumerate}")
        return sc_degrees

    def sc_degrees(self):
        return self.algorithm()

    def sr_degrees(self):
        sc_degrees = self.algorithm()

        print(self.section % "Self-replicating degrees.")

        if self.disp == "latex" and sc_degrees:
            print("\\begin{enumerate}")
            pref2 = "\\item"
        else:
            pref2 = "..."

        for _, (A, A_inv, conds) in sc_degrees.items():
            print(pref2)
            print(self.display("A^{-1} = ", A_inv))
            print("Determinant:")
            print(self.display(A_inv.det(), self.in_sym, self.z))
            print("Conditions on endomorphism (self-coverings):")
            for r in conds:
                for cond in r:
                    if cond.right().is_integer():
                        continue

                    print(self.display(cond.right(), self.in_sym, self.z))
            f = A.charpoly()
            print("Coefficients of the charpoly:")
            for c in f.coefficients():
                print(self.display(c, self.notin_sym, self.z))
            print("Conditions of eigenvalues:")
            for e in f.roots():
                print(self.display(e, self.notin_sym, self.z))

        if self.disp == "latex" and sc_degrees:
            print("\\end{enumerate}")
            print("\\newpage")
