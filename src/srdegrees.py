
from sage.all import (QQ, SR, ascii_art, block_matrix, factor, latex, matrix,
                      solve, var)

from .space_groups import SpaceGroup_gap


class SR_Degrees:
    def __init__(self, group_index, method="ascii"):
        self.disp = method
        self.group_index = group_index
        if method == "latex":
            self._display = latex
            self.title = "\\newpage\n\\subsection{Group %d}" % group_index
            self.section = "\\subsubsection{%s}"
            self.pref = "$$"
        else:
            self._display = ascii_art
            self.title = f"=========================== Group {group_index} =================================="
            self.section = "\n----------------%s-------------------------\n"
            self.pref = "\n"

        self.G = SpaceGroup_gap.from_gap_cryst(group_index, dim=2, change_basis=True)

    def display(self, *args, use_pref=True):
        if use_pref:
            return self.pref + str(self._display(*args)) + self.pref
        else:
            return str(self._display(*args))

    def header(self):
        print(self.title)
        print("Generators of group:")
        print(self.display([matrix(QQ, el) for el in self.G.gap_G.GeneratorsOfGroup()]))
        print("SNoT")
        print(self.display(self.G.snot))

    def construct_congruences(self, A, A_inv):
        G = self.G
        P = self.G.point_group_elements()
        snot = self.G.snot
        a0, a1 = var("a0 a1")
        x = matrix([[a0], [a1]])
        E = matrix(QQ, [[1, 0], [0, 1]])
        conds = []
        variables = []

        for i, g in enumerate(P):
            conj = (A * g * A_inv).simplify_rational()
            # conj = g
            condition = A_inv * G.alpha(g) - G.alpha(conj)
            sym_res = (conj - E) * x

            # print(pref + '\\alpha(g) = ')
            # print(display(alpha[str(g)], use_pref=False) + pref)
            # print(pref + '\\alpha\\left(\\tau\\left(' + display(g) + '\\right)\\right) = ')

            alpha_conj = (E - conj) * x + A_inv() * G.alpha(g)
            # if disp == 'latex':
            #     print(pref + '(A^{-1}, a)' + display(block_matrix([[g, alpha[str(g)]], [0, 1]]), use_pref=False) + '(A, -Aa) = ')
            #     print(display(block_matrix([[conj, alpha_conj], [0, 1]]), use_pref=False) + '=')
            #     print(display(block_matrix([[conj, alpha[str(conj)]], [0, 1]]), use_pref=False) + pref)
            # else:
            #     print(ascii_art('\na_inv\n\n') + ascii_art(' ')
            #           + ascii_art(block_matrix([[g, alpha[str(g)]], [0, 1]])) + ascii_art('\na\n\n') + ascii_art(' ')
            #           + ascii_art('\n=\n\n') + ascii_art(' ') + ascii_art(block_matrix([[conj, alpha_conj], [0,1]])) + ascii_art(' ')
            #           + ascii_art('\n=\n\n') + ascii_art(' ')
            #           + ascii_art(block_matrix([[conj, alpha[str(conj)]], [0, 1]])))
            #     print()

            T = block_matrix([[A_inv, x], [0, 1]])
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
                    + ascii_art(T * snot[i] * T.inverse())
                    + ascii_art(" ")
                    + ascii_art("\n=\n\n")
                    + ascii_art(" ")
                    + ascii_art(
                        block_matrix(
                            [[snot[i].linear_part(), G.alpha(snot[i].linear_part())], [0, 1]]
                        )
                    )
                )
                print()

            conds.append(sym_res[0][0] - (condition[0][0] + var(f"n{i}")))
            cond_coeffs = {}
            # for variable in conds[-1].arguments():
            # cur_coeff = conds[-1].coefficients(variable, sparse=False)
            # cond_coeffs[variable] = cur_coeff[1]
            # coeffs.append(cond_coeffs)

            conds.append(sym_res[1][0] - (condition[1][0] + var(f"m{i}")))
            # cond_coeffs = {}
            # for variable in conds[-1].arguments():
            # cur_coeff = conds[-1].coefficients(variable, sparse=False)
            # cond_coeffs[variable] = cur_coeff[1]
            # coeffs.append(cond_coeffs)

            variables.append(var(f"n{i}"))
            variables.append(var(f"m{i}"))
        return conds, variables

    def solve_congruences(self, conds, variables):
        tmp = [con == 0 for con in conds]
        print("\nequations: ")
        print(self.display(tmp))
        print("\nanswer:")
        res = solve(tmp, *variables)
        print(self.display(*res))
        return res

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

        if self.disp == "latex":
            print("\\begin{enumerate}")
            pref2 = "\\item"
        else:
            pref2 = "..."

        for A_inv in norms:
            print(pref2 + " testing inverse A (should have integral entities):")
            print(self.display(A_inv))

            A = A_inv.inverse()
            if not G.is_symmorphic():
                conds, variables = self.construct_congruences(A, A_inv)
                res = self.solve_congruences(conds, variables)

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

        if self.disp == "latex":
            print("\\end{enumerate}")
