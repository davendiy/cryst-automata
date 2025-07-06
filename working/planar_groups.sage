from self_similar import build_group
from normalizer import normalizers, to_L_basis


# latex, ascii_art, matrix, QQ, block_matrix, solve, var, SR, factor = ...

def permute(space, repeat) -> str:
    """ Get all the words of letters from the given
    space of given length.

    Parameters
    ----------
    space   : Iterable container of str elements of len 1.
    repeat  : int, length of yielded words
    allow_same_neighbours : if False then all words that don't have same elements
                            on the neighbour positions will be returned and all
                            possible words otherwise.
                            Default False
    Yields
    -------
    str, result word
    """
    if repeat == 1:
        for el in space:
            yield [el]
    elif repeat < 1:
        yield []
    else:
        for prev in permute(space, repeat - 1):
            for el in space:
                yield prev + [el]


def all_words(space, max_len=-1):
    """ Get all possible words (infinite set) of elements from the given
    space.

    Parameters
    ----------
    space   : Iterable container of str elements of len 1.
    allow_same_neighbours : if False then all words that don't have same elements
                            on the neighbour positions will be returned and all
                            possible words otherwise.
                            Default False.
    max_len : maximum allowed length of returned words.
              Default None, which means that there is no maximum length

    Yields
    -------
    str, result word
    """
    i = 1
    while True:
        for el in permute(space, repeat=i):
            yield el
        i += 1

        if max_len != -1 and i > max_len:
            break


def build_group(group_index, deep=3):
    e = matrix(QQ, [[1, 0], [0, 1]])
    G = to_L_basis(group_index, dim=2)
    gens = G.GeneratorsOfGroup()
    gens = [matrix(QQ, el) for el in gens]
    gens = [el for el in gens if el[:2, :2] != e]

    e = matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    alpha = {}
    snot = []
    for seq in all_words(gens, max_len=deep): 
        
        el = e 
        for el2 in seq: 
            el *= el2 
        
        g = el[:2, :2]
        t = el[:2, 2]
        if str(g) not in alpha:
            alpha[str(g)] = t
            snot.append(block_matrix([[g, t], [0, 1]]))

    P = G.PointGroup()
    P = [matrix(QQ, el) for el in P.AsList()]

    for mtx in P:
        assert str(mtx) in alpha, 'couldnt build every element of point group'

    return G, P, alpha, snot


class SR_Degrees: 

    def __init__(self, group_index, method='ascii'): 

        self.disp = method
        self.group_index = group_index
        if method == 'latex':
            self._display = latex
            self.title = '\\newpage\n\\subsection{Group %d}' % group_index
            self.section = '\\subsubsection{%s}'
            self.pref = '$$'
        else: 
            self._display = ascii_art
            self.title = f'=========================== Group {group_index} =================================='
            self.section = '\n----------------%s-------------------------\n'
            self.pref = '\n'

        self.G, self.P, self.alpha, self.snot = build_group(self.group_index)
        

    def display(self, *args, use_pref=True): 
        if use_pref:
            return self.pref + str(self._display(*args)) + self.pref
        else: 
            return str(self._display(*args))

    def header(self): 
        print(self.title)
        print('Generators of group:')
        print(self.display([matrix(QQ, el) for el in self.G.GeneratorsOfGroup()]))
        print('SNoT')
        print(self.display(self.snot))


    def algorithm(self): 
        self.header()
        
        G, P, alpha, snot = self.G, self.P, self.alpha, self.snot
        norms = normalizers(self.group_index, to_l_basis=True, 
                            dim=2, verbose=False, ignore_trivial=True)

        a0, a1 = var('a0 a1')
        x = matrix([[a0], [a1]])
        E = matrix(QQ, [[1, 0], [0, 1]])

        print(self.display(norms))
        print(self.section % 'Dilation')   

        if G.IsSymmorphicSpaceGroup(): 
            print(f"Group {self.group_index} is a semi-direct product, therefore the dilation part is trivial and only consists of integral vectors. ")

        if self.disp == 'latex': 
            print('\\begin{enumerate}')
            pref2 = '\\item'
        else: 
            pref2 = '...'


        for A_inv in norms:
            print(pref2 + ' testing inverse A (should have integral entities):')
            print(self.display(A_inv))
            
            A = A_inv.inverse()
            conds = []
            variables = []
            coeffs = []
            if not G.IsSymmorphicSpaceGroup(): 
                for i, g in enumerate(P): 
                    conj = (A * g * A_inv).simplify_rational()
                    # conj = g
                    condition = A_inv * alpha[str(g)] -  alpha[str(conj)]
                    sym_res = (conj - E) * x
                
                    # print(pref + '\\alpha(g) = ')
                    # print(display(alpha[str(g)], use_pref=False) + pref)
                    # print(pref + '\\alpha\\left(\\tau\\left(' + display(g) + '\\right)\\right) = ')
                    
                    alpha_conj = (E - conj) * x  + A_inv() * alpha[str(g)]
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
                    if self.disp == 'latex':
                        print(self.pref + '(A^{-1}, a)' + self.display(block_matrix([[g, alpha[str(g)]], [0, 1]]), use_pref=False) + '(A, -Aa) = ')
                        print(self.display(block_matrix([[conj, alpha_conj], [0, 1]]), use_pref=False) + '=')
                        print(self.display(block_matrix([[conj, alpha[str(conj)]], [0, 1]]), use_pref=False) + self.pref)
                    else: 
                        print(ascii_art('\na_inv\n\n') + ascii_art(' ')
                            + ascii_art(snot[i]) + ascii_art('\na\n\n') + ascii_art(' ')
                            + ascii_art('\n=\n\n') + ascii_art(' ') + ascii_art(T * snot[i] * T.inverse()) + ascii_art(' ')
                            + ascii_art('\n=\n\n') + ascii_art(' ') 
                            + ascii_art(block_matrix([[snot[i][:2, :2], alpha[str(snot[i][:2, :2])]], [0, 1]])))
                        print()
                        
                    conds.append(sym_res[0][0] - (condition[0][0] + var(f'n{i}')))
                    cond_coeffs = {}
                    # for variable in conds[-1].arguments(): 
                        # cur_coeff = conds[-1].coefficients(variable, sparse=False)
                        # cond_coeffs[variable] = cur_coeff[1]
                    # coeffs.append(cond_coeffs)
                    
                    conds.append(sym_res[1][0] - (condition[1][0] + var(f'm{i}')))
                    # cond_coeffs = {}
                    # for variable in conds[-1].arguments(): 
                        # cur_coeff = conds[-1].coefficients(variable, sparse=False)
                        # cond_coeffs[variable] = cur_coeff[1]
                    # coeffs.append(cond_coeffs)
                    
                    variables.append(var(f'n{i}'))
                    variables.append(var(f'm{i}'))
                
                tmp = [con == 0 for con in conds]
                print('\nequations: ')
                print(self.display(tmp))
                print('\nanswer:')
                res = solve(tmp, *variables)
                print(self.display(*res))
            
            print('Simplicity')
            print(self.pref + 'A = \n' + self.display(A.simplify_rational(), use_pref=False) + self.pref)
            print('\neigenvalues:')
            print(self.display([el[0] for el in A.charpoly().roots()]))
            print('charpoly:')
            chp = A.charpoly()(SR('x'))
            chp = factor(chp)
            print(self.display(chp))
            print('\nindex of subgroup:')
            print(self.pref + '[G : H] = \n' + self.display(A_inv.det(), use_pref=False) + self.pref)

        if self.disp == 'latex': 
            print('\\end{enumerate}')
