import os
from collections import deque

import numpy as np
from sage.all import (
    QQ,
    block_matrix,
    copy,
    gap,
    matrix,
    MatrixGroup,
    sign,
    Permutation,
)

from .normalizer import normalizers

MAX_ITERATIONS = 1_000_000


def prepare_gap_env(use_3d_gap=True):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, "el2word.g")

    # change gap version to installed separately. Also use maximum 4Gb memory
    if use_3d_gap:
        os.environ["SAGE_GAP_COMMAND"] = "~/gap/gap -s 4G"
    gap(f'LoadPackage("Cryst");;Read("{file_path}")')


def build_finite_group(gens, trivial, max_iter=MAX_ITERATIONS):
    gens_extended = [trivial] + gens + [el.inverse() for el in gens]
    names = (
        [0]
        + [i for i in range(1, len(gens) + 1)]
        + [-i for i in range(1, len(gens) + 1)]
    )

    found = {str(gen): [name] for gen, name in zip(gens_extended, names)}

    q = deque()
    q.extend(gens_extended)

    _n = 0
    while q:
        if _n > max_iter:
            raise ValueError("seems infinite group")

        x = q.popleft()
        for y, y_name in zip(gens_extended, names):
            tmp1 = x * y
            tmp2 = y * x

            if str(tmp1) not in found:
                found[str(tmp1)] = found[str(x)] + [y_name]
                q.append(tmp1)

            if str(tmp2) not in found:
                found[str(tmp2)] = [y_name] + found[str(x)]
                q.append(tmp2)

    return found


def from_indices_list(gens, triv, seq):
    res = copy(triv)
    for mul in seq:
        if mul:
            gen = gens[abs(mul) - 1]
            gen = gen if mul > 0 else gen.inverse()
            res *= gen
    return res


def _to_L_basis(n, dim=3):
    """Change basis of the crystallographic group to transform the lattice
    L to the Z^dim

    Parameters
    ----------
    n : int
        number of the crystallographic group
    dim : int
        dimension

    Returns
    -------
    MatrixGroup with new generators
    """
    G = gap(f"SpaceGroupOnLeftIT({dim}, {n})")
    gens = [matrix(QQ, el) for el in G.GeneratorsOfGroup()]

    v = matrix(QQ, G.TranslationBasis()).T
    trans = matrix(QQ, [0 for _ in range(dim)]).T
    conj = block_matrix(QQ, [[v, trans], [0, 1]])
    new_gens = [conj.inverse() * el * conj for el in gens]
    return gap.AffineCrystGroupOnLeft(new_gens)


class SpaceGroup_Element:
    def __init__(self, mtx):
        if isinstance(mtx, SpaceGroup_Element):
            mtx = copy(mtx._body)

        assert mtx.dimensions()[0] == mtx.dimensions()[1]
        self._body = mtx

    def inverse(self):
        return type(self)(self._body.inverse())

    def __eq__(self, value):
        return self._body == value._body

    def __mul__(self, other):
        other = SpaceGroup_Element(other)

        if self.dim != other.dim:
            raise ValueError(
                "can't multiply space group elements of different dimensions"
            )

        return SpaceGroup_Element(self._body * other._body)

    def __rmul__(self, other):
        other = SpaceGroup_Element(other)

        if self.dim != other.dim:
            raise ValueError(
                "can't multiply space group elements of different dimensions"
            )

        return SpaceGroup_Element(other._body * self._body)

    @property
    def dim(self):
        return self._body.dimensions()[0]

    def linear_part(self):
        n = self.dim
        return self._body[: n - 1, : n - 1]

    def translation(self):
        n = self.dim
        return self._body[: n - 1, -1]

    @classmethod
    def construct_element(cls, linear, trans):
        n = linear.dimensions()[0] + 1
        res = matrix.identity(QQ, n)
        res[: n - 1, : n - 1] = linear
        res[: n - 1, -1] = trans
        return cls(res)

    @classmethod
    def from_linear(cls, linear):
        n = linear.dimensions()[0] + 1
        res = matrix.identity(QQ, n)
        res[: n - 1, : n - 1] = linear
        return cls(res)
    
    def __pow__(self, n):
        return self.__class__(self._body ** n) 

    def __repr__(self):
        return repr(self._body)

    def __str__(self):
        return str(self._body)

    @classmethod
    def from_translation(cls, trans):
        if isinstance(trans, list):
            trans = np.array(trans).flatten()
        elif isinstance(trans, np.ndarray):
            pass
        else:
            trans = trans.numpy().flatten()
        n = trans.shape[0] + 1
        res = matrix.identity(QQ, n)
        res[: n - 1, -1] = matrix(QQ, trans.reshape(n - 1, 1))
        return cls(res)

    @classmethod
    def from_gap_element(cls, el):
        return cls(matrix(QQ, el))


class SpaceGroup_gap:
    max_iterations = MAX_ITERATIONS

    def __init__(self, gap_G, dim, ita_num):
        self.G = gap_G
        self.dim = dim
        self.ita_num = ita_num

        self.P = self.G.PointGroup()

        self.P_triv = matrix.identity(QQ, self.dim)
        self.G_triv = SpaceGroup_Element(matrix.identity(QQ, self.dim + 1))

        self.G_gens = [
            SpaceGroup_Element.from_gap_element(el) for el in self.G.GeneratorsOfGroup()
        ]
        self.G_nontriv = [el for el in self.G_gens if el.linear_part() != self.P_triv]
        self._P_names = [f'g_{i+1}' for i in range(len(self.G_nontriv))]
        self._L_names = [f'e_{i+1}' for i in range(dim)]

        self.L_gens = [SpaceGroup_Element.from_translation(row) 
                       for row in matrix(QQ, self.G.TranslationBasis())]

        self._gen2name = {}
        self._name2gen = {}
        self._rebuild_names()

        self.G_sorted_gens = self.G_nontriv + self.L_gens
        self.P_gens = [el.linear_part() for el in self.G_nontriv]

        self._P_dict = build_finite_group(self.P_gens, self.P_triv, self.max_iterations)

        self._alpha = {str(self.P_triv) : self.G_triv.translation()}
        self.snot = [self.G_triv]

        for el, seq in self._P_dict.items():
            val = from_indices_list(self.G_nontriv, self.G_triv, seq)
            val = SpaceGroup_Element(val)
            assert str(val.linear_part()) == str(el), f"{val.linear_part()}, {el}"

            self._alpha[str(el)] = val.translation()
            self.snot.append(val)

        self._in_lattice_precalc = None

    def _rebuild_names(self): 
        self._gen2name = {}
        for name, el in zip(self._P_names + self._L_names, 
                            self.G_nontriv + self.L_gens): 
            self._name2gen[name] = el
            self._name2gen[f'{name}^(-1)'] = el.inverse()
            self._gen2name[str(el)] = name
            self._gen2name[str(el.inverse())] = f'{name}^(-1)'
    
    def in_alpha(self, sym):
        return str(sym) in self._alpha

    def point_group_normalizer(self, **kwargs): 
        P = MatrixGroup(self.P_gens)
        return normalizers(P, **kwargs)

    def in_lattice_basis(self):
        """Returns whether L == ZZ^n.""" 
        if self._in_lattice_precalc is None: 
            self._in_lattice_precalc = True
            for i in range(self.dim):
                tr = [0 for _ in range(self.dim)]
                tr[i] = 1
                if self.L_gens[i] != SpaceGroup_Element.from_translation(tr): 
                    self._in_lattice_precalc = False
                    break

        return self._in_lattice_precalc

    def alpha(self, el: SpaceGroup_Element):
        p = el.linear_part() 
        if not self.in_alpha(p):
            raise ValueError(f"element: {el} is not in the group.")
        return self._alpha[str(p)]

    def __contains__(self, el): 
        if not isinstance(el, SpaceGroup_Element): 
            return False 
        
        return self.contains(el)

    def point_group_elements(self): 
        return [el.linear_part() for el in self.snot]

    @staticmethod
    def _is_integral(translation): 
        for row in translation: 
            for el in row: 
                if not el.is_integer(): 
                    return False 
        return True

    def contains(self, el: SpaceGroup_Element): 
        if not self.in_lattice_basis(): 
            raise NotImplementedError()
        
        p = el.linear_part() 
        if not self.in_alpha(p):
            return False 

        tr = el.translation()
        tr_alpha = self.alpha(el)

        diff = tr - tr_alpha

        return self._is_integral(diff)
    
    def as_word(self, el: SpaceGroup_Element, readable=True): 
        if not self.in_lattice_basis(): 
            raise NotImplementedError()

        if el not in self: 
            raise ValueError(f"element {el} is not in the group.")

        p = el.linear_part()
    
        tr = el.translation()
        tr_alpha = self.alpha(el)
        diff = tr - tr_alpha 

        res_word = []

        if readable:
            for idx in self._P_dict[str(p)]: 
                if not idx:
                    continue
                g = self.G_nontriv[abs(idx) - 1]

                g = g if idx > 0 else g.inverse()
                res_word.append(self._gen2name[str(g)])
        else: 
            res_word += self._P_dict[str(p)]

        for i, x in enumerate(self._trans_entities(diff)): 
            if not x:
                continue
                
            if readable:
                res_word.append(f'{self._L_names[i]}^({x})')
            else: 
                idx = len(self.P_gens) + i + 1
                res_word = [sign(x) * idx for _ in range(abs(x))] + res_word

        return res_word

    def random_element(self, max_length=100): 
        import random 

        gens_n = len(self.G_sorted_gens) - 1
        it_seq_n = random.randint(1, max_length) 
        it_seq = [random.randint(-gens_n, gens_n) for _ in range(it_seq_n)]
        return from_indices_list(self.G_sorted_gens, self.G_triv, it_seq)

    def _trans_entities(self, tr): 
        assert tr.dimensions()[0] == self.dim
        return [row[0] for row in tr]

    def set_P_names(self, names): 
        if len(names) != len(self.G_nontriv): 
            raise ValueError("incorrect amount of names.")

        self._P_names = names
        self._rebuild_names()

    def set_L_names(self, names): 
        if len(names) != len(self.L_gens): 
            raise ValueError("incorrect amount of names.")

        self._L_names = names
        self._rebuild_names()

    def is_symmorphic(self):
        return self.G.IsSymmorphicSpaceGroup()

    def self_similar(self, T, verbose=False, safe=True):
        """Construct self-similar action for a crystallographic group
        given element of affine group that is conjugation for virtual
        endomorphism construction.

        Parameters
        ----------
        T : matrix of dimension `dim` + 1
            an element of the affine group A(dim) that can be represented as
            (M + t), where M is a `dim`x`dim` matrix and t is a vector that represents
            translation. Should be given in a matrix form, i.e. (M + t) is a
            block matrix:
                            _           _
                            ||     |     ||
                            ||  M  |  t  ||
                            ||_____|_____||
                            ||  0  |  1  ||
                            -           -
            dim : int
            dimension of euclidean space, where we consider a crystallographic
            group
        verbose : bool
            True to show auxiliary messages
        safe : bool
            if True, then function raises error if `T` doesn't generate
            virtual endomorphism.
        
        Returns
        -------
        dict { (a, i): [j, b] }
            a self-similar action
        """

        if not self.in_lattice_basis(): 
            raise NotImplementedError()

        def phi(_g): 
            return T * _g * T.inverse() 
        
        def phi_inv(_g): 
            return T.inverse() * _g * T 

        if verbose:
            print("=====================================================================")

        # check whether there exists virtual endomorphism
        gens_H = []
        for el in self.G_gens:
            conj = phi_inv(el)
            if verbose:
                print("\nconjugate el:")
                print(conj)
                print("conj in G:", conj in self)
            if conj not in self and safe:
                raise ValueError("Bad matrix T, there is no virtual endomorphism")
            elif conj not in self:
                print("Bad matrix T, there is no virtual endomorphism")
                return -1

            gens_H.append(conj._body)

        # create subgroup as image of virtual endomorphism
        H = self.G.Subgroup(gens_H)

        if verbose:
            print("----------------------------------------------------")
            print("Index of subgroup H:", self.G.Index(H))

        trans = self.G.RightTransversal(H)   # походу треба LeftTransversal
        trans_els = [SpaceGroup_Element(matrix(QQ, el)) for el in trans.AsList()]
        if verbose:
            print("Transversal:")
            print(*trans_els, sep='\n\n')

        # create self-similar action
        res_map = {}
        for a in self.G_gens:

            # NOTE: ENUMERATION STARTS FROM 1
            for i, d_i in enumerate(trans_els, 1):
                adi = a * d_i

                if verbose:
                    print(f"{self._gen2name[str(a)]}d_{i}:")
                    print(adi, end='\n\n')

                # (a * d_i)^{-1} * (a * d_i) = e \in H
                # ==> d_j^{-1} = (a * d_i)^{-1}
                #
                # firstly, find coset for (a * d_i)^{-1}

                # замість цього треба самому написати
                d_j_index = trans.PositionCanonical(adi.inverse()._body)

                # then get d_j^{-1}^{-1}
                d_j = trans_els[int(d_j_index) - 1].inverse()

                if (d_j.inverse() * a * d_i)._body not in H:
                    raise ValueError(f"This shouln't happen--wrong coset: {(a, d_i, d_j)}")

                # conjugation in the right direction, i.e. apply \phi(d_j^{-1} a d_i)
                tmp_res = phi(d_j.inverse() * a * d_i)
                
                if verbose: 
                    print('image:', tmp_res, sep='\n')

                res_map[(self._gen2name[str(a)], i)] = (int(d_j_index), tmp_res)
                                                
        return res_map

    @classmethod
    def from_gap_cryst(cls, ita_num, dim=3, change_basis=True):
        """Construct a crystallographic group, using GAP Cryst package.

        Parameters
        ----------
        ita_num : int
            number of the crystallographic group
        dim : int
            dimension
        change_basis : bool
            Change basis of the crystallographic group to get the lattice
            L equal to ZZ^dim

        """
        G = gap(f"SpaceGroupOnLeftIT({dim}, {ita_num})")

        if change_basis:
            gens = [matrix(QQ, el) for el in G.GeneratorsOfGroup()]
            v = matrix(QQ, G.TranslationBasis()).T
            trans = matrix(QQ, [0 for _ in range(dim)]).T
            conj = block_matrix(QQ, [[v, trans], [0, 1]])
            new_gens = [conj.inverse() * el * conj for el in gens]
            return cls(gap.AffineCrystGroupOnLeft(new_gens), dim, ita_num)
        else: 
            return cls(G, dim, ita_num)

    def to_lattice_basis(self): 
        if self.in_lattice_basis(): 
            return self
        else: 
            return self.__class__.from_gap_cryst(
                self.ita_num, self.dim, change_basis=True
            ) 

    def _wreath_recursion(self, action_dict):

        n = max(k[1] for k in action_dict)
        
        w_el = {}
        w_p = {}
        for k, v in action_dict.items(): 
            el, let = k 
            if el not in w_el:
                w_el[el] = [None for _ in range(n)]
            if el not in w_p: 
                w_p[el] = [None for _ in range(n)]

            w_el[el][let-1] = ''.join(self.as_word(v[1], readable=True))
            w_p[el][let-1] = v[0]

        print('Вінцева рекурсія:')
        for el in w_el.keys(): 
            p = Permutation(w_p[el])
            if p.to_permutation_group_element().order() == 1: 
                p = '()'
            else: 
                p = str(p.cycle_tuples())
            
            print(f'$${el} = {p}{w_el[el]}$$'.replace("'", '')\
                  .replace('[', '(').replace(']', ')')\
                    .replace('(-1)', '{-1}')\
                    .replace('^(1)', ''))

        print('```')
        print(self._name2gen)
        print('```')    
