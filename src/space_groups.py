import os
from collections import deque

import numpy as np
from sage.all import (
    QQ,
    block_matrix,
    copy,
    gap,
    matrix,
    sign,
)


MAX_ITERATIONS = 1_000_000


def prepare_gap_env(use_3d_gap=True):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, "el2word.g")

    # change gap version to installed separately. Also use maximum 4Gb memory
    if use_3d_gap:
        os.environ["SAGE_GAP_COMMAND"] = "~/gap/gap -s 4G"
    gap(f'LoadPackage("Cryst");;Read("{file_path}")')


def to_L_basis(n, dim=3):
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

    def __init__(self, gap_G, dim):
        self.G = gap_G
        self.dim = dim

        self.P = self.G.PointGroup()

        self.P_triv = matrix.identity(QQ, self.dim)
        self.G_triv = SpaceGroup_Element(matrix.identity(QQ, self.dim + 1))

        self.G_gens = [
            SpaceGroup_Element.from_gap_element(el) for el in self.G.GeneratorsOfGroup()
        ]
        self.G_nontriv = [el for el in self.G_gens if el.linear_part() != self.P_triv]
        self._P_names = [f'g_{i+1}' for i in range(len(self.G_nontriv))]
        self._L_names = [f'e_{i+1}' for i in range(dim)]

        self.L_gens = []

        # TODO: wrong assumption that L == ZZ^n
        for i in range(dim):
            tr = [0 for _ in range(dim)]
            tr[i] = 1
            self.L_gens.append(SpaceGroup_Element.from_translation(tr))

        self._gen2name = {}
        self._name2gen = {}
        self._rebuild_names()

        self.G_sorted_gens = self.G_nontriv + self.L_gens
        self.P_gens = [el.linear_part() for el in self.G_nontriv]

        self._P_dict = build_finite_group(self.P_gens, self.P_triv, self.max_iterations)

        self._alpha = {str(self.P_triv) : self.G_triv.translation()}
        self.snot = []

        for el, seq in self._P_dict.items():
            val = from_indices_list(self.G_nontriv, self.G_triv, seq)
            val = SpaceGroup_Element(val)
            assert str(val.linear_part()) == str(el), f"{val.linear_part()}, {el}"

            self._alpha[str(el)] = val.translation()
            self.snot.append(val)
    
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

    def is_normalized(self):
        """Returns whether L == ZZ^n.""" 
        return True

    def alpha(self, el: SpaceGroup_Element):
        p = el.linear_part() 
        if not self.in_alpha(p):
            raise ValueError(f"element: {el} is not in the group.")
        return self._alpha[str(p)]

    def __contains__(self, el): 
        if not isinstance(el, SpaceGroup_Element): 
            return False 
        
        return self.contains(el)

    @staticmethod
    def _is_integral(translation): 
        for row in translation: 
            for el in row: 
                if not el.is_integer(): 
                    return False 
        return True

    def contains(self, el: SpaceGroup_Element): 
        if not self.is_normalized(): 
            raise NotImplementedError()
        
        p = el.linear_part() 
        if not self.in_alpha(p):
            return False 

        tr = el.translation()
        tr_alpha = self.alpha(el)

        diff = tr - tr_alpha

        return self._is_integral(diff)
    
    def as_word(self, el: SpaceGroup_Element, readable=True): 
        if not self.is_normalized(): 
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

        if change_basis:
            return cls(to_L_basis(ita_num, dim=dim), dim)
        else:
            raise NotImplementedError()


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
