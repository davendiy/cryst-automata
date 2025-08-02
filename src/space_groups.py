from __future__ import annotations

import os
from collections import deque
from itertools import product
from typing import Iterable

import numpy as np
from sage.all import (
    QQ,
    MatrixGroup,
    Permutation,
    Subsets,
    ascii_art,
    block_matrix,
    copy,
    gap,
    lcm,
    matrix,
    sign,
    var,
    vector,
)

from .normalizer import normalizers

MAX_ITERATIONS = 1_000_000


# TODO:
#   - [x] possibility to create SpaceGroup from only generators
#   - [x] function to get lattice basis
#   - [x] use lattice basis to check whether g in G for arbitrary G:
#     g in G <=>  L(g) in L(G)
#   - [x] use previous in self-similar method instead of Gap's one
#   - [x] LeftTransversal instead of Gap's RightTransvesal
#   - [x] check simplicity of the given virtual endomorphism
#   - [ ] add method to check whether cryst groups isomorphic
#   - [ ] find min isomorphic subgroup??????


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


def check_div(A):
    """Check whether matrix is simply factorizable.

    Exploits the fact that determinant of (A - xE) matrix (i.e. characteristic polynomial)
    can be found using row/column decomposition. Thus, if a row/column contains
    n-1 zero, then the determinant can be found as x * det(A'), where A' is a
    submatrix.

    Parameters
    ----------
    A : symbolic matrix

    Returns
    -------
    bool
        True if characteristic polynomial is simply factorizable
    """
    x = var("x")

    n = len(list(A))
    test1 = A - x * matrix.identity(n)
    test2 = test1.T()

    # check rows
    for row in test1:
        if list(row).count(0) >= n - 1:
            return True

    # check columns
    for row in test2:
        if list(row).count(0) >= n - 1:
            return True
    return False


def _is_simple_brute(f):
    """Checks every combination of elementary monomials (x - x0) whether
    they have integer coefficients."""
    roots = f.roots()
    x = var("x")
    elementaries = []
    for x0, degree in roots:
        elementaries.extend((x - x0) for _ in range(degree))

    for subel in Subsets(range(len(elementaries))):
        if not subel:
            continue

        pol = 1
        for i in subel:
            pol *= elementaries[i]
        # print(pol)
        if all(coef.is_integer() for coef, _ in pol.coefficients()):
            return False
    return True


def is_simple(A: matrix):
    """Checks whether charpoly of A isn't divisible by a
    monic polynomial with integral entities and, as a result,
    whether A has invariant subgroup of Z^n.
    """
    if A.det() == 0:
        return False

    f = A.charpoly()

    if all(coef.is_integer() for coef in list(f)):
        return False

    if any(lamb.is_integer() for lamb, _ in f.roots()):
        return False

    if f.degree() > 2:
        return _is_simple_brute(f)
    return True


def from_indices_list(gens, triv, seq):
    res = copy(triv)
    for mul in seq:
        if mul:
            gen = gens[abs(mul) - 1]
            gen = gen if mul > 0 else gen.inverse()
            res *= gen
    return res


def gap_cryst_group(n, dim=3):
    return gap(f"SpaceGroupOnLeftIT({dim}, {n})")


# FIXME: wrong
def lattice_cosets(A):
    if A.det() == 0:
        raise NotImplementedError()

    D, U, V = A.smith_form()
    n = D.rank()
    gens = [list(range(D[i, i])) for i in range(n)]
    cosets = []
    for v in product(*gens):
        v = vector(QQ, v)
        cosets.append(V.inverse() * v)
    return cosets


class SpaceGroup_Element:
    def __init__(self, mtx):
        if isinstance(mtx, SpaceGroup_Element):
            mtx = copy(mtx._body)

        assert mtx.dimensions()[0] == mtx.dimensions()[1]
        self._body = matrix(mtx)

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
        return self._body.dimensions()[0] - 1

    def linear_part(self):
        n = self.dim
        return self._body[:n, :n]

    def translation(self):
        n = self.dim
        return self._body[:n, -1]

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
        return self.__class__(self._body**n)

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

    def __init__(
        self, generators: Iterable[SpaceGroup_Element], gap_G=None, dim=2, ita_num=None
    ):
        """Don't use the constructor explicitly. Instead use class method
        `from_gap_cryst` or `from_gens`.
        """
        self.gap_G = gap_G
        self.dim = dim
        self.ita_num = ita_num

        self.gap_P = self.gap_G.PointGroup()

        self.P_triv = matrix.identity(QQ, self.dim)
        self.G_triv = SpaceGroup_Element(matrix.identity(QQ, self.dim + 1))

        self.G_gens = copy(generators)

        self.G_nontriv = [el for el in self.G_gens if el.linear_part() != self.P_triv]
        self._P_names = [f"g_{i + 1}" for i in range(len(self.G_nontriv))]
        self._L_names = [f"e_{i + 1}" for i in range(dim)]

        # TODO: wrong assumption that all the elementary vectors would be in the generating set
        self.L_gens = [el for el in self.G_gens if el.linear_part() == self.P_triv]
        if len(self.L_gens) != dim:
            raise NotImplementedError(f"not enough generators: {len(self.L_gens)}")

        self._tr_basis = matrix(
            QQ, [el.translation().column(0) for el in self.L_gens]
        ).T

        self._gen2name = {}
        self._name2gen = {}
        self._rebuild_names()

        self.G_sorted_gens = self.G_nontriv + self.L_gens
        self.P_gens = [el.linear_part() for el in self.G_nontriv]

        self._P_dict = build_finite_group(self.P_gens, self.P_triv, self.max_iterations)

        self._alpha = {str(self.P_triv): self.G_triv.translation()}
        self.snot = [self.G_triv]
        for el, seq in self._P_dict.items():
            if el == str(self.P_triv):
                continue
            val = from_indices_list(self.G_nontriv, self.G_triv, seq)
            val = SpaceGroup_Element(val)
            assert str(val.linear_part()) == str(el), f"{val.linear_part()}, {el}"

            self._alpha[str(el)] = val.translation()
            self.snot.append(val)

        self._lattice_precalc = None

    def _rebuild_names(self):
        self._gen2name = {}
        self._name2gen = {}
        for name, el in zip(
            self._P_names + self._L_names, self.G_nontriv + self.L_gens
        ):
            self._name2gen[f"{name}^(-1)"] = el.inverse()
            self._name2gen[name] = el

            # firstly assert inverse elements. This guarantees that
            # elements of order 2 will have name x instead of x^(-1)
            self._gen2name[str(el.inverse())] = f"{name}^(-1)"
            self._gen2name[str(el)] = name

    def in_alpha(self, sym):
        return str(sym) in self._alpha

    def point_group_normalizer(self, **kwargs):
        P = MatrixGroup(self.P_gens)
        return normalizers(P, **kwargs)

    def in_lattice_basis(self):
        """Returns whether L == ZZ^n."""
        return self._tr_basis == self.P_triv

    def alpha(self, el):
        if isinstance(el, SpaceGroup_Element):
            p = el.linear_part()
        else:
            p = el
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
            el = self._change_basis([el._body], self._tr_basis)[0]
            return SpaceGroup_Element(el) in self.to_lattice_basis()

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
                res_word.append(f"{self._L_names[i]}^({x})")
            else:
                idx = len(self.P_gens) + i + 1
                res_word = [sign(x) * idx for _ in range(abs(x))] + res_word

        if not res_word:
            res_word = ["e"]

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

    def is_subgroup(self, H: SpaceGroup_gap):
        for el in self.G_sorted_gens:
            if el not in H:
                return False
        return True

    def index(self, H: SpaceGroup_gap):
        if not H.is_subgroup(self):
            raise ValueError("Not a subgroup.")

        return self.gap_G.Index(H.gap_G)

    def cosets(self, H: SpaceGroup_gap, action="left", lattice_only=False):
        if not H.is_subgroup(self):
            raise ValueError(f"Can't build cosets for {H}. Not a subgroup")

        if self.gap_G is None:
            raise NotImplementedError

        H_gap = self.gap_G.Subgroup([el._body for el in H.G_sorted_gens])

        trans = self.gap_G.RightTransversal(H_gap)
        trans_els = [SpaceGroup_Element(matrix(QQ, el)) for el in trans.AsList()]

        if action == "left":
            trans_els = [el.inverse() for el in trans_els]

        if lattice_only:
            trans_els = [self.rechoose_coset(el, H, action=action) for el in trans_els]

        return trans_els

    def rechoose_coset(self, g: SpaceGroup_Element, H: SpaceGroup_gap, action="left"):
        """Choose coset representative which is pure translation."""
        p = g.linear_part().inverse()
        t = H.alpha(p)
        g_inv = SpaceGroup_Element.construct_element(p, t)
        assert g_inv in H

        if action == "left":
            return g * g_inv
        else:
            return g_inv * g

    def find_coset(self, el, transversal, H, action="left"):
        for i, d in enumerate(transversal):
            if action == "left":
                triv = d.inverse() * el
            else:
                triv = el * d.inverse()
            if triv in H:
                return i

        raise ValueError(f"Element {el} is not in G.")

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

    def is_isomorphic(self, other: SpaceGroup_gap):
        # TODO: just go to lattice basis and check isomorphism between point groups
        #        and snot?

        raise NotImplementedError()

    def get_ITA(self):
        if self.ita_num is not None:
            return self.ita_num
        else:
            # TODO: go through every space group and check isomorphism
            raise NotImplementedError()

    def is_symmorphic(self):
        return self.gap_G.IsSymmorphicSpaceGroup()

    def is_simple_virtend(self, T):
        """Checks whether phi(g) = TgT^{-1} forms a simple surjective
        virtual endomorphism from a subgroup T^{-1}GT of finite index."""

        def phi_inv(_g):
            return T.inverse() * _g * T

        for el in self.G_sorted_gens:
            if phi_inv(el) not in self:
                return False

        p = SpaceGroup_Element(T).linear_part()

        return is_simple(p)

    def self_similar(self, T, verbose=False, quiet=False):
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
        quiet : bool
            if False, then function raises error if `T` doesn't generate
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
            print(
                "====================================================================="
            )
            print(f"sorted generators of group G #{self.get_ITA()}:")
            print(ascii_art(self.G_sorted_gens), end="\n\n")
            print("Virtual endomorphism is generated by T:")
            print(ascii_art(T), end="\n\n")

        if not self.is_simple_virtend(T):
            raise ValueError("Bad matrix T, no simple surjective virtual endomorphism")

        gens_H = [phi_inv(el)._body for el in self.G_sorted_gens]
        H = SpaceGroup_gap.from_gens(gens_H)

        if verbose:
            print("----------------------------------------------------")
            print("sorted generators of H:")
            print(ascii_art(gens_H), end="\n\n")
            print("Index of subgroup H:", self.index(H))

        trans_els = self.cosets(H, action="left", lattice_only=True)

        if verbose:
            print("Transversal:")
            print(ascii_art(trans_els), end="\n\n")

        # create self-similar action
        res_map = {}

        for a in self.G_sorted_gens:
            for i, d_i in enumerate(trans_els):
                adi = a * d_i
                name_a = self._gen2name[str(a)]

                if verbose:
                    print(f"{name_a}d_{i + 1}:")
                    print(adi, end="\n\n")

                j = self.find_coset(adi, trans_els, H, action="left")
                d_j = trans_els[j]

                if (d_j.inverse() * a * d_i) not in H:
                    raise ValueError(
                        f"This shouln't happen--wrong coset for triple \n{ascii_art([a, d_i, d_j])}, "
                        f"got: \n{ascii_art(d_j.inverse() * a * d_i)}"
                    )

                tmp_res = phi(d_j.inverse() * a * d_i)

                if verbose:
                    print(f"phi(d_{j + 1}^(-1) {name_a} d_{i + 1}):", tmp_res, sep="\n")

                res_map[(name_a, i + 1)] = (j + 1, tmp_res)

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
            new_gens = cls._change_basis(gens, matrix(QQ, G.TranslationBasis()).T)
            return cls(
                [SpaceGroup_Element(el) for el in new_gens],
                gap.AffineCrystGroupOnLeft(new_gens),
                dim,
                ita_num,
            )
        else:
            return cls(
                [
                    SpaceGroup_Element.from_gap_element(el)
                    for el in G.GeneratorsOfGroup()
                ],
                G,
                dim,
                ita_num,
            )

    @classmethod
    def from_gens(cls, generators: Iterable[matrix]):
        gens = [SpaceGroup_Element(el) for el in generators]
        dim = gens[0].dim
        return cls(gens, gap.AffineCrystGroupOnLeft(generators), dim, None)

    @staticmethod
    def _change_basis(gens, basis: matrix):
        dim = basis.dimensions()[0]
        trans = matrix(QQ, [0 for _ in range(dim)]).T
        conj = block_matrix(QQ, [[basis, trans], [0, 1]])
        new_gens = [conj.inverse() * el * conj for el in gens]
        return new_gens

    def change_basis(self, basis: matrix):
        new_gens = self._change_basis([el._body for el in self.G_sorted_gens], basis)
        # print(ascii_art(new_gens))
        return self.__class__.from_gens(new_gens)

    def naive_self_replicating(self, verbose=False, build_wreath=False):
        denominators = set()
        for el in self._alpha.values():
            for row in el:
                denominators.add(row[0].denominator())

        m = lcm(denominators)

        T = (1 / (m + 1)) * matrix.identity(QQ, self.dim)
        T = SpaceGroup_Element.from_linear(T)

        action = self.self_similar(T, verbose=verbose, quiet=False)
        if build_wreath:
            self._wreath_recursion(action)

        return action

    def minimal_self_replicating(self, verbose=False, build_wreath=False):
        if self.dim != 2:
            raise NotImplementedError()

        T = _PLANAR_GROUPS_MIN_SR[self.get_ITA()]

        action = self.self_similar(T, verbose=verbose, quiet=False)

        if build_wreath:
            self._wreath_recursion(action)
        return action

    def to_lattice_basis(self):
        if self.in_lattice_basis():
            return self

        if self._lattice_precalc is None:
            self._lattice_precalc = self.change_basis(self._tr_basis)
        return self._lattice_precalc

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

            w_el[el][let - 1] = "".join(self.as_word(v[1], readable=True))
            w_p[el][let - 1] = v[0]

        print("Вінцева рекурсія:")
        for el in w_el.keys():
            p = Permutation(w_p[el])
            if p.to_permutation_group_element().order() == 1:
                p = "()"
            else:
                p = str(p.cycle_tuples())

            print(
                f"$${el} = {p}{w_el[el]}$$".replace("'", "")
                .replace("[", "(")
                .replace("]", ")")
                .replace("(-1)", "{-1}")
                .replace("^(1)", "")
            )

        print("```")
        print(self._name2gen)
        print("```")


# fmt: off
_PLANAR_GROUPS_MIN_SR = {
    1: [
        [0, 1/2, 0], 
        [1, 0, 0], 
        [0, 0, 1]
    ], 

    2: [
        [0, 1/2, 0], 
        [1, 0, 0], 
        [0, 0, 1]
    ], 
    
    3: [
        [1/2, 0, 0], 
        [0, 1/2, 0], 
        [0, 0, 1]
    ],

    4: [
        [1/2, 0, 0], 
        [0, 1/3, 0], 
        [0, 0, 1]
    ],

    5: [
        [1/2, 0, 0],
        [0, 1/2, 0],
        [0, 0, 1]
    ],

    6: [
        [0, 1/2, 0],
        [1, 0, 0],
        [0, 0, 1],
    ],

    7: [
        [1/3, 0, 0], 
        [0, 1/2, 0], 
        [0, 0, 1],
    ], 

    8: [
        [0, 1/3, 0], 
        [1, 0, 0], 
        [0, 0, 1],
    ], 

    9: [
        [1/3, 2/3, 0], 
        [-2/3, -1/3, 0], 
        [0, 0, 1],
    ],

    10: [
        [1/2, 1/2, 0], 
        [-1/2, 1/2, 0], 
        [0, 0, 1],
    ], 

    11: [
        [1/2, 1/2, 0], 
        [-1/2, 1/2, 0], 
        [0, 0, 1],
    ], 

    12: [
        [-1/3, 0, 0],
        [0, 1/3, 0],
        [0, 0, 1],
    ], 

    13: [
        [-1/3, 2/3, 0,],
        [1/3, 1/3, 0], 
        [0, 0, 1],
    ], 

    14: [
        [1/2, 0, 0], 
        [0, 1/2, 0], 
        [0, 0, 1]
    ],

    15: [
        [1/2, 0, 0], 
        [0, 1/2, 0], 
        [0, 0, 1]
    ],

    16: [
        [-1/3, 2/3, 0,],
        [1/3, 1/3, 0], 
        [0, 0, 1],
    ], 

    17: [
        [2/3, -1/3, 0], 
        [1/3, 1/3, 0], 
        [0, 0, 1],
    ]
}

for k, el in _PLANAR_GROUPS_MIN_SR.items(): 
    _PLANAR_GROUPS_MIN_SR[k] = matrix(QQ, el)
