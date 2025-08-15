#!/usr/bin/env sage

from collections import defaultdict
from warnings import warn

from sage.all import QQ, matrix, solve, var


def normalize_expressions(exps, allowed=None):
    """Replace all variables except allowed (default: a b c d) with {x0,x1,x2,x3...}"""
    extra_args = set()

    if allowed is None:
        allowed = var('a b c d')

    for exp in exps:
        extra_args |= {el for el in exp.args() if el not in allowed}

    res = []
    for exp in exps:
        for i, el in enumerate(extra_args):
            exp = exp.subs({var(el): var(f'x{i}')})
        res.append(exp)
    return tuple(res)


alphabet = 'abcdefghijklmnopqrstuvwxyz'


def create_symbolic_matrix(dim, use_alphabet=False):
    """Create symbolic matrix of given dimention.

    Parameters
    ----------
    dim : integer
    use_alphabet
        if True, then symbolic matrices are created like
                        [a, b]
                        [c, d]
        if False, then
                    [a_00, a_01]
                    [a_10, a_11]

    Returns
    -------
    Symbolic matrix
    """
    if use_alphabet and dim * dim > len(alphabet):
        raise ValueError(f"Can't use alphabet for matrix {dim}x{dim} due to lack of letters.")
    A = []
    args = []
    for i in range(dim):
        row = []
        for j in range(dim):
            if use_alphabet:
                row.append(var(alphabet[i * dim + j]))
            else:
                row.append(var(f'a_{i}{j}'))
            args.append(row[-1])
        A.append(row)
    A = matrix(A)
    return A, args


def sol2matrix(solution, dim=2, use_alphabet=False):
    """Transform solution into symbolic matrix."""
    if use_alphabet and dim * dim > len(alphabet):
        raise ValueError(f"Can't use alphabet for matrix {dim}x{dim} due to lack of letters.")
    A = []
    sol_dict = {sol.left(): sol.right() for sol in solution}

    for i in range(dim):
        row = []
        for j in range(dim):
            if use_alphabet:
                row.append(sol_dict[var(alphabet[i * dim + j])])
            else:
                row.append(sol_dict[var(f'a_{i}{j}')])

        A.append(row)
    A = matrix(A)
    return A


def els_by_order(group_elements):
    res = defaultdict(list)

    for el in group_elements:
        res[el.order()].append(el)
    return res


def gens_mappings(gens, group_elements):
    orders = els_by_order(group_elements)
    res = []
    for gen in gens:
        res.append(orders[gen.order()])
    return res


def carthesian_wo_duplicates(*spaces):
    for res, _ in _carthesian_wo_duplicates(*spaces):
        yield res


def _carthesian_wo_duplicates(*spaces):
    if not spaces:
        yield [], set()
        return

    for res, used in _carthesian_wo_duplicates(*spaces[:-1]):
        for el in spaces[-1]:
            if str(el) in used:
                continue

            used.add(str(el))
            yield res + [el], used
            used.remove(str(el))


def normalizers(P, verbose=False, use_alphabet=False, normalize_exp=True, to_matrix=True, ignore_trivial=True):
    """Find normalizer of the PointGroup in GL(n, QQ).

    Tries to find normalizer as a solution of
    system of linear equations N*A_i = A_j*N, where A_i and A_j are
    two elements of the PointGroup with respect to a permutation. Since
    Normalizer is a group N, that satisfies condition NA = AN,
    we can check all the automorphisms of A and solve the respective system
    of linear equations.

    Parameters
    ----------
    P : MatrixGroup
        the respective PointGroup

    use_alphabet : bool
        if True, then symbolic matrices are created like
                                [a, b]
                                [c, d]
        if False, then
                            [a_00, a_01]
                            [a_10, a_11]
    normalize_exp : bool
        if True, then solution of the linear system will be
        normalized, i.e. all the independed variables will be
        renamed into x_0, x_1, x_2 ...
    to_matrix : bool
        if True, then the result solutions will be transformed
        into symbolic matrices instead of tuple of Expressions
    verbose : bool
        True to see the results on the fly.
    ignore_trivial : bool
        if True then solutions with zero determinant will be ignored

    Returns
    -------
    If `to_matrix` is True, then list of symbolic matrices is returned.
    Else, list of tuples of expressions like "a_00 == x1" which denotes
         what elements of matrix should be
    """
    P_elements = list(P)
    dim = matrix(QQ, P_elements[0]).dimensions()[0]

    possible_mappings = gens_mappings(P.gens(), P_elements)

    if verbose:
        print('\n====================================================')
        print('point group:', P)

        print('group elements:')
        print(*P_elements, sep='\n')
        print('\n----------------normalizers-------------------------')

    A, args = create_symbolic_matrix(dim, use_alphabet=use_alphabet)
    found_solutions = set()

    for maps_to in carthesian_wo_duplicates(*possible_mappings):
        # if len(set(str(el) for el in maps_to)) < len(maps_to):
        #     continue

        # build conditions AX = YA for every X -> Y due to chosen permutation
        cond_perm = set()
        try:
            homm = P.hom(maps_to)
        except ValueError:
            warn(f'Something went wrong: couldnt create homomorphism for {maps_to}')
            continue

        # we need only automorphisms
        if not homm.kernel().is_trivial():
            continue

        for X in P_elements:
            Y = homm(X)

            X, Y = matrix(QQ, X), matrix(QQ, Y)
            cond_i = A * X - Y * A
            for el in cond_i:
                cond_perm = cond_perm.union(set(el))

        eq = [cond == 0 for cond in cond_perm]
        for res in solve(eq, *args):
            if not res:
                if verbose:
                    print(res)
                continue

            res = tuple(res)
            if normalize_exp:
                res = normalize_expressions(res, allowed=args)

            mtx = sol2matrix(res, dim=dim, use_alphabet=use_alphabet)

            if res not in found_solutions and verbose:
                if to_matrix:
                    if ignore_trivial and mtx.det() == 0:
                        continue
                    print(mtx, end='\n\n')
                else:
                    print(res)

            if mtx.det() != 0 or not ignore_trivial:
                found_solutions.add(res)

    if to_matrix:
        return [sol2matrix(solution, dim=dim, use_alphabet=use_alphabet) for solution in found_solutions]
    else:
        return found_solutions
