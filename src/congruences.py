
from sage.all import (
    column_matrix,
    QQ,
    lcm,
    matrix,
    var,
    ZZ,
)

from dataclasses import dataclass


from .meta import Option


@dataclass
class Solution_QZ_rational:
    part_answer: matrix
    lattice_basis: matrix
    nullbasis: matrix


@dataclass
class Solution_QZ_integral:
    base_variables: list
    free_variables: list
    expressions: list


def solve_qz_integral(M, b, base_variables) -> Option[Solution_QZ_integral]:
    # in most cases there is a less common denominator
    m = lcm(M.denominator(), b.denominator())

    N_conds = M.nrows()

    M_lift = matrix(ZZ, M * m)
    b_lift = matrix(ZZ, b * m)

    E = matrix.identity(ZZ, N_conds)

    # transform system Mx = b mod Z into equivalent -z + Mx = b for integral z
    # then lift it to integers by multiplying on lcm of denominators:   -mz + mMx = mb
    M_aug = (-m * E).augment(M_lift)

    try:
        # one solution for inhomogeneous system
        base_vec = M_aug.solve_right(b_lift, extend=False)
    except ValueError:
        return Option[Solution_QZ_integral].error("can't solve over integers: no part answer")

    # solutions space for homogeneous system. The result solution will be base_vec + lattice
    # Notes:
    #     right_kernel returns echelon form of the basis vectors. We placed auxiliary variables at the beginning of
    #     the list of variables. Therefore, the first r columns of the lattice will be the standard e1 e2 ... er vectors
    #     and the respective first r variables will be free, which are all auxiliary variables + free base variables.
    #
    #     By multiplying variables*lattice we get linear combinations for all the dependent variables, which will be
    #     our answer.
    lattice = M_aug.right_kernel(algorithm='pari').basis_matrix()
    aug_variables = [var(f'y{i}') for i in range(N_conds)]
    variables = aug_variables + list(base_variables)

    # exploit echelon form
    free_n = lattice.nrows()
    free_vars = variables[:free_n]
    free_vars = matrix(free_vars)

    # get only base variables, don't use auxiliary
    result = (free_vars * lattice + base_vec.T)
    result = result[0][N_conds:]

    return Option[Solution_QZ_rational].success(
        Solution_QZ_integral(
            base_variables,
            free_vars,
            result,
        )
    )


def solve_qz_rational(M, b) -> Option[Solution_QZ_rational]:
    # in most cases there is a less common denominator
    m = lcm(M.denominator(), b.denominator())

    N_conds, N_vars = M.nrows(), M.ncols()

    M_lift = matrix(ZZ, M * m)
    b_lift = matrix(b * m)

    S, U, V = M_lift.smith_form()

    rank = 0
    for i in range(min(N_conds, N_vars)):
        if S[i, i] != 0:
            rank = i+1
    right = U * b_lift
    for i in range(right.dimensions()[0]):
        for j in range(right.dimensions()[1]):
            right[i, j] = right[i, j] % m

    part_answer = matrix(QQ, [0 for i in range(N_vars)]).T
    lattice_basis = []
    null_basis = []

    # every solution for the dx = b  mod M  has the form
    # x = b/d + nM/d, where n is arbitrary integral
    for i in range(rank):
        if right[i] != 0:
            part_answer[i, 0] = right[i, 0] / S[i, i]

        # the lattice nM/d
        zv = [0 for i in range(N_vars)]
        zv[i] = m / S[i, i]
        lattice_basis.append(zv)

    for i in range(rank, N_conds):
        if right[i] != 0:
            return Option[Solution_QZ_rational].error(f'no solutions: {right[i]} != 0')

    # the null basis consists of all the arbitrary rational vectors which
    # can be placed in the equations of kind 0x = 0
    for i in range(rank, N_vars):
        rv = [0 for i in range(N_vars)]
        rv[i] = 1
        null_basis.append(rv)

    if lattice_basis:
        lattice_basis = V * column_matrix(QQ, lattice_basis)
    if null_basis:
        null_basis = V * column_matrix(QQ, null_basis)

    return Option[Solution_QZ_rational].success(
        Solution_QZ_rational(
            V * part_answer,
            lattice_basis,
            null_basis,
        )
    )
