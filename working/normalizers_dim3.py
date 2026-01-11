from sage.all import MatrixGroup, var, matrix, Integer
from src.space_groups import SpaceGroup_gap, prepare_gap_env, check_div

# from src.srdegrees import _factorize

prepare_gap_env()

# __f = open('normalizers_dim3_changed_basis', 'w')


def tau(A):
    row1, row2, row3 = A
    return _tau(*row1, *row2, *row3)


def _tau(y0, y1, y2, y3, y4, y5, y6, y7, y8):
    return y0 * y4 + y0 * y8 + y4 * y8 - y1 * y3 - y2 * y6 - y5 * y7


def nprint(*args, **kwargs):
    print(*args, **kwargs)
    # print(*args, **kwargs, file=__f)


nontrivial = [16, 20, 22, 23, 47, 63, 69, 71, 195, 196, 197, 202, 204, 207, 209, 211, 216, 217, 229]


def rawdog():
    found = set()
    point_groups = set()

    for i in nontrivial:
        nprint(f'##################### {i} ############################')
        G = SpaceGroup_gap.from_gap_cryst(i, dim=3, change_basis=True)
        els = list(sorted(MatrixGroup(G.P_gens)))
        if str(els) in point_groups:
            continue
        point_groups.add(str(els))
        for A in G.point_group_normalizer():
            if tau(A) == 0:
                continue
            if str(A) in found:
                continue
            if check_div(A):
                continue
            found.add(str(A))
            nprint(A)
            nprint('tau(A) = ', tau(A).simplify_rational())
            nprint('\n')


# __f.close()

half = Integer(1) / Integer(2)
quar = Integer(1) / Integer(4)
x0, x1, x2, x3 = var('x0 x1 x2 x3')

# fmt: off
cached_matrices = {

    16: [
        [
            [0, 0, x2],
            [x1, 0, 0],
            [0, x0, 0],
        ],
        [
            [0, x0, 0],
            [0, 0, x1],
            [x2, 0,  0],
            ],
    ],

    20: [
        [
            [0, 0, -2*x2],
            [x1, 0, x2],
            [half*x0, x0, 0],
        ],
        [
            [half*x1, x1, 0],
            [-quar*x1, -half*x1, x2],
            [x0, 0, 0],]
        ],

    22: [
        [
            [-x2, -x2, -2 * x2],
            [-2 * x1 + x2, 0, 0],
            [x1, x0, x2],
        ],
        [
            [0, x0 - 2 * x1, 0],
            [-x0, -x0, -2 * x0],
            [x2, x1, x0],
        ],
    ],

    23: [
        [
            [-x0, 0, -2 * x0],
            [x1, 0, x0],
            [half * x0 + half * x2, x2, x0],
        ],
        [
            [half * x1, x1, 0],
            [x2, -half * x1, half * x1 + 2 * x2],
            [x0, -half * x1, 0],
        ],
    ],

    47: [
        [
            [0, x1, 0],
            [0, 0, x2],
            [x0, 0, 0],
        ],
    ],

    63: [
        [
            [0, 0, -2 * x1],
            [x0, 0, x1],
            [half * x2, x2, 0],
        ],
        [
            [x0, 2 * x0, 0],
            [-half * x0, -x0, x1],
            [x2, 0, 0],
        ],
    ],

    69: [
        [
            [-x0, -x0, -2 * x0],
            [x0 - 2 * x2, 0, 0],
            [x2, x1, x0],
        ],
    ],

    71: [
        [
            [x1, 2 * x1, 0],
            [x2, -x1, x1 + 2 * x2],
            [x0, -x1, 0],
        ],
        [
            [-2 * x1 + x2, 0, -4 * x1 + 2 * x2],
            [x0, 0, 2 * x1 - x2],
            [x1, x2, 2 * x1 - x2],
        ],
    ],

    195: [
        [
            [0, x0, 0],
            [0, 0, -x0],
            [-x0, 0, 0],
        ],
        [
            [0, 0, x0],
            [x0, 0, 0],
            [0, -x0, 0],
        ],
        [
            [0, 0, x0],
            [-x0, 0, 0],
            [0, x0, 0],
        ],
        [
            [0, x0, 0],
            [0, 0, x0],
            [x0, 0, 0],
        ],
        [
            [0, x0, 0],
            [0, 0, x0],
            [-x0, 0, 0],
        ],
        [
            [0, 0, x0],
            [x0, 0, 0],
            [0, x0, 0],
        ],
        [
            [0, x0, 0],
            [0, 0, -x0],
            [x0, 0, 0],
        ],
        [
            [0, 0, x0],
            [-x0, 0, 0],
            [0, -x0, 0],
        ],
    ],

    196: [
        [
            [0, x0, 0],
            [x0, x0, 2 * x0],
            [-x0, -x0, -x0],
        ],
        [
            [-x0, -x0, -2 * x0],
            [-x0, 0, 0],
            [x0, x0, x0],
        ],
        [
            [half * x0, half * x0, x0],
            [half * x0, 0, 0],
            [-half * x0, 0, -half * x0],
        ],
        [
            [0, -half * x0, 0],
            [half * x0, half * x0, x0],
            [-half * x0, 0, -half * x0],
        ],
        [
            [-x0, -x0, -2 * x0],
            [x0, 0, 0],
            [0, x0, x0],
        ],
        [
            [0, x0, 0],
            [x0, x0, 2 * x0],
            [0, -x0, -x0],
        ],
    ],

    197: [
        [
            [-x0, -2 * x0, 0],
            [x0, x0, x0],
            [x0, x0, 0],
        ],
        [
            [-x0, 0, -2 * x0],
            [0, 0, x0],
            [x0, x0, x0],
        ],
        [
            [-x0, 0, -2 * x0],
            [x0, 0, x0],
            [x0, x0, x0],
        ],
        [
            [-x0, -2 * x0, 0],
            [x0, x0, x0],
            [0, x0, 0],
        ],
        [
            [x0, 2 * x0, 0],
            [0, -x0, x0],
            [-x0, -x0, 0],
        ],
        [
            [-x0, 0, -2 * x0],
            [x0, 0, x0],
            [0, -x0, x0]
        ],
    ],

    198: [
        [
            [x0, x0, 2 * x0],
            [-x0, 0, 0],
            [0, -x0, -x0],
        ],
        [
            [0, x0, 0],
            [-x0, -x0, -2 * x0],
            [x0, 0, x0],
        ],
        [
            [0, -x0, 0],
            [-x0, -x0, -2 * x0],
            [x0, x0, x0],
        ],
        [
            [-x0, -x0, -2 * x0],
            [-x0, 0, 0],
            [x0, 0, x0],
        ],
    ],

    203: [
        [
            [x0, 2 * x0, 0],
            [-x0, -x0, -x0],
            [-x0, -x0, 0],
        ],
    ],

    205: [
        [
            [0, 0, -x0],
            [x0, 0, 0],
            [0, -x0, 0],
        ],
        [
            [0, 0, -x0],
            [x0, 0, 0],
            [0, x0, 0],
        ],
    ],

    208: [
        [
            [x0, x0, 2 * x0],
            [x0, 0, 0],
            [-x0, -x0, -x0],
        ],
        [
            [0, -x0, 0],
            [x0, x0, 2 * x0],
            [-x0, 0, -x0],
        ],
    ],

    210: [
        [
            [x0, 0, 2 * x0],
            [0, 0, -x0],
            [-x0, -x0, -x0],
        ],
    ],

    212: [
        [
            [x0, x0, 2 * x0],
            [x0, 0, 0],
            [-x0, 0, -x0],
        ],
    ],

    217: [
        [
            [-x0, -2 * x0, 0],
            [0, x0, -x0],
            [x0, x0, 0],
        ],
        [
            [x0, 0, 2 * x0],
            [-x0, 0, -x0],
            [-x0, -x0, -x0],
        ],
    ],

    229: [
        [
            [half * x0, 0, x0],
            [-half * x0, 0, -half * x0],
            [0, half * x0, -half * x0],
        ],
        [
            [half * x0, x0, 0],
            [-half * x0, -half * x0, -half * x0],
            [0, -half * x0, 0],
        ],
    ],
}

for gr_num in cached_matrices:
    cached_matrices[gr_num] = [matrix(el) for el in cached_matrices[gr_num]]


for gr_num in cached_matrices:
    print(f'################# {gr_num} ######################')
    for a in cached_matrices[gr_num]:
        # if tau(a) == 0:
        #     continue
        # print(a)
        print("equation:", (a.det() - a.trace()).simplify_rational())
        # print('tau(A):', tau(a))
