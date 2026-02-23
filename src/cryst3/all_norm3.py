from sage.all import Integer, var

int0 = Integer(0)
int1 = Integer(1)
int2 = Integer(2)
int3 = Integer(3)
int4 = Integer(4)

x0, x1, x2, x3, x4, x5 = var('x0 x1 x2 x3 x4 x5')
all_matrices = [
    [
        [x1, int0, x2],
        [int0, x4, int0],
        [x3, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [int0, -int2 * x0, int0],
        [x0, int2 * x0, int0],
        [int1 / int2 * x1, x0, x1],
    ],
    [
        [-x1, -int2 * x1, int0],
        [x1, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, x1, x0],
    ],
    [
        [int0, int2 * x1, int0],
        [x1, int0, int0],
        [int1 / int2 * x0, -x1, x0],
    ],
    [
        [-x1, -int2 * x1, int0],
        [int0, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, x1, x0],
    ],
    [
        [-int2 * x1, -int2 * x1, int0],
        [x1, int2 * x1, int0],
        [int1 / int2 * x0 + x1, x1, x0],
    ],
    [
        [-x1, int0, int0],
        [x1, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, int0, x0],
    ],
    [
        [-int2 * x1, -int2 * x1, int0],
        [x1, int0, int0],
        [int1 / int2 * x0 + x1, x1, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int2 * x0, int0],
        [int0, -x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, -x0, x1],
    ],
    [
        [-x1, int0, int0],
        [x1, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [-x1, -int2 * x1, int0],
        [x1, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, x1, x0],
    ],
    [
        [x0, int0, -int2 * x1],
        [-int1 / int2 * x0 + int1 / int2 * x3, x3, x1],
        [x2, int0, x4],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int0],
        [
            -int1 / int4 * x0 + int1 / int2 * x1,
            -int1 / int2 * x0,
            x1,
        ],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [x0, int0, int0],
        [-x0, -x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [-x1, -int2 * x1, int0],
        [x1, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, x1, x0],
    ],
    [
        [-x0, x0, int0],
        [x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, -x0, int0],
        [x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, x0, int0],
        [-x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, x0, int0],
        [x0, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int0],
        [
            -int1 / int4 * x0 + int1 / int2 * x1,
            -int1 / int2 * x0,
            x1,
        ],
    ],
    [
        [int0, -int2 * x1, int0],
        [x1, int2 * x1, int0],
        [int1 / int2 * x0, x1, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, x0, x1],
    ],
    [
        [int0, int2 * x1, int0],
        [x1, int0, int0],
        [int1 / int2 * x0, -x1, x0],
    ],
    [
        [x0, int0, int0],
        [-x0, -x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [-int2 * x1, -int2 * x1, int0],
        [x1, int2 * x1, int0],
        [int1 / int2 * x0 + x1, x1, x0],
    ],
    [
        [-int2 * x1, -int2 * x1, int0],
        [x1, int0, int0],
        [int1 / int2 * x0 + x1, x1, x0],
    ],
    [
        [x1, x0 - x1, int0],
        [-x0 + x1, x0, int0],
        [int0, int0, x2],
    ],
    [
        [-x2, x1 + x2, int0],
        [x1, x2, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int3 * x0 + int3 * x2, int0],
        [x2, x0, int0],
        [
            int2 / int3 * x0 + int2 / int3 * x1,
            int2 * x0 - int2 * x2,
            x1,
        ],
    ],
    [
        [x0 - int3 / int2 * x1, x0 - int3 / int2 * x1 - x2, int0],
        [
            -int1 / int3 * x0
            + int1 / int2 * x1
            + int1 / int3 * x2,
            x2,
            int0,
        ],
        [x1, -int2 / int3 * x0 + x1 + int2 / int3 * x2, x0],
    ],
    [
        [x0, -x0 + x2, int0],
        [x0 - x2, x2, int0],
        [int0, int0, x1],
    ],
    [
        [x2, x1 - x2, int0],
        [x1, -x2, int0],
        [int0, int0, x0],
    ],
    [
        [-x2, int3 * x0 - int3 * x2, int0],
        [x0, x2, int0],
        [
            int2 / int3 * x1 + int2 / int3 * x2,
            -int2 * x0 + int2 * x2,
            x1,
        ],
    ],
    [
        [x1 - int3 / int2 * x2, -x0 + x1 - int3 / int2 * x2, int0],
        [
            int1 / int3 * x0
            - int1 / int3 * x1
            + int1 / int2 * x2,
            x0,
            int0,
        ],
        [x2, int2 / int3 * x0 - int2 / int3 * x1 + x2, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [x0, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, -x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, x0, int0],
        [int0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, -x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-int1 / int2 * x1, -int3 / int2 * x1, int0],
        [int0, int1 / int2 * x1, int0],
        [int2 / int3 * x0 + int1 / int3 * x1, x1, x0],
    ],
    [
        [-int2 * x1, -int3 * x1, int0],
        [x1, int2 * x1, int0],
        [int2 / int3 * x0 + int4 / int3 * x1, int2 * x1, x0],
    ],
    [
        [-x1, -int3 * x1, int0],
        [x1, int2 * x1, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, int2 * x1, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int2 / int3 * x0 - int2 / int3 * x1, int0, x0],
    ],
    [
        [-x1, int0, int0],
        [x1, x1, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, int0, x0],
    ],
    [
        [-x0, -int3 / int2 * x0, int0],
        [int1 / int2 * x0, int1 / int2 * x0, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, x0, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [x0, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, -x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, int0, x0],
        [int0, x2, int0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x1],
        [x2, int0, int0],
    ],
    [
        [int0, x0, int0],
        [x2, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, int0, int0],
        [int0, x1, int0],
        [int0, int0, x2],
    ],
    [
        [int0, int0, x1],
        [x0, int0, int0],
        [int0, x2, int0],
    ],
    [
        [int0, int0, x0],
        [int0, x2, int0],
        [x1, int0, int0],
    ],
    [
        [-int1 / int2 * x1, -int3 / int2 * x1, int0],
        [int0, int1 / int2 * x1, int0],
        [int2 / int3 * x0 + int1 / int3 * x1, x1, x0],
    ],
    [
        [-int1 / int2 * x1, -int3 / int2 * x1, int0],
        [int1 / int2 * x1, x1, int0],
        [int2 / int3 * x0 + int1 / int3 * x1, x1, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-int2 / int3 * x0 + int2 / int3 * x1, int0, x1],
    ],
    [
        [-x1, -int3 / int2 * x1, int0],
        [int1 / int2 * x1, int1 / int2 * x1, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, x1, x0],
    ],
    [
        [-x0, -int3 / int2 * x0, int0],
        [int1 / int2 * x0, x0, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, x0, x1],
    ],
    [
        [-x1, int0, int0],
        [x1, x1, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, int0, x0],
    ],
    [
        [-x0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, x1, int0],
        [-x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [x0, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-int2 * x0, -int3 * x0, int0],
        [x0, x0, int0],
        [int4 / int3 * x0 + int2 / int3 * x1, int2 * x0, x1],
    ],
    [
        [-int1 / int2 * x1, -int3 / int2 * x1, int0],
        [int0, int1 / int2 * x1, int0],
        [int2 / int3 * x0 + int1 / int3 * x1, x1, x0],
    ],
    [
        [-x1, -int3 / int2 * x1, int0],
        [int1 / int2 * x1, x1, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, x1, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-int2 / int3 * x0 + int2 / int3 * x1, int0, x1],
    ],
    [
        [-x1, -int3 * x1, int0],
        [x1, int2 * x1, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, int2 * x1, x0],
    ],
    [
        [-x1, int0, int0],
        [x1, x1, int0],
        [int2 / int3 * x0 + int2 / int3 * x1, int0, x0],
    ],
    [
        [x1, x0, int0],
        [-x0, x0 + x1, int0],
        [int0, int0, x2],
    ],
    [
        [x2, x1, int0],
        [x1 + x2, -x2, int0],
        [int0, int0, x0],
    ],
    [
        [x1, x0 - x1, int0],
        [x0, -x1, int0],
        [int0, int0, x2],
    ],
    [
        [x2, -x1, int0],
        [x1, -x1 + x2, int0],
        [int0, int0, x0],
    ],
    [
        [x2, x1, int0],
        [-x1, x1 + x2, int0],
        [int0, int0, x0],
    ],
    [
        [x2, x1, int0],
        [x1 + x2, -x2, int0],
        [int0, int0, x0],
    ],
    [
        [x1, x1, int0],
        [-x1, int2 * x1, int0],
        [int0, int0, x0],
    ],
    [
        [-int1 / int2 * x0, x0, int0],
        [-x0, int1 / int2 * x0, int0],
        [int0, int0, x1],
    ],
    [
        [-x1, -x1, int0],
        [-int2 * x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, int0, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-int2 * x0, x0, int0],
        [-x0, int2 * x0, int0],
        [int0, int0, x1],
    ],
    [
        [-int2 * x0, x0, int0],
        [-x0, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, -x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [-int1 / int2 * x0, x0, int0],
        [int1 / int2 * x0, int1 / int2 * x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [-int2 * x1, x1, int0],
        [-x1, int2 * x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, x1, int0],
        [int2 * x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, -x1, int0],
        [x1, -int2 * x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [x0, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [int2 * x0, -x0, int0],
        [x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, -x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, int2 * x0, int0],
        [x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [int1 / int2 * x1, -x1, int0],
        [x1, -int1 / int2 * x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, x0, int0],
        [int0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, -x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, x1, int0],
        [-x1, int2 * x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -x0, int0],
        [-int2 * x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [-x1, x1, int0],
        [-x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [-int2 * x1, x1, int0],
        [-x1, int2 * x1, int0],
        [int0, int0, x0],
    ],
    [
        [-int1 / int2 * x0, x0, int0],
        [-x0, int1 / int2 * x0, int0],
        [int0, int0, x1],
    ],
    [
        [-x1, int0, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, x1, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [x1, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x1, int0],
        [-x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [-int2 * x1, x1, int0],
        [-x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-int1 / int2 * x0, x0, int0],
        [int1 / int2 * x0, int1 / int2 * x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [x0, int0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, -x0],
        [int0, x0, int0],
        [x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, -x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, -x0],
        [x0, int0, int0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, -x0],
        [int0, x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, int0, -x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, int0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [x0, int0, x0],
    ],
    [
        [x0, x0, int2 * x0],
        [-x0, int0, int0],
        [int0, -x0, -x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [-x0, int0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [-x0, int0, int0],
        [x0, x0, x0],
    ],
    [
        [int1 / int2 * x0, int0, int0],
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [-int1 / int2 * x0, int0, -int1 / int2 * x0],
    ],
    [
        [x0, int0, int0],
        [int0, -x0, int0],
        [int0, x0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, x0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [x0, x0, int2 * x0],
        [-x0, int0, int0],
        [int0, int0, -x0],
    ],
    [
        [int0, x0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [x0, x0, int2 * x0],
        [int0, -x0, int0],
        [int0, int0, -x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, x0, x0],
    ],
    [
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [int0, int1 / int2 * x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [int0, x0, int0],
        [x0, x0, int2 * x0],
        [-x0, -x0, -x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [int0, int1 / int2 * x0, int0],
        [int0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, x0],
        [x0, x0, int0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, int0],
        [int0, x0, -x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, int0, x0],
        [int0, -x0, x0],
    ],
    [
        [x0, int0, int0],
        [-x0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [int0, x0, -x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, x0, x0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [x0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [x0, x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, int0, x0],
        [int0, -x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, x0, x0],
        [x0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [int0, int0, -int1 / int2 * x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [int0, -int1 / int2 * x0, int0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
        [int0, -int1 / int2 * x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [-x0, -x0, int0],
    ],
    [
        [x0, int0, int0],
        [x1, x0 + int2 * x1, int0],
        [int0, int0, x2],
    ],
    [
        [x0, int0, int0],
        [-int1 / int2 * x0, int0, x2],
        [int1 / int2 * x1, x1, int0],
    ],
    [
        [x1, int2 * x1, int0],
        [x2, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x1, x1, int0],
        [-int1 / int4 * x1, -int1 / int2 * x1, x2],
        [x0, int0, int0],
    ],
    [
        [int0, int0, x2],
        [x1, int2 * x1, -int1 / int2 * x2],
        [x0, int0, int0],
    ],
    [
        [int0, int0, -int2 * x1],
        [x0, int0, x1],
        [int1 / int2 * x2, x2, int0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [x0, int0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, -x0],
        [int0, x0, int0],
        [x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, -x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, -x0],
        [x0, int0, int0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, -x0],
        [int0, x0, int0],
        [-x0, int0, int0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, -x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [int1 / int2 * x0, int0, int0],
        [-int1 / int2 * x0, int0, -int1 / int2 * x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, int0, int0],
        [x0, x0, x0],
    ],
    [
        [int0, -x0, int0],
        [x0, x0, int2 * x0],
        [int0, int0, -x0],
    ],
    [
        [-x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [-x0, int0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, -x0, int0],
        [int0, x0, x0],
    ],
    [
        [x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, -x0, int0],
        [int0, x0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [x0, int0, int0],
        [int0, x0, x0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, x0, x0],
    ],
    [
        [x0, int0, int0],
        [x0, x0, int2 * x0],
        [-x0, -x0, -x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, -x0, -int2 * x0],
        [int0, x0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, x0, int2 * x0],
        [-x0, -x0, -x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [int0, int0, -int1 / int2 * x0],
        [int0, int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [-int1 / int2 * x0, int0, -int1 / int2 * x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [x0, int0, int0],
        [-x0, int0, -x0],
        [int0, x0, int0],
    ],
    [
        [x0, int0, int0],
        [-x0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [int0, x0, -x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, int0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, -x0],
        [x0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [x0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [int0, int0, -int1 / int2 * x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [x0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [x0, x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, x0, x0],
        [x0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, -x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [int0, -int1 / int2 * x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [-x0, -x0, int0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [-int1 / int2 * x0, int0, -int1 / int2 * x0],
        [int0, int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [x0, int0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [int0, -x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, -x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, -x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, -x0, int0],
        [x0, int0, int0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, -x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [int1 / int2 * x0, int0, int0],
        [-int1 / int2 * x0, int0, -int1 / int2 * x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, int0, -x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [-x0, int0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, -x0, -x0],
    ],
    [
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [-int1 / int2 * x0, int0, int0],
        [int0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, x0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [x0, x0, int2 * x0],
        [int0, -x0, int0],
        [-x0, int0, -x0],
    ],
    [
        [int0, x0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, x0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, x0, int2 * x0],
        [-x0, -x0, -x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [int0, int1 / int2 * x0, int0],
        [int0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, x0],
        [x0, x0, int0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [int0, int0, -int1 / int2 * x0],
        [int0, int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, -x0],
        [x0, x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, int0, x0],
        [int0, -x0, x0],
    ],
    [
        [x0, int0, int0],
        [-x0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [int0, x0, -x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [x0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [int0, int0, -int1 / int2 * x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [x0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [x0, x0, int0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, x0, x0],
        [x0, int0, x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, -x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [int0, -int1 / int2 * x0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [x0, int0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [int0, -x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, -x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, -x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, -x0, int0],
        [x0, int0, int0],
    ],
    [
        [x0, x0, int2 * x0],
        [x0, int0, int0],
        [-x0, int0, -x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [x0, x0, int2 * x0],
        [int0, x0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, -x0, int0],
        [int0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, -x0, -x0],
    ],
    [
        [x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int1 / int2 * x0, x0],
        [-int1 / int2 * x0, int0, int0],
        [int0, -int1 / int2 * x0, -int1 / int2 * x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, x0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [x0, x0, int2 * x0],
        [int0, -x0, int0],
        [-x0, int0, -x0],
    ],
    [
        [int0, x0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [x0, x0, int2 * x0],
        [int0, -x0, int0],
        [int0, int0, -x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, -x0, -x0],
    ],
    [
        [x0, x0, int2 * x0],
        [x0, int0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int0],
        [int0, int0, -x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, int0, x0],
        [int0, -x0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, x0],
        [int0, x0, int0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, int0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, -x0],
        [x0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [x0, int0, x0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [int0, int0, -int1 / int2 * x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [x0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [x0, x0, int0],
    ],
    [
        [int1 / int2 * x0, int0, x0],
        [int0, int1 / int2 * x0, -int1 / int2 * x0],
        [int0, int0, -int1 / int2 * x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, int0, x0],
        [int0, -x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, x0, x0],
        [x0, int0, x0],
    ],
    [
        [x0, int2 * x0, int0],
        [-x0, -x0, int0],
        [int0, -x0, x0],
    ],
    [
        [x0, int2 * x0, int0],
        [int0, -x0, int0],
        [int0, -x0, x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
        [int0, -int1 / int2 * x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [-x0, -x0, int0],
    ],
    [
        [-x2, -x2, -int2 * x2],
        [-int2 * x1 + x2, int0, int0],
        [x1, x0, x2],
    ],
    [
        [-x2, -x2, -int2 * x2],
        [int0, -int2 * x1 + x2, int0],
        [x0, x1, x2],
    ],
    [
        [int0, x1 - int2 * x2, int0],
        [-int2 * x0 + x1, int0, int0],
        [x0, x2, x1],
    ],
    [
        [-int2 * x0 + x2, int0, int0],
        [-x2, -x2, -int2 * x2],
        [x0, x1, x2],
    ],
    [
        [x1 - int2 * x2, int0, int0],
        [int0, -int2 * x0 + x1, int0],
        [x2, x0, x1],
    ],
    [
        [int0, x0 - int2 * x1, int0],
        [-x0, -x0, -int2 * x0],
        [x2, x1, x0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [x0, int0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, -x0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, x0],
        [int0, x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [x0, int0, int0],
        [int0, -x0, int0],
    ],
    [
        [int0, -x0, int0],
        [int0, int0, x0],
        [-x0, int0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, -x0, int0],
        [-x0, int0, int0],
    ],
    [
        [int0, -x0, int0],
        [int0, int0, x0],
        [x0, int0, int0],
    ],
    [
        [int0, int0, -x0],
        [x0, int0, int0],
        [int0, x0, int0],
    ],
    [
        [int0, int0, x0],
        [int0, -x0, int0],
        [x0, int0, int0],
    ],
    [
        [x0, x0, int2 * x0],
        [int0, x0, int0],
        [int0, -x0, -x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [x0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, int0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, -x0, -x0],
    ],
    [
        [x0, int0, int0],
        [x0, x0, int2 * x0],
        [-x0, -x0, -x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x0],
    ],
    [
        [x0, x0, int2 * x0],
        [int0, x0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [int0, -x0, int0],
        [-x0, -x0, -int2 * x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int2 * x0],
        [int0, int0, -x0],
    ],
    [
        [x0, int0, int0],
        [x0, x0, int2 * x0],
        [-x0, int0, -x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, x0, int2 * x0],
        [int0, -x0, -x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [x0, int0, int0],
        [int0, x0, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [-x0, int0, int0],
        [x0, x0, x0],
    ],
    [
        [x0, x0, int2 * x0],
        [x0, int0, int0],
        [-x0, int0, -x0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [int0, -x0, -x0],
    ],
    [
        [-x0, int0, int0],
        [int0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [-x0, -x0, -x0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [x0, x0, x0],
    ],
    [
        [int0, x0, int0],
        [-x0, -x0, -int2 * x0],
        [int0, int0, x0],
    ],
    [
        [-x0, int0, int0],
        [int0, -x0, int0],
        [x0, int0, x0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int0],
        [int0, int0, -x0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [int0, -x0, int0],
    ],
    [
        [x0, int0, int2 * x0],
        [-x0, int0, -x0],
        [int0, x0, -x0],
    ],
    [
        [-x0, int0, int0],
        [int0, int0, -x0],
        [x0, x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [int0, int0, x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, -int1 / int2 * x0],
        [-int1 / int2 * x0, -int1 / int2 * x0, int0],
    ],
    [
        [-x0, int0, int0],
        [x0, int0, x0],
        [x0, x0, int0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [-int1 / int2 * x0, -int1 / int2 * x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
    ],
    [
        [-x0, int0, int0],
        [x0, x0, int0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, x0, x0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, x0],
        [int0, x0, int0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, -x0],
        [x0, x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, int0, x0],
        [x0, x0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, int0],
        [int0, x0, -x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, int0, x0],
        [int0, -x0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [int0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [int0, -int1 / int2 * x0, int1 / int2 * x0],
        [int0, -int1 / int2 * x0, int0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x0, x0, x0],
        [x0, int0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [int0, int0, x0],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x0, x0, int0],
        [x0, x0, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [int0, -x0, x0],
        [x0, int0, x0],
    ],
    [
        [-int2 * x0 + x2, int0, int0],
        [x1, int0, -int2 * x0 + int2 * x1 + x2],
        [x0, x2, int0],
    ],
    [
        [int1 / int2 * x1, int0, x1],
        [x2, int1 / int2 * x1 + int2 * x2, -int1 / int2 * x1],
        [x0, int0, -int1 / int2 * x1],
    ],
    [
        [x1 - int2 * x2, int0, int0],
        [x0, int2 * x0 + x1 - int2 * x2, int0],
        [x2, int0, x1],
    ],
    [
        [-x1, -int2 * x1, int0],
        [x2, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, x1, x0],
    ],
    [
        [int1 / int2 * x1, x1, int0],
        [x2, -int1 / int2 * x1, int1 / int2 * x1 + int2 * x2],
        [x0, -int1 / int2 * x1, int0],
    ],
    [
        [x0 - int2 * x1, int0, int2 * x0 - int4 * x1],
        [x2, int0, -x0 + int2 * x1],
        [x1, x0, -x0 + int2 * x1],
    ],
    [
        [x0, int0, int0],
        [int0, x1, int0],
        [int0, int0, x2],
    ],
    [
        [int0, x0, int0],
        [x2, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, x0],
        [int0, x2, int0],
        [x3, int0, x4],
    ],
    [
        [x0 - int2 * x1, int0, int0],
        [x1, x0, int0],
        [int0, int0, x2],
    ],
    [
        [int1 / int2 * x1, x1, int0],
        [x2, -int1 / int2 * x1, int0],
        [int0, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x2, int0, int0],
        [-int1 / int2 * x2, int1 / int2 * x1, x1],
    ],
    [
        [x2, int0, int0],
        [int0, x0 - int2 * x1, int0],
        [int0, x1, x0],
    ],
    [
        [int0, x0 - int2 * x1, int0],
        [x0 - int2 * x2, int0, int0],
        [x2, x1, x0],
    ],
    [
        [x1 - int2 * x2, int0, int0],
        [int0, -int2 * x0 + x1, int0],
        [x2, x0, x1],
    ],
    [
        [x0, int0, int0],
        [x2, x0 + int2 * x2, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [-x0, -int2 * x0, int0],
        [x1, x0, int0],
        [int1 / int2 * x0 + int1 / int2 * x2, x0, x2],
    ],
    [
        [int0, int0, x2],
        [x1, int0, int0],
        [int0, x0, int0],
    ],
    [
        [x0, int0, int0],
        [int0, int0, x2],
        [int0, x1, int0],
    ],
    [
        [int0, x0, int0],
        [int0, int0, x1],
        [x2, int0, int0],
    ],
    [
        [int0, x0, int0],
        [x2, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, int0, int0],
        [int0, x1, int0],
        [int0, int0, x2],
    ],
    [
        [int0, int0, x0],
        [int0, x2, int0],
        [x1, int0, int0],
    ],
    [
        [x4, int0, -int2 * x0],
        [int1 / int2 * x2 - int1 / int2 * x4, x2, x0],
        [x1, int0, x3],
    ],
    [
        [x4, int0, x3],
        [int0, x0, int0],
        [x1, int0, x2],
    ],
    [
        [x1, int2 * x1, int0],
        [x2, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, -int2 * x1, int0],
        [int1 / int2 * x1, x1, x2],
        [x0, int0, int0],
    ],
    [
        [x1, int0, int0],
        [x2, x1 + int2 * x2, int0],
        [int0, int0, x0],
    ],
    [
        [int0, int0, x0],
        [x2, int2 * x2, -int1 / int2 * x0],
        [x1, int0, int0],
    ],
    [
        [int0, int0, -int2 * x0],
        [x2, int0, x0],
        [int1 / int2 * x1, x1, int0],
    ],
    [
        [-int2 * x0, int0, int0],
        [x0, int0, x2],
        [int1 / int2 * x1, x1, int0],
    ],
    [
        [-x2, -x2, -int2 * x2],
        [int0, -int2 * x0 + x2, int0],
        [x1, x0, x2],
    ],
    [
        [int0, x0 - int2 * x1, int0],
        [x0 - int2 * x2, int0, int0],
        [x2, x1, x0],
    ],
    [
        [-x0, -x0, -int2 * x0],
        [x0 - int2 * x2, int0, int0],
        [x2, x1, x0],
    ],
    [
        [x1 - int2 * x2, int0, int0],
        [-x1, -x1, -int2 * x1],
        [x2, x0, x1],
    ],
    [
        [x1 - int2 * x2, int0, int0],
        [int0, -int2 * x0 + x1, int0],
        [x2, x0, x1],
    ],
    [
        [int0, x0 - int2 * x1, int0],
        [-x0, -x0, -int2 * x0],
        [x2, x1, x0],
    ],
    [
        [-x0, int0, -int2 * x0],
        [x1, int0, x0],
        [int1 / int2 * x0 + int1 / int2 * x2, x2, x0],
    ],
    [
        [x0, int2 * x0, int0],
        [x1, -x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x2, -x0, x2],
    ],
    [
        [int1 / int2 * x0, x0, int0],
        [x1, -int1 / int2 * x0, int1 / int2 * x0 + int2 * x1],
        [x2, -int1 / int2 * x0, int0],
    ],
    [
        [x1 - int2 * x2, int0, int0],
        [x0, int2 * x0 + x1 - int2 * x2, int0],
        [x2, int0, x1],
    ],
    [
        [x1, int0, int2 * x1],
        [x2, x1 + int2 * x2, -x1],
        [x0, int0, -x1],
    ],
    [
        [x1, int0, int0],
        [-int1 / int2 * x1 + int1 / int2 * x2, int0, x2],
        [int1 / int2 * x0 - int1 / int2 * x1, x0, int0],
    ],
    [
        [x0, x2, int0],
        [-x2, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x2, x1, int0],
        [x1, -x2, int0],
        [int0, int0, x0],
    ],
    [
        [x1 - int2 * x2, int2 * x0 + int2 * x1 - int4 * x2, int0],
        [x0, -x1 + int2 * x2, int0],
        [x2, -x0 - x1 + int2 * x2, x1],
    ],
    [
        [x0, x0 - x1, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, x1, int0],
        [
            -int1 / int2 * x0 + int1 / int2 * x2,
            -int1 / int2 * x0 + int1 / int2 * x1,
            x2,
        ],
    ],
    [
        [x3, int0, -int2 * x4],
        [int1 / int2 * x1 - int1 / int2 * x3, x1, x4],
        [x0, int0, x2],
    ],
    [
        [x2, x1, int0],
        [-x1, x2, int0],
        [int0, int0, x0],
    ],
    [
        [x2, x1, int0],
        [x1, -x2, int0],
        [int0, int0, x0],
    ],
    [
        [x1, x0, int0],
        [-int1 / int2 * x0, -x0 + x1, int0],
        [
            -int1 / int2 * x1 + int1 / int2 * x2,
            -int1 / int2 * x0,
            x2,
        ],
    ],
    [
        [x2, int2 * x0 + int2 * x2, int0],
        [x0, -x2, int0],
        [int1 / int2 * x1 - int1 / int2 * x2, -x0 - x2, x1],
    ],
    [
        [x0, x2, int0],
        [x2, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, -x2, int0],
        [x2, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1 - int2 * x2, -int2 * x0, int0],
        [x0, int2 * x0 + x1 - int2 * x2, int0],
        [x2, x0, x1],
    ],
    [
        [-x0 + int1 / int2 * x1, x1, int0],
        [x0, x0 - int1 / int2 * x1, int0],
        [
            int1 / int2 * x0
            - int1 / int4 * x1
            + int1 / int2 * x2,
            -int1 / int2 * x1,
            x2,
        ],
    ],
    [
        [-x0, x0, int0],
        [x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, -x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x0, x0, int0],
        [-x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, -x1, int0],
        [x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, x1, int0],
        [x1, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [-x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, int2 * x0, int0],
        [int0, -x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, -x0, x1],
    ],
    [
        [x0, int0, int0],
        [int0, x0, int0],
        [-int1 / int2 * x0 + int1 / int2 * x1, int0, x1],
    ],
    [
        [int0, int2 * x0, int0],
        [x0, int0, int0],
        [int1 / int2 * x1, -x0, x1],
    ],
    [
        [x1, x1, int0],
        [-int1 / int2 * x1, int0, int0],
        [
            int1 / int2 * x0 - int1 / int2 * x1,
            -int1 / int2 * x1,
            x0,
        ],
    ],
    [
        [int0, -int2 * x1, int0],
        [x1, int2 * x1, int0],
        [int1 / int2 * x0, x1, x0],
    ],
    [
        [-x1, -int2 * x1, int0],
        [x1, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, x1, x0],
    ],
    [
        [-int2 * x1, -int2 * x1, int0],
        [x1, int2 * x1, int0],
        [int1 / int2 * x0 + x1, x1, x0],
    ],
    [
        [-x1, int0, int0],
        [x1, x1, int0],
        [int1 / int2 * x0 + int1 / int2 * x1, int0, x0],
    ],
    [
        [int0, x0, int0],
        [x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, int0, int0],
        [int0, -x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, x0, int0],
        [-x0, x0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, -x1, int0],
        [x1, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x1, int0, int0],
        [int0, x1, int0],
        [int0, int0, x0],
    ],
    [
        [x0, -x0, int0],
        [-x0, -x0, int0],
        [int0, int0, x1],
    ],
    [
        [int0, x0, int0],
        [-x0, int0, int0],
        [int0, int0, x1],
    ],
    [
        [x1, x1, int0],
        [x1, -x1, int0],
        [int0, int0, x0],
    ],
]
