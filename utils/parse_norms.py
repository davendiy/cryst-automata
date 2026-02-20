
import sys


def normalize_matrice(m: str):
    m = m.replace(' + ', '+').replace(' - ', '-')
    res = []
    for line in m.splitlines():
        tokens = line.split()
        norm_line = ' '.join(tokens)
        norm_line = norm_line.replace('[ ', '[').replace(' ]', ']').replace(']', '],')
        norm_line = norm_line.replace(' ', ', ')
        res.append('    ' + norm_line)
    return '\n'.join(res)


filename = sys.argv[1]

with open(filename) as file:
    data = file.read()

matrices = data.split('\n\n')
res_m = [normalize_matrice(mat) for mat in matrices]

print('[')
print(*res_m, sep='\n],\n\n[\n')
print(']')
