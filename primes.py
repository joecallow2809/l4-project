def tri_numbers(n):
    tri_number = 0
    for i in range(n+1):
        tri_number += i
    return tri_number

print tri_numbers(3)
