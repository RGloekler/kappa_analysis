import numpy as np
import sys

def eqn(r):
    return (-0.0025 * r ** 2) + 0.0743*r + 0.9936

def main():
    radius = 4.25
    k_ij =  0.4194620036963194

    F = eqn(radius)

    k = k_ij / F
    print(k)
    return

if __name__ == "__main__":
    main()
