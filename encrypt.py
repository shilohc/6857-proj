import cv2
import numpy as np
import rsa
import secrets

def update_xyz(x, y, z, r):
    # TODO: implement
    # x, y, z are complex numbers!
    # use equations 1, 2, 3
    pass

def gen_ciphertexts(pkb):
    # randomly select three messages (encoded as bytes)
    messages = [secrets.token_bytes(53) for i in range(3)]
    return messages, [rsa.encrypt(m, pkb) for m in messages]

def enc(img, pkb):
    rows, cols = img.shape
    M, N = img.shape
    r = 0
    for i in range(rows):
        for j in range(cols):
            # TODO: what's going on with image channels??
            # using just the blue channel for now
            pixel_b, pixel_g, pixel_r = img[i][j]
            r += ((pixel_b + i + j)**2)**(1/5)
    m, c = gen_ciphertexts(pkb)
    x_0, y_0, z_0 = [1/(abs(m[i] - c[i]) + r) for i in range(3)]
    x, y, z = x_0, y_0, z_0
    for i in range(500):
        x, y, z = update_xyz(x, y, z, r)

    # Encryption round
    for rk in range(5):
        xyzs = [[x, y, z]]
        for i in range(M*N):
            xyzs.append(update_xyz(*xyzs[-1], r))
        xk_primes = []
        for k in range(1, M+1):
            xk_primes.append((np.floor(xyzs[k+1][0] * 1e14)) % (N+1))
        yk_primes = []
        for k in range(1, N+1):
            yk_primes.append((np.floor(xyzs[k+1][1] * 1e14)) % (M+1))
        zk_primes = []
        for k in range(1, M+1):
            zk_primes.append((np.floor((xyzs[k+1][2] * 0.6 + xyzs[k+1][0] * 0.4) * 1e14)) % (N+1))
        wk_primes = []
        for k in range(1, N+1):
            wk_primes.append((np.floor((xyzs[k+1][2] * 0.6 + xyzs[k+1][1] * 0.4) * 1e14)) % (M+1))
        sk_primes = []
        for k in range(1, (M*N)+1):
            sk_primes.append((np.fix(sum(xyzs[k+1])) * 1e14) % 256)

        # TODO: row permutations
        # TODO: column permutations
        # TODO: discrete cosine transform coefficient matrix
        # TODO: row permutations
        # TODO: column permutations
        # TODO: inverse discrete cosine transform coefficient matrix
        # TODO: generate encrypted image for round rk
    # TODO: output ciphertexts, r, and encrypted image

def dec(img, ciphertexts, r, pkb, skb):
    # TODO: implement; should be inverse of enc
    pass

