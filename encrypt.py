import cv2
import numpy as np
import rsa
import secrets
import sys
import struct
import time
#import cmath, math
#from mpmath import mp

def update_xyz(xyz, mod_n, r=3.99, beta=6):
    """
    Updates x, y, z according to the quantum logistic map.  Default values of the parameters beta and r are given on page 6 of the paper; the quantum logistic map is given in equations 1-3. Note that x, y, z are complex numbers.

    Calculates all values mod (N+1)(M+1) (where M and N are the dimensions of the input image).
    """
    x, y, z = xyz
    x_new = mod(r * (x - abs(x)**2) - r * y, mod_n)
    y_new = mod(-y * np.exp(-2*beta) + np.exp(-beta) * r * ((2 - x - np.conjugate(x)) * y - x * np.conjugate(z) - np.conjugate(x) * z), mod_n)
    z_new = mod(-z * np.exp(-2*beta) + np.exp(-beta) * r * (2 * (1 - np.conjugate(x)) * z - 2 * x * y - x), mod_n)
    return x_new, y_new, z_new

def dct(u, v, img):
    """
    Computes the discrete cosine transform at (u, v).  Assumes a grayscale
    (one-channel) image.
    """
    M, N = img.shape
    sigma_u = np.sqrt(2/M)
    sigma_v = np.sqrt(2/N)
    if u == 0:
        sigma_u = np.sqrt(1/M)
    if v == 0:
        sigma_v = np.sqrt(1/N)

    # TODO(shiloh): should be elementwise
    cos_array = np.array([[(np.cos(u * (2*i + 1) * np.pi / (2*M)) \
            * np.cos(v * (2*j + 1) * np.pi / (2*N))) \
            for j in range(N)] for i in range(M)])
    dct_mat = sigma_u * sigma_v * img * cos_array
    dct = np.sum(dct_mat)
    return sigma_u * sigma_v * dct

def inverse_dct(u, v, dct):
    """
    Computes the inverse discrete cosine transform at (u, v).  Assumes a
    grayscale (one-channel) image.
    """
    M, N = img.shape
    sigma_u = np.sqrt(2/M)
    sigma_v = np.sqrt(2/N)
    if u == 0:
        sigma_u = np.sqrt(1/M)
    if v == 0:
        sigma_v = np.sqrt(1/N)

    # TODO(shiloh): should be elementwise
    cos_array = np.array([[(np.cos(u * (2*i + 1) * np.pi / (2*M)) \
            * np.cos(v * (2*j + 1) * np.pi / (2*N))) \
            for j in range(N)] for i in range(M)])
    return sigma_u * sigma_v * dct * cos_array

def gen_ciphertexts(pkb):
    """
    Randomly select three messages, encoded as bytes, then use the messages and the public key to generate three ciphertexts. The messages can be at most (length of public key in bytes - 11) bytes long. 

    Returns the messages and ciphertexts. 
    """
    messages = [secrets.token_bytes(53) for i in range(3)]
    return messages, [rsa.encrypt(m, pkb) for m in messages]

def read_keys(public_filename, secret_filename):
    """  takes any two public and secret key files 
    returns rsa.key.PublicKey and rsa.key.PrivateKey types
    easily converted to string if need be"""
    #TESTED and working!

    with open(public_filename, mode='rb') as public_file:
        key_data = public_file.read()
        public_key = rsa.PublicKey.load_pkcs1_openssl_pem(key_data)

    with open(secret_filename, mode='rb') as secret_file:
        key_data = secret_file.read()
        secret_key = rsa.PrivateKey.load_pkcs1(key_data)

    return public_key, secret_key

def sigma(x, L):
    """ Helper function called by enc_channel."""
    if x==0:
        return np.sqrt(1/L)
    return np.sqrt(2/L)

def enc(img, pkb, verbose=False):
    """ Calls enc_channel on each channel in the image. """
    if verbose:
        print("Starting encryption...")

    image_dims = img.shape
    # If image has one channel (grayscale)
    if len(image_dims)==2:
        return enc_channel(img, pkb)

    # Else if image has 3 or 4 channels (RGB or RGBA) 
    if len(image_dims)==3 and (image_dims[2]==3 or image_dims[2]==4): 
        cipher, r, enc_img = ([], [], [])
        for channel in range(image_dims[2]):
            cipher_channel, r_channel, enc_img_channel = enc_channel(img[:,:,channel], pkb, verbose=verbose)
            cipher.append(cipher_channel)
            r.append(r_channel)
            enc_img.append(enc_img_channel.T) # (rows, cols) => (cols, rows)
        enc_img = np.array(enc_img).T # (channels, cols, rows) => (rows, cols, channels)
        return (cipher, r, enc_img) 

    # Else invalid image
    return (None, None, None)

def mod(a, n):
    """
    Returns `a mod n`, where `a` can be a complex or real number. 

    Implementation for modular arithmetic over complex numbers taken from this stack exchange thread: https://codegolf.stackexchange.com/questions/196122/gaussian-integer-division-reminder
    """
    # If a is not complex
    if not isinstance(a, complex):
        return a % n

    # If a is complex
    # Compute p = a * conj(b)
    a1, a2 = (a.real, a.imag)
    n1, n2 = (n.real, n.imag)
    p1 = a1 * n1 + a2 * n2
    p2 = - a1 * n2 + a2 * n1

    # Compute the norm-squared of n
    n_nsq = n1 ** 2 + n2 ** 2
    n_nsq_half = n_nsq // 2

    # Do a symmetric mod of d1 and d2 by n_nsq
    # That is, into the interval [-n_nsq_half, +n_nsq_half)
    d1 = (p1 + n_nsq_half) % n_nsq - n_nsq_half
    d2 = (p2 + n_nsq_half) % n_nsq - n_nsq_half

    # Compute r = d / n
    r1 = (d1 * n1 - d2 * n2) // n_nsq
    r2 = (d1 * n2 + d2 * n1) // n_nsq
    return complex(r1, r2)

def xor(x, y):
    """
    Calculates the XOR of numpy arrays `x` and `y`. They must contain values of type np.float32, and be the same length.

    Implementation taken from this stack exchange thread: https://codegolf.stackexchange.com/questions/192862/floating-point-xor
    """
    xor_vals = (x.view("i")^y.view("i")).view("f")

    # If xor with floats returns NAN
    # return xor of int values
    nan_inds = np.argwhere(np.isnan(xor_vals))
    for i in nan_inds:
        xor_vals[i] = int(x[i]) ^ int(y[i])

    return xor_vals

def enc_channel(img, pkb, verbose=False):
    """ Takes in one channel of the plain image and the public key, and returns the ciphertexts, r, and the encrypted image. """

    if verbose:
        print("Encrypting new channel...")
    # Calculate r
    M, N = img.shape
    r = np.sum((np.concatenate(list(np.arange(j,N+j) for j in range(M))) + img.flatten())**(2/5))
    if verbose:
        print("Calculated r...")
    mod_n = (M+1)*(N+1) if (M+1)*(N+1)>256 else (M+1)*(N+1)*256

    # Calculate m, c
    m, c = gen_ciphertexts(pkb)
    m_int = [int.from_bytes(m_i, sys.byteorder) for m_i in m]
    c_int = [int.from_bytes(c_i, sys.byteorder) for c_i in c]
    if verbose:
        print("Calculated m, c...")

    # Calculate xyz_500
    xyz = [mod(1/(abs(m_int[i] - c_int[i]) + r), mod_n) for i in range(3)]
    for i in range(500):
        xyz = update_xyz(xyz, mod_n, r)

    # Encryption round variables
    # X_rk holds the value of the image at start of round rk; X_0 = plain image
    X_rk = np.array(img, dtype=np.float32) 
    # xyz_prev holds the value xyz_{500 + rk(MN)} (xyz_500 before rk loop)
    xyz_prev = xyz 

    # Encryption round
    for rk in range(5):
        if verbose:
            print("In encryption round {}...".format(rk))
        # xyzs = [xyz_{500+rk(MN)}, xyz_{500+rk(MN)+1}, ..., xyz_{500+rk(MN)+MN}
        xyzs = [list(xyz_prev)] 
        for i in range(M*N):
            xyzs.append(update_xyz(xyzs[-1], mod_n, r))
        xyz_prev = xyzs[-1] # update xyz_prev to last xyzs array

        # Calculate x_k', y_k', z_k', w_k', and s_k
        xk_primes = [] # xk_primes = [x_1', x_2', ..., x_M']
        for k in range(M):
            xk_primes.append((np.floor(xyzs[k+1][0] * 1e14)) % (N+1))
        yk_primes = [] # yk_primes = [y_1', y_2', ..., y_N']
        for k in range(N):
            yk_primes.append((np.floor(xyzs[k+1][1] * 1e14)) % (M+1))
        zk_primes = [] # zk_primes = [z_1', z_2', ..., z_M']
        for k in range(M):
            zk_primes.append((np.floor((xyzs[k+1][2] * 0.6 + xyzs[k+1][0] * 0.4) * 1e14)) % (N+1))
        wk_primes = [] # wk_primes = [w_1', w_2', ..., w_M']
        for k in range(N):
            wk_primes.append((np.floor((xyzs[k+1][2] * 0.6 + xyzs[k+1][1] * 0.4) * 1e14)) % (M+1))
        sk = [] # sk = [s_1, s_2, ..., s_{MN}]
        for k in range(M*N):
            sk.append(np.fix(sum(xyzs[k+1]) * 1e14) % 256)

        if verbose:
            print("Starting row and column permutations...")
        # Row permutations
        Xrk_prime = np.zeros((M,N))
        for i in range(M): # rows
            for j in range(N): # cols
                Xrk_prime[i][j] = X_rk[i][int((j+xk_primes[i]-1) % N)]

        # Column permutations
        Xrk_dp = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                Xrk_dp[i][j] = Xrk_prime[int((i+yk_primes[i]-1) % M)][j]

        # Discrete cosine transform coefficient matrix
        F = np.zeros((M,N))
        for u in range(M):
            for v in range(N):
                for i in range(M):
                    for j in range(N):
                        F[u][v] += Xrk_dp[i][j] * np.cos(((2*i+1)*np.pi*u)/(2*M)) * np.cos(((2*j+1)*np.pi*v)/(2*N))
                F[u][v] *= sigma(u,M) * sigma(v,N)

        # Row permutations
        F_prime = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                F_prime[i][j] = F[i][int((j+zk_primes[i]-1) % N)]

        # Column permutations
        F_dp = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                F_dp[i][j] = F_prime[int((i+wk_primes[i]-1) % M)][j]

        if verbose:
            print("Starting inverse discrete cosine transform coefficient matrix...")
        # Inverse discrete cosine transform coefficient matrix
        G = np.zeros((M,N), dtype=np.float32)
        for i in range(M):
            for j in range(N):
                for u in range(M):
                    for v in range(N):
                        G[i][j] += sigma(u,M) * sigma(v,N) * F_dp[i][j] * np.cos(((2*i+1)*np.pi*u)/(2*M)) * np.cos(((2*j+1)*np.pi*v)/(2*N))

        if verbose:
            print("Generating encrypted image for round {}...".format(rk))
        # Generate encrypted image for round rk
        xor_values = xor(xor(G.flatten(), np.array(sk, dtype=np.float32)), X_rk.flatten())
        for ind in range(M*N):
            val = xor_values[ind]
            i = int(ind // N)
            j = int(ind % N)
            X_rk[i][j] = val

    # Output ciphertexts, r, and encrypted image
    return (c, r, X_rk)

def dec(img, ciphertexts, r, pkb, skb):
    # TODO: implement; should be inverse of enc
    pass

if __name__ == "__main__":
    print("Starting...")
    start = time.time()
    pk, sk = read_keys("rsa-keys/public.pem", "rsa-keys/private.pem")
    print("Read keys...")
    c, r, enc_img = enc(cv2.imread("images/testimage1_32x24.jpg"), pk, verbose=True)
    print(enc_img)
    end = time.time()
    print("Time:", end - start)
