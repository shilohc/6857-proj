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

def dct(img):
    """
    Computes the discrete cosine transform of the array img.  Assumes that img
    is grayscale (single-channel).
    """
    M, N = img.shape
    out = np.zeros((M,N))
    i_array = np.array([[i for j in range(N)] for i in range(M)])
    j_array = np.array([[j for j in range(N)] for i in range(M)])
    sigma_array = np.array([[sigma(i, M) * sigma(j, N) for j in range(N)] \
            for i in range(M)])
    for u in range(M):
        for v in range(N):
            cos_array = np.cos(u * (2*i_array + 1) * np.pi)/(2*M) * \
                        np.cos(v * (2*j_array + 1) * np.pi)/(2*N)
            out[u][v] = np.sum(img * cos_array)
    return out * sigma_array

def inverse_dct(img):
    """
    Computes the inverse discrete cosine transform of the array img.  Assumes
    that img is grayscale (single-channel).
    """
    M, N = img.shape
    out = np.zeros((M,N), dtype=np.float32)
    u_array = np.array([[u for v in range(N)] for u in range(M)])
    v_array = np.array([[v for v in range(N)] for u in range(M)])
    sigma_array = np.array([[sigma(u, M) * sigma(v, N) for v in range(N)] \
            for u in range(M)])
    for i in range(M):
        for j in range(N):
            cos_array = np.cos(u_array * (2*i + 1) * np.pi)/(2*M) * \
                        np.cos(v_array * (2*j + 1) * np.pi)/(2*N)
            out[i][j] = np.sum(sigma_array * cos_array)
    return img * out

def permute_rows(img, primes):
    M, N = img.shape
    out = np.zeros((M, N))
    for i in range(M):
        for j in range(N):
            out[i][j] = img[i][int((j + primes[i] - 1) % N)]
    return out

def permute_cols(img, primes):
    M, N = img.shape
    out = np.zeros((M, N))
    for i in range(M):
        for j in range(N):
            out[i][j] = img[int((i + primes[j] - 1) % M)][j]
    return out

def encryption_round(X_rk, xyz_prev, mod_n, r, verbose=False):
    M, N = X_rk.shape

    # xyzs = [xyz_{500+rk(MN)}, xyz_{500+rk(MN)+1}, ..., xyz_{500+rk(MN)+MN}
    xyzs = [list(xyz_prev)] 
    for i in range(M*N):
        xyzs.append(update_xyz(xyzs[-1], mod_n, r))

    # Calculate x_k', y_k', z_k', w_k', and s_k
    xk_primes = [(np.floor(xyzs[k+1][0] * 1e14)) % (N+1) for k in range(M)]
    yk_primes = [(np.floor(xyzs[k+1][1] * 1e14)) % (M+1) for k in range(N)]
    zk_primes = [(np.floor((xyzs[k+1][2] * 0.6 + xyzs[k+1][0] * 0.4) * 1e14)) \
            % (N+1) for k in range(M)]
    wk_primes = [(np.floor((xyzs[k+1][2] * 0.6 + xyzs[k+1][1] * 0.4) * 1e14)) \
            % (M+1) for k in range(N)]
    sk = np.array([np.fix(sum(xyzs[k+1]) * 1e14) % 256 for k in range(M*N)], \
            dtype=np.float32)

    if verbose:
        print("Starting row and column permutations...")
    Xrk_permuted = permute_cols(permute_rows(X_rk, xk_primes), yk_primes)
    if verbose:
        print("Starting DCT coefficient matrix...")
    F = dct(Xrk_permuted)
    if verbose:
        print("Starting row and column permutations...")
    F_permuted = permute_cols(permute_rows(F, zk_primes), wk_primes)
    if verbose:
        print("Starting inverse DCT coefficient matrix...")
    G = inverse_dct(F_permuted)
    if verbose:
        print("Generating encrypted image...")
    xor_values = xor(xor(G.flatten(), sk), X_rk.flatten())
    return np.reshape(xor_values, (M, N)), xyzs[-1]

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
    """
    Helper function called by dct and inverse_dct.
    """
    if x==0:
        return np.sqrt(1/L)
    return np.sqrt(2/L)

def mod(a, n):
    """
    Returns `a mod n`, where `a` can be a complex or real number. 

    Implementation for modular arithmetic over complex numbers taken from this
    stack exchange thread: https://codegolf.stackexchange.com/questions/196122/gaussian-integer-division-reminder
    """
    # TODO: check this again too
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
    Calculates the absolute value of the XOR of two arrays `x` and `y` by converting both to ints and returning their magnitude.
    """
    # TODO: look at this again, may be causing the encryption bug
    # why does smaller image look encrypted?  not just smaller, the
    # color distribution looks better
    # try to make a 64x48 test image, various other sizes of test image
    return np.abs(np.array(x, dtype=int) ^ np.array(y, dtype=int))

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

    # X_rk holds the value of the image at start of round rk; X_0 = plain image
    X_rk = np.array(img, dtype=np.float32) 
    # xyz_prev holds the value xyz_{500 + rk(MN)} (xyz_500 before rk loop)
    xyz_prev = xyz 

    for rk in range(5):
        if verbose:
            print("In encryption round {}...".format(rk))
        X_rk, xyz_prev = encryption_round(X_rk, xyz_prev, mod_n, r, verbose)

    # Output ciphertexts, r, and encrypted image
    return (c, r, X_rk)

def dec(img, ciphertexts, r, pkb, skb, verbose=False):
    """ Calls dec_channel on each channel in the img. """
    if verbose:
        print("Starting encryption...")

    image_dims = img.shape
    # If image has one channel (grayscale)
    if len(image_dims)==2:
        return dec_channel(img, ciphertexts, r, pkb, skb)

    # Else if image has 3 or 4 channels (RGB or RGBA) 
    if len(image_dims)==3 and (image_dims[2]==3 or image_dims[2]==4): 
        dec_img = []
        for channel in range(image_dims[2]):
            dec_img_channel = dec_channel(img[:,:,channel], ciphertexts[channel], r[channel], pkb, skb, verbose=verbose)
            dec_img.append(dec_img_channel.T) # (rows, cols) => (cols, rows)
        dec_img = np.array(dec_img).T # (channels, cols, rows) => (rows, cols, channels)
        return dec_img

    # Else invalid image
    return (None, None, None)

def dec_channel(img, ciphertexts, r, pkb, skb, verbose=False):
    if verbose:
        print("Encrypting new channel...")
    # Calculate r
    M, N = img.shape
    if verbose:
        print("Calculated r...")
    mod_n = (M+1)*(N+1) if (M+1)*(N+1)>256 else (M+1)*(N+1)*256

    # Calculate m, c
    c = ciphertexts
    m = [rsa.decrypt(i, skb) for i in c]
    m_int = [int.from_bytes(m_i, sys.byteorder) for m_i in m]
    c_int = [int.from_bytes(c_i, sys.byteorder) for c_i in c]
    if verbose:
        print("Calculated m, c...")

    # Calculate xyz_500
    xyz = [mod(1/(abs(m_int[i] - c_int[i]) + r), mod_n) for i in range(3)]
    for i in range(500):
        xyz = update_xyz(xyz, mod_n, r)

    # C_rk is the encrypted image at start of round rk
    # C_0 = original encrypted image
    C_rk = np.array(img, dtype=np.float32) 
    # xyz_prev holds the value xyz_{500 + rk(MN)} (xyz_500 before rk loop)
    xyz_prev = xyz 

    for rk in range(5):
        if verbose:
            print("In decryption round {}...".format(rk))
        C_rk, xyz_prev = encryption_round(C_rk, xyz_prev, mod_n, r, verbose)

    # Output decrypted image
    return C_rk

if __name__ == "__main__":
    #img_filename = "testimage1_128x96.jpg"
    #img_filename = "testimage1_64x48.jpg"
    img_filename = "testimage1_96x72.jpg"
    #img_filename = "testimage1_32x24.jpg"
    print("Starting...")
    start = time.time()
    pk, sk = read_keys("rsa-keys/public.pem", "rsa-keys/private.pem")
    c, r, enc_img = enc(cv2.imread("images/" + img_filename), pk, verbose=True)
    print("Time to encrypt:", time.time() - start)
    cv2.imwrite("encrypted/" + img_filename, enc_img)

    start = time.time()
    dec_img = dec(enc_img, c, r, pk, sk, verbose=False)
    print("Time to decrypt:", time.time() - start)
    cv2.imwrite("results/" + img_filename, dec_img)
