import cv2
import numpy as np
import rsa
import secrets
import cmath, math
import sys
from mpmath import mp

def update_xyz(x, y, z, r=3.99, beta=6):
    """
    Updates x, y, z according to the quantum logistic map.  Default values of the parameters beta and r are given on page 6 of the paper; the quantum logistic map is given in equations 1-3. Note that x, y, z are complex numbers.

    Uses mp.dps to set precision.
    """
    mp.dps = 100 # set precision
    x_new = r * (x - abs(x)**2) - r * y
    y_new = -y * mp.exp(-2*beta) + mp.exp(-beta) * r * ((2 - x - mp.conjugate(x)) * y - \
            x * mp.conjugate(z) - mp.conjugate(x) * z)
    z_new = -z * mp.exp(-2*beta) + mp.exp(-beta) * r * (2 * (1 - mp.conjugate(x)) * z - \
            2 * x * y - x)
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
    # randomly select three messages (encoded as bytes)
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

def enc(img, pkb):
    """ Calls enc_channel on each channel in the image. """
    image_dims = img.shape
    # if image has one channel (grayscale)
    if len(image_dims)==2:
        return enc_channel(img, pkb)
    # else if image has 3 or 4 channels (RGB or RGBA) 
    if len(image_dims)==3 and (image_dims[2]==3 or image_dims[2]==4): 
        cipher, r, enc_img = ([], [], [])
        for channel in range(image_dims[2]):
            cipher_channel, r_channel, enc_img_channel = enc_channel(img[:,:,channel], pkb)
            cipher.append(cipher_channel)
            r.append(r_channel)
            enc_img.append(enc_img_channel)
        enc_img = np.concatenate(enc_img, axis=2)
        return (cipher, r, enc_img)
    # else invalid image
    return (None, None, None)

def enc_channel(img, pkb):
    """ Takes in one channel of the plain image and the public key, and returns the ciphertexts, r, and the encrypted image. """

    # calculate r
    M, N = img.shape
    r = np.sum((np.concatenate(list(np.arange(j,N+j) for j in range(M))) + img.flatten())**(2/5))

    # calculate m, c
    m, c = gen_ciphertexts(pkb)
    m_int = [int.from_bytes(m_i, sys.byteorder) for m_i in m]
    c_int = [int.from_bytes(c_i, sys.byteorder) for c_i in c]

    # calculate xyz
    x_0, y_0, z_0 = [1/(abs(m_int[i] - c_int[i]) + r) for i in range(3)]
    x, y, z = x_0, y_0, z_0
    for i in range(500):
        x, y, z = np.array(update_xyz(x, y, z, r))

    # Encryption round
    X_rk = img # X_0 = plain image
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

        # row permutations
        Xrk_prime = np.zeros((M,N))
        for i in range(M): # rows
            for j in range(N): # cols
                Xrk_prime[i][j] = X_rk[i][(j+xk_primes[i]-1) % N]
        # column permutations
        Xrk_dp = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                Xrk_dp[i][j] = Xrk_prime[(i+jk_primes[i]-1) % M][j]
        # discrete cosine transform coefficient matrix
        F = np.zeros((M,N))
        for u in range(M):
            for v in range(N):
                for i in range(M):
                    for j in range(N):
                        F[u][v] += Xrk_dp[i][j] * np.cos(((2*i+1)*np.pi*u)/(2*M)) * np.cos(((2*j+1)*np.pi*v)/(2*N))
                F[u][v] *= sigma(u,M) * sigma(v,N)
        # row permutations
        F_prime = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                F_prime[i][j] = F[i][(j+zk_primes[i]-1) % N]
        # column permutations
        F_dp = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                F_dp[i][j] = F_prime[(i+wk_primes[i]-1) % M][j]
        # inverse discrete cosine transform coefficient matrix
        G = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                for u in range(M):
                    for v in range(N):
                        G[i][j] += sigma(u,M) * sigma(v,N) * F_dp[i][j] * np.cos(((2*i+1)*np.pi*u)/(2*M)) * np.cos(((2*j+1)*np.pi*v)/(2*N))
        # generate encrypted image for round rk
        for ind in range(M*N):
            val = G.flatten()[ind] ^ sk_primes[ind] ^ X_rk.flatten()[ind]
            i = ind // N
            j = ind % N
            X_rk[i][j] = val
    # output ciphertexts, r, and encrypted image
    return (c, r, X_rk)

def dec(img, ciphertexts, r, pkb, skb):
    # TODO: implement; should be inverse of enc
    pass

if __name__ == "__main__":
    pk, sk = read_keys("rsa-keys/public.pem", "rsa-keys/private.pem") #test with valid RSA key pair
    #update_xyz(1+3j, 2-4j, 4+2j, r=3.99, beta=6) #random complex number test
    c, r, enc_img = enc(cv2.imread("test_image_1.jpg"), pk)
