import cv2
import numpy as np
import rsa
import secrets
import cmath, math

def update_xyz(x, y, z, r=3.99, beta=6):
    """
    Updates x, y, z according to the quantum logistic map.  Default values of
    the parameters beta and r are given on page 6 of the paper; the quantum
    logistic map is given in equations 1-3.
    Note that x, y, z are complex numbers.
    """

    #X
    x_new = r*(x - abs(x)**2) - r*y

    #Y
    y_new = -y*math.exp(-2*beta) + math.exp(-beta)*r*((2-x-x.conjugate())*y - \
            x*z.conjugate() - x.conjugate()*z)
    
    #Z
    z_new = -z*math.exp(-2*beta) + math.exp(-beta)*r*(2*(1-x.conjugate())*z - \
            2*x*y - x)
    
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

def read_keys(public_filename, private_filename):
    """  takes any two public and private key files 
    returns rsa.key.PublicKey and rsa.key.PrivateKey types
    easily converted to string if need be"""
    #TESTED and working!

    with open(public_filename, mode='rb') as public_file:
        key_data = public_file.read()
        public_key = rsa.PublicKey.load_pkcs1_openssl_pem(key_data)

    with open(private_filename, mode='rb') as private_file:
        key_data = private_file.read()
        private_key = rsa.PrivateKey.load_pkcs1(key_data)

    return public_key, private_key

def sigma(x, L):
    """ Helper function called by enc(img, pkb)."""
    if x==0:
        return np.sqrt(1/L)
    return np.sqrt(2/L)

def enc(img, pkb):
    # Assumes a grayscale (one-channel) image.  I think we're supposed to
    # encrypt each channel separately.
    M, N = img.shape
    r = 0
    for i in range(M): # rows
        for j in range(N): # cols
            r += ((img[i][j] + i + j)**2)**(1/5)
    m, c = gen_ciphertexts(pkb)
    x_0, y_0, z_0 = [1/(abs(m[i] - c[i]) + r) for i in range(3)]
    x, y, z = x_0, y_0, z_0
    for i in range(500):
        x, y, z = update_xyz(x, y, z, r)

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

        # TODO(ray)
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
        # TODO: generate encrypted image for round rk
        #X_rk = np.zeros((M,N))
        #for i in range(M):
    # TODO: output ciphertexts, r, and encrypted image

def dec(img, ciphertexts, r, pkb, skb):
    # TODO: implement; should be inverse of enc
    pass


if __name__ == "__main__":
    read_keys("rsa-keys/public.pem", "rsa-keys/private.pem") #test with valid RSA key pair
    update_xyz(1+3j, 2-4j, 4+2j, r=3.99, beta=6) #random complex number test

