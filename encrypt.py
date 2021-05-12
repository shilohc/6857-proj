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
    Updates x, y, z according to the quantum logistic map.  Default values of
    the parameters beta and r are given on page 6 of the paper; the quantum
    logistic map is given in equations 1-3. Note that x, y, z are complex
    numbers.

    Calculates all values mod (N+1)(M+1) (where M and N are the dimensions of
    the input image).
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
    sk = np.array([[np.fix(sum(xyzs[M*j+i+1]) * 1e14) % 256 for j in range(N)] \
            for i in range(M)], dtype=np.int64)

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
    xor_values = np.bitwise_xor(np.bitwise_xor(G.astype('int64'), sk),
            X_rk.astype('int64'))
    return xor_values, xyzs[-1]

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
    """
    # If a is not complex
    if not isinstance(a, complex):
        return a % n

    return complex(a.real % n, a.imag % n)

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
    img_filename = "testimage1_64x48.jpg"
    #img_filename = "testimage1_96x72.jpg"
    #img_filename = "testimage1_32x24.jpg"
    #img_filename = "testimage2_cropped.jpg"
    #img_filename = "testimage2_cropped_big.jpg"
    #img_filename = "testimage2_cropped_less_big.jpg"
    #img_filename = "testimage2_cropped_medium.jpg"
    print("Starting...")
    start = time.time()
    pk, sk = read_keys("rsa-keys/public.pem", "rsa-keys/private.pem")
    img = cv2.imread("images/" + img_filename)
    c, r, enc_img = enc(img, pk, verbose=True)
    print("Time to encrypt:", time.time() - start)
    cv2.imwrite("encrypted/" + img_filename, enc_img)
    print(c, r)

    """
    c = [[b"KUJ\x8d]\xb4\x8cy\x9e\xec_\xc1N\xb8\x8c\xa7\x8ds\xd35\xf1\xf4\xb2o\x8c`\xd6\xfb\x8e8\xd7\x14f\xa5\xccI?\xb1%\xf5\xd7\xc3\xa3\x9b\xeei\xcd\x14c\x1c}y}=O\xdc\x03\x9d\xa7\x88E\x12EE.4Vu#R'ILcb\xa9\xb7o\x1d\x16\xbe\x01S1><}\x18\xf2\xea\xcc\xb3\x84\xe1Z\xe4\xde\x1d\x0c\x1c\x88\xfd\x9b\xd3\xed\xe7\xac`|\xc7\xc1\x07<-\xff\x90\xb5\xd6\xf7\xb7\x1c-\x94\xdd\x189Z\xdb", b'w\xea\xee\x16\x08\xca\xf1\x16\xd6\x8e\x0b\xb7\xce\xb8\x14g\x9d\xcebT\xec\n\xc7t2\x13\xad\xf0\x06AM\x1cz\xc1k\x91=\xe91/\x8dv\r^\xcf`Q\xe2[1\x11`\xe6}\x19\xc1\xa4\xe2\xb7\xf2\x9c\xa3\x1e=p-\xd4\x19@\xc5\x14t\xa0\xf7(3!\xba\xd9\xd6\xde\x9et\xe5\x85\x82\xd2\xd8\xd5\xa6\x02*\x9a\x8cz\x07\x9f\xcdRE\xab\xad\x9f7\xc38C\xdc~\xab\xa6_~\xa6\x9d\x04\x91\xcf\xf1\xe8\xc9\x12H\xa6\xc5\x858\xc1', b'A\x8f\x88R\x00\x99\xae\x96\x0b\xd9\x03\xb4R\x00-\x1e\x13\x8c\xd9\xa8eJ\x865&eM#\x8f\xe3\xdc\xacU\xc7|\x1a\xb9\x9e\xb2(\xe5\xc7\x85\xe0\x9b\xf9\xbdaD\xfd\x1e\x1a\xbdI\xd9Q\xe0\xf7`\xc6"\xff\xcbm\xb3\xc4E\xa0,`\x9d\xa7Z\xec|\x07\xd7\xf64\x1a\xc4#_\xe0j~\xbf2K\xeb,n\x10GSb\x0c\x1a\xe1\xb13\x1e\xb0\xe9j\x00O\x1e,]\x14Bo7\x19$/\xa8z\xad\x19\x91\x98\xcc\xea\xda\xaa\x03'], [b'4,s\x07"?\xae!\xc6RD\xb7\xbf\x17\xde\xba\xa8\x81EC\xa0\x1d\xafHe\x06bp\xc4\x9d\xf7"`0\xf4_\x07\xda\\y\xa4m\x89e\x9eX\x86R|\x0e\xcd\xdd\x1e\xf6~\x7f.\x9b\xdfx \x04\x8a\xe8K\x80\xc4\x10\x1b\xf5\x9a\xa9/ \x82\x0e&-/\x9e\x8ci\xe8\x02\x0b\xf4dE\x92\xf2?!T\x1c\xee\t\x13}\x89\xecv\xa4\xf1\xd3)\x96l\xbd\xe2|\x0c\x82\xf8\xe5t\x08T\x17?\x87\xf4\xd9AX\x00 z\xbd', b"\x95+\xd0@\x8f\x04h'\x8e\x90b\x90\xa4F\x0cmm\xa7W\x12\xac\x18\xf3\x9a\xedS\x94\xb6\xac\xb3\x1d\xc6\xc7\xf7z9\xba\x89v\xc3\x86$ko\x11'd\xd3\xd8\xc1\x95Sq\xcd\x0fH\xdc\xdb\xea\xfa\xd5\xd4%$R8ePj\xc1\xc6\xdc\xd39\xe0t\xcbGF\x96U6h\xa3}\x02\x11\x98]\xb9*y\xf9\x0b\xec\xeb\x82\x94\x8cj\x8a\xde\xe8\xfb\xeeN=\x9eY\x16I\x8f\xacJS_\xc0+\xf4\x162\xea\x85<a;:\x0c", b'\x1c\x9cK\x95u\x93\x94\t^&M3S\xc3\xe7\xde\x10\xf7\x0c\xfcC\x8e\x89\xff\xc8\xbd\x88[\xf8\x12a\xd3\x00MS\x1f\xbbI\x7f\x15xd\xd6\xbc\n\xdbN:4\xb1\t\xf2\xa0\xec\xc1[\xef$\x10h\x98\xf0\x9d\x195\xf9\x1c\xca^\x93\x17\xd8\x95\xb6\xf0\xcdy\xc2\x10\xdc\xd8\x07\x03R\x9aH@\x15\xc4#*|\n}y\xa3\xcd\xf2\xef\xeb\xf4\xc6\x15c\xa1\xd9\xe1\xb0\x8b\x01>\x8e\xb5>\xbfm[\xcf\x94G\x7f\xeb\xe1Rj|gA'], [b'\t\x8f\xd3\xe2+*\x80\xaa\xd9\xce\x8a\xd5\x99]\x95;>H\xec\x9a\xcf"\x92`\xf9\xa8\xc2\xfe\xfe\xb5S\xdf\x17\x17\t\x03i\x89\xa6\xaeJ\x0cM\xba\xb1\xf9\x85Xk|Bq\xf5`\x9c{\xb3\xc3\xc7\n\x86\xc7b\xeb\x93\x16$0\xbe\xed\xe1%\xf4>\xb5|\xae\xa69\xef\xe2\xe4\x04\x99\xf2\xa2\xd3\x92\xeb\xfe}\x1ei\x12\xc0I\xe5o\x15\xdb\xa8Ci\x08r\x8a.*\xdd\x866\xb4\xdf\x01\xa1@S!\xdb\xea\xfaC\xe3\r\x1f\x15\x85S', b'e\xdc\xba\xecZ\xd8\xf60\xf8\xcc\xdd\xc69\xb4\xbe\x86\xb0\\\x12\x11\xd2\xb7\xd6\x92c\xc2I\x06\xb3\r\x9e/q\xa0\x83\xa8\xddc\x96pI\xcd\xad\x03\xf4\xf3\xd0\x95(#\xed\x7fY4B<\x98~:\xa3\xf1\xf7\xb1[\xaa\xb9\xb6\xed=\rt\xe6\xe9@\x05(G\xfc\xbd\xf5o\xe0\x01Lb\xe0\xb1\xc4\x8f\x9b\xde\x0c\x01\xc9\x9b\x1c7$Qp\xcb\x80\xff4\xd2\xbc)\xb6\xe2n\xb0X\x92\ty\x90P\xceEC!6\x99\x90TqZ\xc5', b'lG\x8cj\x9a|\x802\xf9De\x8fw\xdd\x14d\x9d\x18\x93x\x85\x12\x80\xf1\xd8\xfc\x18\xf8d\xc0&`m\x04w\x9fl\xf81JS\xa7dC\xb9\xe1\x16j\xfa\xd5\x1f%\xac\xa4\xe9X\xee\xccu\xcb*\x89\xf6\xcc\x99[\xe4\xb4\xf3q\xa5K\xec\xbb\xf0\xb28Nb\xff_\xac\xf1m\x80\x17i\x006v\xb7\x86\x19n]\x11\x03\xdc\xecL\xb9\xe9\xc0\xbf\x819\xe6\x041Z\xa3\xe9\xa9D\x01\xa9\x87\x0f\r\x86\xd0x\xfb \xeb[\xf7F']]
    r = [21670.49893146257, 22007.01037841568, 21804.01480001869]

    enc_filename = "testimage1_64x48_blur_only.jpg"
    enc_img = cv2.imread("encrypted/" + enc_filename)
    start = time.time()
    dec_img = dec(enc_img, c, r, pk, sk, verbose=False)
    print("Time to decrypt:", time.time() - start)
    cv2.imwrite("results/" + enc_filename, dec_img)
    """
