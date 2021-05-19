# To run this code, download Paillier Homomorphic Encryption algorithm from
# https://github.com/chronarchitect/Homomorphic-Image-Encryption
import Paillier as P
import ImageCryptography as IC
from PIL import Image
import numpy as np

def noisy_encryption(fn):
    # get image and keys
    img = Image.open(fn)
    pk, sk = P.generate_keys()
    # generate cipher and add noise
    cipher = IC.ImgEncrypt(pk, img)
    noise = np.random.choice([-1, 0, 1], p=[0.1, 0.8, 0.1], size=cipher.shape)
    noisy_cipher = cipher + noise
    # decrypt cipher and noisy cipher and display
    decoded_noisy = IC.ImgDecrypt(pk, sk, noisy_cipher)
    decoded = IC.ImgDecrypt(pk, sk, cipher)
    decoded_noisy.show()
    decoded.show()
    # return decrypted images
    return (decoded, decoded_noisy)

