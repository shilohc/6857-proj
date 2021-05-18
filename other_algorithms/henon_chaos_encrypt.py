from PIL import Image, ImageOps
import os
import numpy as np

def dec(bitSequence):
    decimal = 0
    for bit in bitSequence:
        decimal = decimal * 2 + int(bit)
    return decimal

#Takes size of the image as an input
def genTransformationMatrix(m):
    #Replace Hardcoded Pixel Values
    #Serves as the initial Parameter and also the symmetric secret key
    x = 0.1
    y = 0.1
    sequenceSize = m * m * 8 #Total Number of bitSequence produced
    bitSequence = []    #Each bitSequence contains 8 bits
    byteArray = []      #Each byteArray contains m( i.e 512 in this case) bitSequence
    TImageMatrix = []   #Each TImageMatrix contains m*n byteArray( i.e 512 byteArray in this case)
    for i in range(sequenceSize):
        #Henon Map formula
        xN = y + 1 - 1.4 * x**2
        yN = 0.3 * x

        # New x = xN and y = yN
        x = xN
        y = yN

        # Each Value of xN is converted into 0 or 1 based on the threshold value
        if xN <= 0.3992:
            bit = 0
        else:
            bit = 1
        #bits are inserted into bitSequence
        bitSequence.append(bit)
        #Each bitSequence is converted into a decimal number
        #This decimal number is inserted into byteArray
        if i % 8 == 7:
            decimal = dec(bitSequence)
            byteArray.append(decimal)
            bitSequence = []
        #ByteArray is inserted into TImageMatrix
        byteArraySize = m*8
        if i % byteArraySize == byteArraySize-1:
            TImageMatrix.append(byteArray)
            byteArray = []
            
    return TImageMatrix

def decryptHenonImage(imageMatrix):
    transformationMatrix = genTransformationMatrix(len(imageMatrix))

    henonDecryptedImage = []
    for i in range(len(imageMatrix)):
        row = []
        for j in range(len(imageMatrix)):
            row.append(imageMatrix[i][j] ^ transformationMatrix[i][j])
        henonDecryptedImage.append(row)

    return henonDecryptedImage

def pixelManipulation(size, imageMatrix):
    transformationMatrix = genTransformationMatrix(size)

    #Performing Ex-Or Operation between the transformation Matrix and ImageMatrix
    #Storing the result in resultantMatrix
    resultantMatrix = []
    for i in range(size):
        row = []
        for j in range(size):
            row.append(transformationMatrix[i][j] ^ imageMatrix[i][j])
        resultantMatrix.append(row)

    return resultantMatrix

def combine_channels(image_matrix):
    """
    Take in 3D array of shape (3, width, height) and convert to (width, height, 3). Can't use numpy b/c it breaks the encryption for some reason.
    """
    channels, width, height = np.array(image_matrix).shape
    im_t = []
    for i in range(width):
        row = []
        for j in range(height):
            row.append([image_matrix[c][i][j] for c in range(channels)])
        im_t.append(row)
    return im_t

def convert_to_array(im, channel=None):
    """
    Take in PIL Image and convert it to an array. If channel is not none, index into that channel.
    """

    pix = im.load()
    image_shape = im.size # Get the width and height of the image for iterating 

    image_matrix = []
    for width in range(int(image_shape[0])):
        row = []
        for height in range(int(image_shape[1])):
            if channel is None:
                row.append((pix[width,height]))
            else:
                row.append((pix[width,height][channel]))
        image_matrix.append(row)

    return image_matrix

def get_image_matrix(im):
    """
    Take in PIL Image and convert it to an array. 

    If image is a 3D array of shape (width, height, 3), convert it to an image matrix of shape (3, width, height). 

    Can't use numpy because it breaks the encryption for some reason.
    """
    pix = im.load()
    img_shape = im.size

    if len(img_shape)==3:
        for channel in range(channels):
            channel_matrix = convert_to_array(im, channel)
            encr_channel = pixelManipulation(size, channel_matrix)
            encr_channels.append(encr_channel)

def array_to_image(imageMatrix):
    """
    Convert list to PIL Image.
    """
    image_shape = np.array(imageMatrix).shape

    # if 3D
    if len(image_shape)==3:
        width, height, channels = image_shape
        im = Image.new("RGB", (width, height))
        pix = im.load()
        for x in range(width):
            for y in range(height):
                pix[x, y] = tuple(imageMatrix[x][y])
        return im

    # if 2D
    if len(image_shape)==2:
        width, height, channels = image_shape
        im = Image.new("L", (width, height))
        pix = im.load()
        for x in range(width):
            for y in range(height):
                pix[x, y] = imageMatrix[x][y]
        return im

    # else invalid
    return None

def encrypt(im):
    """
    Takes in a PIL image and returns an encrypted Henon Transformed Image as a PIL Image of shape (size, size).
    """
    # get image and shape
    pix = im.load()
    image_shape = np.array(im).shape

    # set size to min of height and width
    # the output encrypted image will be of shape (size, size)
    size = min(image_shape[0], image_shape[1])

    # if color image
    if len(image_shape) == 3:
        encr_channels = []
        # individually encrypt each channel
        for channel in range(3):
            channel_matrix = convert_to_array(im, channel)
            encr_channels.append(pixelManipulation(size, channel_matrix))
        encr_img = combine_channels(encr_channels)
        return array_to_image(encr_img)

    # if 2D grayscale image
    if len(image_shape) == 2:
        encr_img = pixelManipulation(size, convert_to_array(im))
        return array_to_image(encr_img)
    
    # else invalid image
    return None

def decrypt(im):
    """
    Takes in encrypted Henon Transformed Image as PIL Image and decrypts it, then returns the decrypted image as a PIL Image.
    """
    # get image and shape
    pix = im.load()
    image_shape = np.array(im).shape

    # set size to min of height and width
    # the output decrypted image will be of shape (size, size)
    size = min(image_shape[0], image_shape[1])

    # if color image
    if len(image_shape) == 3:
        decr_channels = []
        # individually decrypt each channel
        for channel in range(3):
            channel_matrix = convert_to_array(im, channel)
            decr_channels.append(decryptHenonImage(channel_matrix))
        decr_img = combine_channels(decr_channels)
        return array_to_image(decr_img)

    # if 2D grayscale image
    if len(image_shape) == 2:
        decr_img = decryptHenonImage(convert_to_array(encr_img))
        return array_to_image(decr_img)
    
    # else invalid image
    return None

def getImage(img_path):
    return Image.open(img_path)
