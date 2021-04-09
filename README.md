# 6.857 final project

## Image compression

The Instagram compression algorithm is proprietary and seems to be platform-dependent.  For better consistency and tunable parameters, we can use this ImageMagick command:

```
convert feldroy-512x512-unoptimized.jpg \
-sampling-factor 4:2:0 \
-strip \
-quality 85 \
-interlace Plane \
-gaussian-blur 0.05 \
-colorspace RGB \
feldroy-512x512-combined.jpg 
```
([source](https://dev.to/feldroy/til-strategies-for-compressing-jpg-files-with-imagemagick-5fn9))
