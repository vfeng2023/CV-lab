Lab name : l051.cpp

This lab requires you to implement the canny edge detection algorithm for an image stored in a ppm.

Please read he following document:

CannyEdgeDetection.pdf

Also the following 2 links will be useful to read:

https://en.wikipedia.org/wiki/Grayscale

https://en.wikipedia.org/wiki/Sobel_operator

---------------------------------------------------------------------

Using the convert tool offered by Linux:

You can use the convert utility to change a .ppm (P3) file to , say, .jpg format like this:

    convert infile.ppm outfile.jpg

However, if you try to convert a .jpg file to .ppm file by typing

    convert infile.jpg outfile.ppm

you will get a P6-formatted .ppm file, which you won't be able to read with a text-editor (other than the first line).

In order to convert a .jpg file to .ppm (P3) format (which can be inspected with a text editor), type

    convert -compress none infile.jpg outfile.ppm

----------------------

Create a file named l051.cpp in which you create part1() method that (your main has only a part1() call):

a) that reads the file image.ppm (a p3 format, size may be variable) and converts it into 2 files:
a) imageg.ppm = a gray version of the image implementing the algorithm discussed in class
b) imagem.ppm = an image after you applied the sobell operator and a threshold

Make sure your name is inside the cpp file you upload. Name your cpp file l051.cpp (lower case L and digits zero five one).

In the assignment upload the image you used (name it image.ppm), the cpp file you created and the filled document below.
Your code may be run against other images for testing, so you should not hard code the size of the image since that information is already in the initial ppm file

The document to fill for this project is:

 Project 5 Part 1 Canny Edge Detection with Sobel and single threshold.docx