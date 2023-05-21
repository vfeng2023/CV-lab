Write the method part3() inside the l053.cpp file that reads the file image.ppm (a p3 format, size may be variable) and generates:

a) greyscale image called imageg.ppm

b) also it generates 3 binary files:

    1) image1.ppm= a binary image after you applied the hystheresis algorithm (double threshold) we discussed in class

    2) image2.ppm = a binary version of the image implementing the non-maximum suppression algorithm we discussed in class

    3) imagef.ppm = a binary image after both algorithms above are applied

c) your main should only call part3()


Make sure your name is inside the cpp file you upload. Name your cpp file l053.cpp (lower case L and numbers 0, 5 and 3).
Upload here the image you used (name it image.ppm) and the cpp file you created. Your code may be run against other images for testing, so you should not hard code the size of the image since that information is already in the initial ppm file

----------In order to test your code we will continue the usage of parameters as follows:

1) you code will accept runtime parameters. Assume when you compile you obtain the executable l053, then:

     i) if I launch with no parameters like:

                 ./l053

       then your code will read the image.ppm and will use for the low and high threshold the values that you taylored for your

       image. If I run this and the image you pciked is in the same folder I should obtain the result you obtained 

       (your values are considered default values. So default value for input image is image.ppm, default values for low and high threshold values are the values you taylored to your image, the image with the edges you create is image1.ppm)

     ii) if I launch with parameters likethis:

 

                 ./l053 -f myimg.ppm -lt 150 -ht 200 -fg gray.ppm -of out.ppm -f2 myimg2.ppm -ff myFinal.ppm

      then your code will use the comand line arguments above as follow:

         a) myimg.ppm will be used instead of image.ppm (-f option)

         b) for low and high threshold the values 150 and 200 will be used (-lt and -ht options)

         c) out.ppm will be created instead of image1.ppm (-of option)

         d) gray.ppm will be created instead of imageg.ppm (-fg option)

         e) myimg2.ppm will be created instead of image2.ppm ( -f2 option)

         f) myFinal.ppm is created instead of imagef.ppm (-ff option)

      (The low and high threshold parameters will be applied assuming that your code used for magnitude the complete formula. If you did not apply the square root than you will have to adjust in your code such that 150 and 200 are used as if you did apply the square root)

   Note:: the command line parameters can be in any orders and if any of them is missing, the default values are used (see part i)) 

Turn in the following document after you fill it:

Project 5 Part 3 Canny Edge Detection complete.docx

Project due Friday 01/20/2023