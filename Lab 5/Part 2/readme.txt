In this lab we will replace the single threshold aproach with Hysteresis&Double threshold approach in identifying the edges and compare the results obtained by the 2 algorithms.

Write the method part2() inside the l052.cpp file that:

a) reads the file image.ppm or the command line argument if there is one (a p3 format, size may be variable)

b) generates image1.ppm= a binary image after you applied the hysteresis algorithm (with double threshold) we discussed in class. You will use for low and high threshold the values that you picked for your image or the command line arguments if there are ones.

c) Hysteresis must be implemented using a recursive method.

d) Your main should have only a call to part2 (part1 is commented out).

----------In order to test your code wewill introduce the usage of parameters as follows:

1) you code will accept runtime parameters. Assume when you compile you obtain the executable l052, then:

     i) if I launch with no parameters like:

                 ./l052

       then your code will read the image.ppm and will use for the low and high threshold the values that you taylored for your

       image. If I run this and the image you pciked is in the same folder I should obtain the result you obtained 

       (your values are considered default values. So default value for input image is image.ppm, default values for low and high threshold values are the values you taylored to your image, the image with the edges you create is image1.ppm)

     ii) if I launch with parameters like this:

 

                 ./l052 -f myimg.ppm -lt 150 -ht 200 -of out.ppm

      then your code will instead read the file myimg.ppm, the ourput file will be out.ppm

      and will use for low and high threshold the values 150 and respctively

      200 (the parameters will be applied assuming that your code used for magnitude the complete formula. If you did not apply

     the square root than you will have to adjust in your code such that 150 and 200 are used as if you did apply the square root

   Note:: the command line parameters can be in any orders and if any of them is missing, the default values are used (see part i)) 

    

Also complete the attached document with the results (Project 5 Part 2 Canny Edge Detection Hysteresis.docx)

