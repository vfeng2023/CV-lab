Name offile l061.cpp

Your main should call only the folllowing method: part1(...) (with or without parameters.

Create part1(...) method that does the following:

1) read the file image.ppm (image will contain an image of coins)

2) create the file imagef.ppm (and all intermediary ppms from 5.3) by applying the complete canny edge detection from previous lab (optimize the values to work best for the easy file I provided) NOT Optional: you also create the imageg.ppm. image1.ppm and image2.ppm for debugging purposes.

3) use the edges you found and gradient direction (not rounded!!!) to implement voting and create the imagev.ppm that displays the result of the phase 1 voting for centers.

4) use the results of voting and a threshold value to pick good candidates for circle centers and create imageCC.ppm file that will display the original image and on the original image will do a filled circle of radius 4 of color red for each candidate for circle center. (a filled circle of radius 4 can be obtained by drawing 4 circles with that center of radius 1,2,3,4)

5) allow the tester to change the lower and upper threshold by using command line arguments as follows:

Ex1: 

  ./l061 -f myimg.ppm -lt 150 -ht 200 -fg mygray.ppm -f1 myimage1.ppm -f2 myimage2.ppm -ff myimagef.ppm -TC 300 -fv myimagev.ppm -fcc myimageCC.ppm

        - notice only the last 3 parameters are new (all other ones are from previous project 5.3)

        - the line above will use 150 and 200 for low and high threshold for th canny edge and the 300 will be the threshold value used for deciding if a pixel is a candidate for a center or not (so if a pixel received >=300 votes it means is a high chance that the pixel is a center) Also myimagev.ppm will be the image for the votes and myimageCC.ppm will be the image for the canddates for center

Ex2: ./l061

The line above has no parameters which means it will use the values for lower and higher threshold that you optimized for your image, also will use the threshold value you optimized for the first image and it will apply it to the file image.ppm.

***********************************************************

 

Submit the following document after your code runs on the easy image I provided (which you will transform in p3 ppm file using the convert tool on gnu linux):

Project 6 Coin Detection Part 1.docx

Your code will be tested against any piece of the easy image provided not just against the entire image.