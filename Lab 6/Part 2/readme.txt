Name of file l062.cpp

Your main should call only the folllowing method part2(...). (feel free to copy paste part1 or call part1 in part2)

Create part2(...) method that does the following:

1) read the file image.ppm (image will contain an image of coins). make yoru code work on the easy image I provided

2) create the file imagef.ppm by applying the complete canny edge detection from previous lab (optimize the values to work best for the easy file I provided, rest of files may or may not be created by your code, is your choice, since not creating ppm files will boost the speed)

3) use the edge you found and gradient direction (not rounded) to implement voting and create the imagev.ppm that displays the result of the phase 1 voting for centers (.this is part 6.1)

4) use the results of voting and a threshold value to pick good candidates for circle centers and create imageCC.ppm file that will display the original image and on the original image will do a filled circle of radius 4 of color red for each candidate for circle center. (a filled circle of radius 4 can be obtained by drawing 4 circles with that center of radius 1,2,3,4)

5) For each candidatefor circle center count the edges that are at distance rmin, then at distance rmin+1 all the way to rmax (where rmin and rmax are values you decide to be the min radius of a coin and max radius of a coin for the image. hese values can be taylored for the image, or may be a percentage of the image)

6) Set a threshold (a percentage of total pixels; in the process of your application you may use a set value to start coding/testing) of when a certain count/percentage is a circle  and obtain a list of centers and radii.

7) implement some algorithm that tries to eliminate duplicates circles (For ex: if 2 circles have the center very close by and their radiuses also very close, than keep he one with a higher count/percentage only and discard the other one). You may use any approach you wish to eliminate duplicates. How you define very close is your choice.

8) Create imageCircles.ppm that displays over the original image for each circle found (center = (xc,yc), radius=r} a filled circle of radius 4 with he center (xc,yc) for he center and a bold circle of radius r (can be done by drawing multiple circles with raduses for example r, r+1, r+2 with the center (xc,yc). Both should have the color red. Check visually if you ar happy with your results, if not attempt to fix this part.

8) Once you have the list of (center, radius) of circles you found, implement an algorithm to classify each circle as a penny, nickel, dime, quarter, silver dollar, half dollar based on the radius value and maybe color of pixels inside the coin (use as many pixels from inside the coin as you consider you need).

8) create the file imageCoins.ppm obtained by drawing over the original image

  a) a filled red circle for each center

  b) and a bold outline of the coins you found by using the following colors: red for penny, yellow for nickel, blue for dime, purple for quarters, green for half dollar and pink for silver dollar.. So draw a circle around each coin detected using the correct color.

  Check visually if you are happy with this result, if not go back and change parameters or fix your code.

9) display on the screen and in the file results.txt a summary of the coins you found ( Ex: 10 quarters, 5 dimes,... Total sum: $10.50)

10) allow the tester to change the lower and upper threshold by using command line arguments, but also the value of the threshold you are using for candidates for center and also the threshold value you are using for the step of finding centers as follows:

A) Ex1: ./l062 -f myimg.ppm -lt 150 -ht 200 -ff myimagef.ppm -TC 300 -fv myimagev.ppm -fcc myimageCC.ppm -TCircle 55 -fCi myimageCircles.ppm -fCo myimageCoins.ppm -fR myresults.txt

        - the line above will use 110 and 200 for low and high threshold for te canny edge; the 300 will be the threshold value used for deciding if a pixel is a candidate for a center or not (so if a pixel received >=300 votes it means is a high chance that the pixel is a center), and the thhreshold value of 55 which is 55% is used to identify if a candidate for a center and a certain radius is a circle or not. Rest of parameters allow the customization of the output files' names and they are self explanatory. (in case you use other variable values feel free to add them at the end in a similar fashion and specify them in the document so I know how to use them) Any or all the command line parameters may be missing.

B) Ex2: ./l062

The line above has no parameters which means it will use the values for lower and higher threshold that you optimized for the easy image, also will use the threshold value you optimized for the first image and it will apply it to the file image.ppm. The output file names are the oe specified above and are set.

 

11) complete the document below.

You must run the code on the easy images I provided, and your goal is to optimize the parameters in such a way to identify most/all coins and identify them correctly.

Accuracy of this algorithm depends greatly on the results of your canny edge detection algorithm.

The speed of this algorithm depends on 2 things;

  - the number of candidates for center you obtained in part1

  - the range you are using for rmin and rmax... the bigger this interval, the slower your algorithm is

Hints: use Bresenham'salgorithm to vote on the direction of the radiant. Also a good improvement would be to calculate the intersection points on the extremes so the line is more close to reality compared to using 2 very close points in which case the line will be off.

Submit the following document after your code runs on the easy image I provided (which you will transform in p3 ppm file using theconvert tool on gnu linux:

Project 6 Coin Detection Part 2.docx

Your code will be tested against any piece of the easy image.

 

Hint: while you code keep in mind in 6.3 you will have to scale your code to work on all 3 images. So while coding keep this in mind.

DUE DATE Monday 02/27/2023

If you finish the project faster, try to make your code  work against the other 2 images by only modifying the values in the command line arguments (they have to be the ones I specified + maybe some you created in your implementation). Maybe you need extra arguments?