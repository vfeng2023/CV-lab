Turn in the file called l071.cpp,for reading in a .jpg file of a collection of coins and outputting the dollar value of the coins.

In order to visualize and test the result of your application do the following:
- read from the file called image.jpg, create the grayscale imageg.jpg, detect edges using the cannyedge and create the imagef.jpg (feel free to use any tools that OpenCV offeres you in order to obtain better results)

- use command line parameters like before with whatever parameters you need in such a way that I can compile one time your code and run your code with different parameters for the 3 images i provided (you decide what are the parameters and how you write them, since you will write it in the document you attach to your project) (one of the parameters MUST BE the name of the file so I can run your code on any file I want!!)

- identify the centers and radiuses using OpenCV and hough transform (feel free to use any tools that OpenCV offeres you in order to obtain better results)

- create imageCircles.jpg : the original image with the circles you identified drawn in red over the image (you draw 4 circles with radiuses r, r+1, r+2, r+3 to make it more bold and easier to see on the image).

 - create the file imageCoins.ppm obtained by drawing over the original image the correct color for each coin, after you identified each coin. Use different colors for the coins you identified (the same colors described in project 6.2). You aim to have a circle drawn for each coin. (Again draw the circles bold by drawing a r, r+1, r+2 and r+3 using the correct color for the coin you identified)

- display on the screen and also save in the file results.txt a summary of the coins you found (Ex: 10 quarters, 5 dimes,... Total sum: $10.50)


After running your application 5 files will be created imageg.jpg, imagef.jpg, imageCircles.jpg , imageCoins.jpg and results.txt containing the information explained above

Also complete the following file:

Project 7 Coin Detection with OpenCV.docx

and display in it the results you obtained by applying your algorithm to all 3 images provided.

Your code must work on all 3 images.

You may test your code against any piece/part of the 3 images

Tentative due date: 03/20/2022.