The goal is to generate 60 random points in the unit square, then create a convex hull using 2 algorithms.

Here are the links to the 2 algorithms that need to be implemented for lab 4 to obtain a convex hull:

Part 1: https://en.wikipedia.org/wiki/Quickhull

Part 2: https://en.wikipedia.org/wiki/Graham_scan

Once the convex hull is obtained create the quickhull.ppm image in which you will draw a circle of radius 3 (bold) for each point and also will draw the convexhull by connecting the points that form the convex hull.

Part 1 is due Monday 12/05

 a) name your file l041.cpp (lower case L followed by the digits 0,4,1)

 b) in the main you should have only a call to part1 method (with arguments or not, is your choice)

 c) generate 60 random points in the unit square (save them in a data structure of your chosing) and save them in points.txt in the same format as 3.4 (each line has the x and y coordinate of a point separated by 2 spaces, use at least 20 digits precision)

 d) apply the QuickHull algorithm to find a convex hull that encloses all 60 points and the vertices are points from the set of 60 points you generated

 e) create a ppm of size 400x400 named quickhull.ppm in which you display a circle of radius 3 for each scaled point you generated, as well as connect the vertices you obtained for the convex hull (so in the image I should see 60 points and a convex hull with vertices circles of radus 3 that are conected forming a convex hull)

  f) complete the attached form Project 4 Part 1 QuickHull.docx and turn in this document as well as the cpp file

Part 2 is due Monday 12/12

--details soon