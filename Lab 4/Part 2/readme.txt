The goal is to create a convex hull using Graham Scan algorithm.

Part 2: https://en.wikipedia.org/wiki/Graham_scan

Once the convex hull is obtained create the grahamscan.ppm image in which you will draw a circle of radius 3 for each point and also will draw the cconvex hull by connecting the points that form the convex hull.

Part 2 is due Monday 12/12

 a) name your file l042.cpp (lower case L followed by the digits 0,4,2)

 b) in the main you should have only a call to part2 method (with arguments or not, is your choice). Make sure part1 is commented out!

 c) read the points from points.txt that has the same format as 3.4 (each line has the x and y coordinate of a point separated by 2 spaces, use at least 20 digits precision)

 d) apply the GrahamScan algorithm to find a convex hull that encloses all the points read from the file; the vertices are points from the set of the points you read

 e) create a ppm of size 400x400 named grahamscan.ppm in which you display a circle of radius 3 for each scaled point you generated, as well as connect the vertices you obtained for the convex hull (so in the image I should see 60 points and a convex hull with vertices circles of radius 3 that are conected forming a convex hull)

 f) complete the attached form Project 4 Part 2 Graham Scan Hull.docx and turn in this document as well as the cpp file