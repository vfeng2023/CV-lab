Rename your lab l021.cpp into l022.cpp.

Create a method part2(), in the main you should comment out calling part1, and just call part2 method.

In method part2 you should:
1) Read 4 points from the file points.txt; (sample you can use for testing: points.txt, do not hard code the precision, but use spaces, comas and parentheses to delimitate the coordinates)

2) Then create all possible squares that have the 4 points on its side or the extension of its sides.
You may assume the problem is well defined (or you may deal with this situation too if you wish).

3) Create the file output.ppm of size 800x800 in which:

  a) Display the 4 points by creating a circle for each point with radius = 3 and center on each point

  b) Draw the square with the minimum area with all sides extended. (for the sample above, you should obtain: output.ppm)


4) Save in a file called output.txt the coordinates of the 4 points and the coordinates of the 4 vertices for ALL the squares you found. it should have the format of the file sampleoutput.txt 

(for the sample above, you should obtain:output.txt , but your order of the squares may differ)

To find the squares you may use any of the solutions presented at
https://www.cut-the-knot.org/Curriculum/Geometry/GeoGebra/SquareFromFourPoints.shtml
or
http://kirkmcd.princeton.edu/examples/4point.pdf

You will upload the following:

1) l022.cpp you created (remember to comment out part1 call in the main so you read the points from the file)

2) a completed document of:Project 2 Part 2 document.docx