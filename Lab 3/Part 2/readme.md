Create a copy of your previous lab and rename it l032.cpp

Part 2: create a method part2() that
Create and submit an application that solves the closest-pair problem in the unit square using the "brute force" approach discussed in class as well as the preliminary recursive algorithm.

For the recursive approach you must use at least one vector to store the points to have O(1) access to any element. Using an iterator with the vector is optional. Vectors give O() access to their elements based on the index.

In part2() method implement the following:

1) Before you start the preliminary recursive call you should read the points from the file points.txt into a vector (this part should not be timed).

2) Implement the preliminary recursive algorithm:

a) Sort the vector based on the x coordinate; O(nlogn). This step is outside the recursive method and should be done only once. Next steps describe the recursive method.

    (---recursive method - discussed in class)

b) Divide the vector into 2 parts and recur for each part (logical division preferred since is faster); O(1) since is sorted

c) If there are 3 or 2 points just simply return the minimum distance since it can be done in constant time O(1)

d) When the 2 recursive calls return with 2 minimum distances (let's say (d1,p1,p2) and (d2,p3,p4)), calculate d = min(d1,d2)

e) create a strip of distance d to the left of the middle and d to the right of middle (middle value)

f) for each point on the left of the strip calculate the distance to eachpoint on the right side of the strip and if you find a smaller distance than d, then update the d and the points that represent the closest pair of points. O(n^2) since it may have almost all points inside the strip (is almost a brute force..you may also brute force the strip but is a worse performance)

g) return the minimum distance you obtained and the 2 points that have that distance

 

In the main you should:

a) call part0()

a) call part1()

b) call part2()

c) display on the screen and in the results.txt the 2 points and minimum distance obtained for both approaches also the time to complete each approach (you may do this either by creating global variables or by making both part1 and part2 to return some result)

Please complete the document: Project 3 Part 2.docx

 

Obs:

- when you are timing the bruteforcedo not include the time to create the list, just the time to find the 2 closest points

- when you are timing the preliminary recursive method do not include the time to create the vector, but do include the time to sort the vector and to find the 2 closest points. (yes the time to sort is a price you pay for this approach , is the overhead)

- when you compare the 2 methods make sure you test also for more points (go up to 1000, or even more if you wish)

- testing on terminals is faster since some of the machines have gpu's and stronger processors