Part 3: create a method part3() that

Create and submit an application for the recursive, O(n*log^2 n) algorithm presented in class (a variation of the algorithm in the paper:closestPairRecursive.pdf).

Reading from the file points.txt and the creation of the vector should not be timed, only the algorithm should be timed.

The complete recursive algorithm (this assumes the vector is already filled with the points):
a) Sort the vector based on the x coordinate; O(nlogn). This step is outside the recursive method and should be done only once. Next steps describe the recursive method.
b) Divide the vector into 2 parts and recur for each part (logical division preferred since is faster); O(1) since is sorted
c) If there are 3 or 2 points just simply return the minimum distance since it can be done in constant time O(1)
d) When the 2 recursive calls return with 2 minimum distances and the points(let's say d1,p1,p2 and d2,p3,p4), calculate d = min(d1,d2) and sae the points that correspond to the minimum
e) create a strip of distance d (disyance d on x coordinate) to the left of the middle and d to the right of middle (middle value). Save all the points in the strip in a vector (let's call it vectorStrip, but you may use any name you wish)

f) sort the vectorStrip by y coordinate O(nlogn)

g) for each point in vectorStrip calculate the distance to max of 15 following points to see if you can find a smaller distance than d, then update the d and the points. This is O(15n)=O(n)

g) return the minimum distance you obtained and the 2 points that have that distance

In the main you should:
a) call part2(...)
b) call part3(..)
c) display on the screen and in the results.txt the 2 points and minimum distance obtained for both approaches also the time to complete each approach (you may do this either by creating global variables or by making both part2 and part3 to return some result)

Please turn in the cpp fle and the following file filled:

Project 3 Part 3.docx

************************************************

1) You may use for part 3 the following files to test (care for many points the brute force will be almost impossible to test)

points1k.txt

points10k.txt

points100k.txt

points1m.txt

2) You may also attempt (after completing 3.3) this extremely efficient algorithm:closestPairKhuller.pdf (this is actually very similar to 3.4 which we will implement in part 4)