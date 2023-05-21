The next lab is to implement the following algorithm to find the closest pair of points:

ClosestPairRandomized.pdf

Also a visual simulation Created by Dr. Torbert:

Closest Pair.mp4

The reading from the file and creation of the vectors used in part3 and part4 should not be included in the timing part.

In part 4:

a) shuffle the vector using the knuth shuffle we discussed in class

b) Complete the randomized algorithm described in the document above

c) Make sure your timing for part4 includes the shuffle and the randomized algorithm

 

In the main you should:
a) call part3(...)
b) call part4(..)
c) display on the screen and in the results.txt the 2 points and minimum distance obtained for both approaches also the time to complete each approach (you may do this either by creating global variables or by making both part3 and part4 to return some result)

Please turn in to Mr. Jurj a printout of the following file filled:

Project 3 Part 4.docx

*************************************************

A) You may use the following files to test:

points1k.txt

points10k.txt

points100k.txt

points1m.txt

B) use the following link if you need help to make an unordered map of pairs:

https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-user-defined-class-in-cpp/

https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-pairs-in-c/

 

C) for the pair (x,y) of the grid use "unsigned long long" so you don't run over the limit of int