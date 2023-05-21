Your goal is to create a movie of a cube and after that of a tetrahdron rotating (for example rotating around z axis, or around any line; you may apply more than one transformation) and apply a perspective projection to it (simplified version where the screen is a plane defined by a vector a and a vector n as we discussed in class). The eye is defined by a vector such that the eye is on the other side of the planed you defined compared to the vertices of the platonic solid you rotate

Create the file l083.cpp which:

a) you will have 3 command line arguments as follows (they are optional; if any is missing ue the default values you picked):

  -a (x,y,z) == sets the value for vector a, which is a point on the viewing plane 

  -n (x,y,z) == sets the value of the normal vectoron the viewing plane

  -e (x,y,z) == sets the value of the eye vector, position of the eye

Example: ./l083 -a (700.5,200.1,100.12) -n (1,2,3.12) -e (1000,900.23,100.1156)

  will use for the viewing plane the points abobe and for the eye will use (1000,900.23,100.1156). Notice the values have comas but do not have spaces!!

a) in the main you will call method part3 (with or without parameters, maybe send command line arguments) only

b) in part3 you will:

   1) define the plane using 2 vectors a, n (this will be overriden if the command line arguments exist)

   2) define the position of the eye (on other side of the plane compared to the vertices of the platonic solid) (this will be overriden if the command line arguments exist)

   3) will create the file rotation.avi with a cube rotating, then another platonic solid of your choice rotating using a perspective projection using the plane and the eye you defined

(Just like in previous project,:

 - all transformations are applied to the 3D coordinates of the vertices (save the values for the 3d coordinates for the vertices for first 4 frames in coordinates.txt --only for first 4 frames)

- after each transformation, in order to obtain the 2D image in the video, you will apply a perspective projection to the 3D coordinates that was presented in class (save the coordinates of the vertices for first 4 frames in the file coordinates2d.txt)

- then you do the rendering, the creation of the actual image, so you will shift the point in such a way that if I have for example the 2D point (0,0) the point will show up exactly in the center of the image, or a point (100,100) will show up in the Q1 if you would divide the screen with axis like in math with the rigin in the center) Each image in the video has the size 800x600.

   2) in each image you create, the (0,0) has to be in the center of the image, the size fo each image is 800x600 

   3) in the file log.txt you will save all the information that define the problem, follow the exact format of the sample provided. The information you will include are (but not limited to): what is a, what is n, what is e, what is p0 and what is w1 and w2 (remember w1 p0 w2 are the 3 vectors that define the coordinate system in the plane defined by A and n, they depend on where the eye is and where the cube is), what were the first 4 sets of the edges in 3D coordinates ad 2D coordinates for the first 4 images in your video for the cube, each row containing the vertices for the cube in 3D like , see log.txt provided as a sample

 

c) complete the following file: Project 8_3 Movie with rotating platonic solid with OpenCV.docx

 

Each frame showing the cube will be a wire frame cube. For vertices use a filled circle, preferably red, for each edge just connect the 2 vertices with a line.

Due date: 05/22/2023, Monday