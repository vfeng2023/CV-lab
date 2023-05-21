Read the following:

Chapter-16---The-three-dimensional-world_2018_Computer-Vision.pdf (you had to write an essay on it)

Chapter-17---Tackling-the-perspective-n-point-problem_2018_Computer-Vision.pdf (you had to write an essay on it)

Lookup homogeneouscoordinates,homography matrix / perspective projection matrix (Wikipediais a good source)

Platonic solids coordinates:PlatonicVertices.pdf

Your goal is to create a movie of a cube and after that of a tetrahdron rotating (for example rotating around z axis, or around any line; you may apply more than one transformation) and apply an orthographic projection to it.

Create the file l081.cpp which;

a) in the main you will call method part1 (with or without parameters) only

b) in part1 you will:

   1) will create the file rotation.avi with a cube rotating, then another platonic solid of your choise (Ex: tetrahedron) rotating using an ortographic projection

Just like i explained in class, all transformations are applied to the 3d coordinates of the vertices,  after each transformation, in order to obtain the 2d image in the video, you will apply an ortographic projection to the 3d coordinates and also you will shift the point in such a way that if I have for example the 3d point (0,0,0) the point will show up exactly in the center of the image). Also quartiles are exactly like in math, in other words point (100,100) will show up in quartile 1 which is the top right corner of your image. Each image in the video has the size 800x600.

   2) in each image you create, the (0,0) has to be in the center of the image, the size fo each image is 800x600 

   3) in the file coordinates.txt you will save the first 4 sets of 3D coordinates for the first 4 images in your video for the cube, each row containing the vertices for the cube in 3d like 

(1,1,1), (1,1,-1), (1,-1,1) , (1,-1,-1), (-1,1,1), (-1,1,-1), (-1,-1,1) , (-1,-1,-1)

....

Obs: your vertices have to be in the same order as above (first is the vertex corresponding to the vertex obtained after transforming (1,1,1) and so on

c) complete the following file: Project 8_1 Movie with rotating platonic solid with OpenCV.docx

 

Each frame showing the cube will be a wire frame cube. For vertices use a filled circle if possible, for each edge just connet the 2 vertices with a line.