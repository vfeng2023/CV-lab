

Your goal is to create a movie of a cube and after that of another platonic solid of our choice rotating (for example rotating around z axis, or around any line; you may apply more than one transformation) and apply a perspective projection to it (simplified version where the screen is a plane defined by x=a_number and the eye is on the other side of the viewing plane compared to the vertices). The viewing plane should be picked such that the rotating object's vertices are always on teh ame side of the viewing plane and the eye is always on the other side.

Create the file l082.cpp which;

a) in the main you will call method part2 (with or without parameters) only

b) in part2 you will:

   1) will create the file rotation.avi with a cube rotating, then another platonic solid rotating using a perspective projection

(Just like in previous project, all transformations are applied to the 3D coordinates of the vertices,  after each transformation, in order to obtain the 2D image in the video, you will apply an perspective projection to the 3D coordinates and also you will shift the point in such a way that if I have for example the 2D point (0,0) the point will show up exactly in the center of the image).  Also quartiles are exactly like in math, in other words point (100,100) will show up in quartile 1 which is the top right corner of your image. Each image in the video has the size 800x600.

   2) in each image you create, the (0,0) has to be in the center of the image, the size fo each image is 800x600 

   3) in the file coordinates.txt you will save the first 4 sets of 3D coordinates for the first 4 images in your video for the cube, each row containing the vertices for the cube in 3d like 

(1,1,1), (1,1,-1), (1,-1,1) , (1,-1,-1), (-1,1,1), (-1,1,-1), (-1,-1,1) , (-1,-1,-1)

....

Obs: your vertices have to be in the same order as above (first is the vertex corresponding to the vertex obtained after transforming (1,1,1) and so on

  4) in the file coordinates2d.txt you will save the first 4 sets of 2d coordinates for the first 4 frames in your video. The 2d coordinates are obtained after applying the perspective rendering to the 3d coordinates.

Obs: your vertices should match the vertices order from your coordinates.txt

....//rest of lines in similar fashion

c) complete the following file: Project 8_2 Movie with rotating platonic solid with OpenCV.docx

 

Each frame showing the cube will be a wire frame cube. For vertices use a filled circle so is visibleeven if the vertex ends up on an edge, for each edge just connect the 2 vertices with a line/segment. You may use different colors for vertices/edges such that when it rotates you see it goes at the far end.

Project 8.2 is due Monday 05/01. The due date was extended to Monday so that students who are ahead have a chance to complete the 2 extensions, while students who need more time to ctahc up have a chance to do so.

Optional possible improvements:

1) create a rotation that for some frames some edges or vertices go out of the viewing screen (still all ertices are on the other side of the viewing plane compared to the eye)

2) ccreate such a rotation that some vertices cross over the viewing plane and end up on the same side with the eye.