
#define _USE_MATH_DEFINES // for C++
#include <iostream>
#include <fstream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <cmath>



using namespace std;
using namespace cv;

/**
* @brief returns transformation matrix
* @param scaleFactor - amount to scale coordinates by
* @param x - degrees to rotate around xaxis
* @param y - degrees to rotate around yaxis
* @param z - degrees to rotate around zaxis
* @param transMatrix - matrix to output obtained transformation matrix in
*/
void findtransform(double scaleFactor, double x, double y, double z, Mat &transMatrix) {
	//multiply by scale factor
	//rotate by some amount
	//return the combined matrix
	double scale[4][4] = { {scaleFactor,0,0,0},
							{0, scaleFactor,0,0},
							{0,0,scaleFactor,0},
							{0,0,0,1}};
	double rotX[4][4] = { {1,0,0,0},
						{0, cos(x),-sin(x),0},
						{0,sin(x),cos(x),0},
						{0,0,0,1} };
	double rotY[4][4] = { {cos(y),0,sin(y),0},
						{0, 1,0,0},
						{-sin(y),0,cos(y),0},
						{0,0,0,1} };
	double rotZ[4][4] = { {cos(z),-sin(z),0,0},
						{sin(z), cos(z),0,0},
						{0,0,1,0},
						{0,0,0,1} };
	Mat scaleMatrix = Mat(4, 4, CV_64F, scale);
	Mat rotXMatrix = Mat(4, 4, CV_64F, rotX);
	Mat rotYMatrix = Mat(4, 4, CV_64F, rotY);
	Mat rotZMatrix = Mat(4, 4, CV_64F, rotZ);
	transMatrix = rotZMatrix * rotYMatrix * rotXMatrix * scaleMatrix;
	//cout << "Error 42\n";
}
/**
*@brief project coordinates
*/
void project(Mat &coordinates, Mat &projected) {
	//do projections
	Mat projMatrix = (Mat_<double>(4, 4) << 1, 0, 0, 0, 
											0, 1, 0, 0,
											0, 0, 0, 0,
											0, 0, 0, 1);
	projected = projMatrix * coordinates;
	//cout << "Error 54\n";
}
/**
* @brief make video
*/
void makeVideo(double scale, Mat& coordinates, vector<vector<int>> &edges, int rows, int cols, bool writeCoordinates, vector<Mat> &vidFrames) {
	//vector<Mat> vidFrames;
	//rotate around x and y axis, scale by a factor of 250, create 360 frames (rotate by 1 degree)
	//30 frames per second
	//double scale = 250;
	double x = 0;
	double y = 0;
	double z = 0;
	//vector<Mat> framevertices;
	int count = 0;
	ofstream coordsfile;
	if (writeCoordinates) {
		coordsfile.open("coordinates.txt");
		coordsfile.precision(17);
	}
	
	for (int i = 0; i < 360; i++) {
		Mat tranformationMat;
		findtransform(scale, x, y, z, tranformationMat);
		//cout << tranformationMat.size()<<endl;
		Mat result = tranformationMat * coordinates;
		//cout << "line 71\n";
		Mat projected;
		project(result, projected);
		//create video frame
		//move points to center
		if (count < 4 && writeCoordinates) {
			for (int idx = 0; idx < result.cols; idx++) {
				double x = result.at<double>(0, idx);
				double y = result.at<double>(1, idx);
				double z = result.at<double>(2, idx);
				coordsfile << "(" << x << ", " << y << ", " << z << ")";
				if (idx < result.cols - 1) {
					coordsfile << " , ";
				}
				else {
					coordsfile << "\n";
				}
			}
			count += 1;
		}
		for (int i = 0; i < projected.cols; i++) {
			projected.at<double>(0, i) += double(cols / 2);
			projected.at<double>(1, i) = double(rows/2) - projected.at<double>(1, i);

		}
		//draw vertices
		Mat frame = Mat(rows, cols, CV_8UC3, Scalar(255, 255, 255));;
		for (int p = 0; p < result.cols; p++) {
			int xval = cvRound(projected.at<double>(0, p));
			int yval = cvRound(projected.at<double>(1, p));
			circle(frame, Point(xval, yval), 3, Scalar(0, 0, 0), -1);
			circle(frame, Point(xval, yval), 5, Scalar(0, 0, 255), 3);
			

		}
		//cout << "line 92";
		//connect edges
		for (int i = 0; i < edges.size(); i++) {
			int e1 = edges[i][0];
			int e2 = edges[i][1];
			int x1 = cvRound(projected.at<double>(0, e1));
			int y1 = cvRound(projected.at<double>(1, e1));

			int x2 = cvRound(projected.at<double>(0, e2));
			int y2 = cvRound(projected.at<double>(1, e2));
			line(frame, Point(x1, y1), Point(x2, y2), Scalar(0, 0, 255), 3);
			

		}
		x += M_PI/180;
		y += M_PI/180;
		vidFrames.push_back(frame);
		

	}
	if (writeCoordinates) {
		coordsfile.close();
	}
	
	


}
void part1(int shape) {
	//calculations
	//Mat coordinates;
	//vector<vector<int>> edges;
	//if (shape == 0) {
		Mat cubecoordinates = (Mat_<double>(4, 8) << 1,1,1,1,-1,-1,-1,-1,
													 1,1,-1,-1,1,1,-1,-1,
													 1,-1,1,-1,1,-1,1,-1,
													 1, 1, 1, 1,1,1,1,1);
		vector<vector<int>> cubeedges = { {0,1},{1,5},{1,3},{0,2},{0,4},{4,5},{4,6},{7,5},{6,7},{6,2},{3,2},{3,7} };

	//}
	//else {
		double sqrt3 = sqrt(3);
		Mat tetracoordinates = (Mat_<double>(4,4) << 1/sqrt3, 0, -sqrt3/6,-sqrt3/6,
											0,       0, 0.5,     -0.5,
											0,	2/sqrt(6),0, 0,
											1, 1, 1 ,1);
		vector<vector<int>> tetraedges = { {0,1},{0,3},{0,2},{1,3},{1,2},{2,3} };
	//}
	//making the frames
	vector<Mat> vidFrames;
	//makeVideo("tetrahedron.avi", coordinates, edges, 600, 800);
	makeVideo(150, cubecoordinates, cubeedges, 600, 800, true, vidFrames);
	makeVideo(250, tetracoordinates, tetraedges, 600, 800, false, vidFrames);
	//write to file
	VideoWriter outputvid;
	int codec = VideoWriter::fourcc('M', 'J', 'P', 'G');  // select desired codec (must be available at runtime)
	double fps = 12;                          // framerate of the created video stream
	string filename = "rotation.avi";             // name of the output video file
	outputvid.open(filename, codec, fps, vidFrames[0].size(), true);
	for (int c = 0; c < vidFrames.size(); c++) {
		outputvid << vidFrames[c];
	}
	outputvid.release();
}
int main() {
	part1(1);
}
