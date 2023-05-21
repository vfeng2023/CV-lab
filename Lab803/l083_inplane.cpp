
#define _USE_MATH_DEFINES // for C++
#include <iostream>
#include <fstream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <cmath>
#include <stdlib.h>
#include <time.h> 
#include <unordered_set>


using namespace std;
using namespace cv;


//class Vector {
//private:
//	vector<double> mypoints;
//
//public:
//	Vector() {
//		mypoints = vector<double>(0.0, 3);
//	}
//	Vector(vector<double> v) {
//		mypoints = v;
//	}
//	Vector(const Vector& other) {
//		mypoints = vector<double>(other.mypoints);
//	}
//
//
//};
namespace vecops {

	double dot(Mat u, Mat v) {
		double res = 0;
		for (int i = 0; i < u.rows; i++) {
			res += u.at<double>(i, 0) * v.at<double>(i, 0);
		}
		return res;

	}
	double magnitude(Mat u, int col) {
		double size = 0;
		for (int i = 0; i < u.rows; i++) {
			size += u.at<double>(i,col) * u.at<double>(i,col);
		}
		return sqrt(size);
	}

	void scalarmult(Mat u, int col, double scalar) {
		//Mat res = u.cos
		for (int i = 0; i < u.rows; i++) {
			u.at<double>(i, col) *= scalar;
		}
	}
	Mat subtract(Mat& a, Mat& b) {
		// vector<double> res = vector<double>(a.rows,0);
		Mat res = Mat(a.rows, 1, CV_64F);
		for (int i = 0; i < a.rows; i++) {
			res.at<double>(i,0) = a.at<double>(i, 0) - b.at<double>(i, 0);
		}
		return res;
	}
};


/**
* @brief returns transformation matrix
* @param scaleFactor - amount to scale coordinates by
* @param x - degrees to rotate around xaxis
* @param y - degrees to rotate around yaxis
* @param z - degrees to rotate around zaxis
* @param transMatrix - matrix to output obtained transformation matrix in
*/
void findtransform(double scaleFactor, double x, double y, double z, Mat& transMatrix) {
	//multiply by scale factor
	//rotate by some amount
	//return the combined matrix
	double scale[4][4] = { {scaleFactor,0,0,0},
							{0, scaleFactor,0,0},
							{0,0,scaleFactor,0},
							{0,0,0,1} };
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
*@brief project single point onto the plane defined by a &n
* Precondition - 0 < plane < eye
* eye is location of eye, plane is location of plane for perspective projection
*/
void project(Mat v, Mat& result, Mat& a, Mat& normal, Mat& eye) {


	//use algebra owo

	//for (int i = 0; i < coordinates.cols; i++) {
		//Rect roi =
		//Mat v = coordinates(Range::all(), Range(i, i + 1));
		Mat ve = v-eye;
		Mat ae = a-eye;
		double t = vecops::dot(ae, normal) / vecops::dot(ve, normal);
		Mat x;
		ve.convertTo(result, -1, t);
		result += eye;
		//cout << "line 122\n";


	//}
	



}

/**
* @brief applies the gram schmit process to two vectors
*/
void gramschmit(Mat a, Mat b, Mat& c, Mat& d) {
	c = a / vecops::magnitude(a, 0);
	d = b - (c.dot(b) / c.dot(c)) * c;

	d /= vecops::magnitude(d,0);
}
/**
*@brief get the vectors p0, 01, and p2
* 
*/
void getBasis(Mat& coordinates, Mat& a, Mat& normal, Mat& eye, Mat &w1,Mat& w2,Mat &p0) {
	//Mat p0;
	double origin[3][1] = { {0},{0}, {0} };
	Mat center = Mat(3, 1, CV_64F, origin);
	//cout << center.rows << " " << center.cols << " size of center\n";
	project(center, p0, a, normal, eye);
	//cout << p0.rows << " " << p0.cols << " size of p0\n";

	Mat p1;
	//cout << coordinates.col(0).rows << " " << coordinates.col(0).cols << " size of p0\n";
	project(coordinates(Range(0, p0.rows), Range(0, 1)), p1, a, normal, eye);
	Mat p2;
	project(coordinates(Range(0, p0.rows), Range(1, 2)), p2, a, normal, eye);
	Mat p3;
	project(coordinates(Range(0, p0.rows), Range(2, 3)), p3, a, normal, eye);
	//create a and b
	Mat aprime = p1 - p2;
	Mat b = p3 - p2;

	//get normalized basis
	//Mat w1;
	//Mat w2;
	gramschmit(aprime, b, w1, w2);
	/*cout << "applied gram schmit\n";
	cout << w1.rows << " " << w1.cols << "size of w1\n";
	cout << vecops::magnitude(w1, 0) << vecops::magnitude(w2, 0) << endl;
	cout << w1.dot(w2) << "dot product w1 w2";*/
}
/**
* @brief project a single point
*/
void projectPoint(Mat& point, Mat &sol, Mat& a, Mat& normal, Mat& eye, Mat& w1, Mat& w2, Mat &p0) {
	//convert to u,v coordinates by solving systems of equations
	Mat pv;
	project(point, pv, a, normal, eye);
	Mat constcol = pv - p0;
	//combine the values into a matrix
	// 
	Mat eqns = Mat(pv.rows, 2, CV_64F);
	w1.copyTo(eqns.col(0));
	w2.copyTo(eqns.col(1));
	//constcol.copyTo(eqns.col(2));
	//perform row reduction
	/*cout << eqns.rows << " " << eqns.cols << "\n";
	cout << constcol.rows << " " << constcol.cols << "\n";*/
	//Mat sol;
	solve(eqns, constcol, sol, DECOMP_SVD);
	//cout << "solved " << endl;
}
/**
* @brief project the coordinates onto the plane
*/
void projectShape(Mat& coordinates, Mat& proj, Mat &p0, Mat &a, Mat &normal, Mat &eye, Mat &w1, Mat &w2) {
	//identify the points defining the plane in terms of origin (p0) and pv1 and pv2 and pv3
	//place in the matrix projected
	proj = Mat(coordinates.rows, coordinates.cols, CV_64F);
	for (int i = 0; i < coordinates.cols; i++) {
		//identify the projection coordinate
		Mat pv;
		project(coordinates(Range(0,p0.rows),Range(i,i+1)), pv, a, normal, eye);
		//convert to u,v coordinates by solving systems of equations
		Mat constcol = pv - p0;
		//combine the values into a matrix
		// 
		Mat eqns = Mat(pv.rows, 2, CV_64F);
		w1.copyTo(eqns.col(0));
		w2.copyTo(eqns.col(1));
		//constcol.copyTo(eqns.col(2));
		//perform row reduction
		/*cout << eqns.rows << " " << eqns.cols << "\n";
		cout << constcol.rows << " " << constcol.cols << "\n";*/
		Mat sol;
		solve(eqns, constcol, sol, DECOMP_SVD);
		//cout << "solved " << endl;
		proj.at<double>(0, i) = sol.at<double>(0, 0);
		proj.at<double>(1, i) = sol.at<double>(1, 0);
	}





}
//class edge to use in set 
class Edge {
private:
	int v1;
	int v2;
public:
	Edge() {
		v1 = 0;
		v2 = 0;
	}
	Edge(int ver1, int ver2) {
		v1 = ver1;
		v2 = ver2;
	}Edge(const Edge& other) {
		v1 = other.v1;
		v2 = other.v2;
	}
	int getv1() const{
		return v1;
	}
	int getv2() const{
		return v2;
	}
	bool operator==(const Edge &other) const{
		return (v1 == other.v1 && v2 == other.v2) || (v2 == other.v1 && v1 == other.v2);
	}

	class HashFunc {
	public:
		size_t operator()(const Edge& other) const {
			return (other.getv1()) ^ (other.getv2());
		}
	};


};


/**
* @brief make video
*/
void makeVideo(double scale, Mat& coordinates, vector<vector<int>>& edges, int rows, int cols, bool writeCoordinates, vector<Mat>& vidFrames, double plane, double eye) {
	//vector<Mat> vidFrames;
	//rotate around x and y axis, scale by a factor of 250, create 360 frames (rotate by 1 degree)
	//30 frames per second
	//double scale = 250;
	double x = 0;
	double y = 0;
	double z = 0;

	//define the normal parameters
	Mat a = (Mat_<double>(3, 1) << 0,plane,0);
	Mat normal = (Mat_<double>(3, 1) << 0, 1, 0);
	Mat eyepos = (Mat_<double>(3, 1) << 0, eye, 0);
	//vector<Mat> framevertices;
	int count = 0;
	ofstream coordsfile;
	ofstream coords2dfile;
	if (writeCoordinates) {
		coordsfile.open("coordinates.txt");
		coordsfile.precision(17);
		coords2dfile.open("coordinates2d.txt");
		coords2dfile.precision(17);
	}
	//make edges random colors
	vector<Scalar> edgeColor = vector<Scalar>();
	for (size_t q = 0; q < edges.size(); q++) {
		int r = rand() % 256;
		int g = rand() % 256;
		int b = rand() % 256;
		edgeColor.push_back(Scalar(b, g, r));
	}
	Mat w1, w2, p0;
	getBasis(coordinates, a, normal, eyepos, w1, w2, p0);
	for (int i = 0; i < 360; i++) {
		Mat tranformationMat;
		findtransform(scale, x, y, z, tranformationMat);
		//cout << tranformationMat.size()<<endl;
		Mat result = tranformationMat * coordinates;
		//cout << "line 71\n";
		Mat projected;
		
		projectShape(result, projected, p0, a, normal,eyepos,w1,w2);
		//create video frame
		//move points to center
		if (count < 4 && writeCoordinates) {
			for (int idx = 0; idx < result.cols; idx++) {
				double x = result.at<double>(0, idx);
				double y = result.at<double>(1, idx);
				double z = result.at<double>(2, idx);

				double u = projected.at<double>(0, idx);
				double v = projected.at<double>(1, idx);
				coordsfile << "(" << x << ", " << y << ", " << z << ")";
				coords2dfile << "(" << u << ", " << v << ")";
				if (idx < result.cols - 1) {
					coordsfile << " , ";
					coords2dfile << " , ";
				}
				else {
					coordsfile << "\n";
					coords2dfile << "\n";
				}
			}
			count += 1;
		}
		//for (int i = 0; i < projected.cols; i++) {
		//	projected.at<double>(0, i) += double(cols / 2);
		//	projected.at<double>(1, i) = double(rows / 2) - projected.at<double>(1, i);

		//}
		//draw vertices
		Mat frame = Mat(rows, cols, CV_8UC3, Scalar(255, 255, 255));;

		//if plane is straddling the figure
		vector<vector<Point>> projpoints = vector<vector<Point>>(coordinates.cols, vector<Point>()); // list of potential locations for vertices
		vector<vector<int>> projpointedges = vector<vector<int>>(coordinates.cols, vector<int>()); //list of corresponding endpoints for vertices
		//iterate through all edges
		for (int i = 0; i < edges.size(); i++) {
			//identify the projection vertices in two lists -- mod and the corresponding edge for the mode
			int e1 = edges[i][0];
			int e2 = edges[i][1];

			Mat v1 = result(Range(0, 3), Range(e1, e1 + 1));
			Mat v2 = result(Range(0, 3), Range(e2, e2 + 1));

			double dot1 = normal.dot(v1 - a);
			double dot2 = normal.dot(v2 - a);

			if (dot1 <= 0 && dot2 <= 0) {
				int x1 = cvRound(projected.at<double>(0, e1) + double(cols / 2));
				int y1 = cvRound(double(rows/2) - projected.at<double>(1, e1));

				int x2 = cvRound(projected.at<double>(0, e2) + double(cols / 2));
				int y2 = cvRound(double(rows / 2) - projected.at<double>(1, e2));
				projpoints[e1].push_back(Point(x1, y1));
				projpointedges[e1].push_back(e2);
				projpoints[e2].push_back(Point(x2, y2));
				projpointedges[e2].push_back(e1);
			}
			else if(dot1 > 0 and dot2 > 0){
				//do nothing
			}
			else {
				Mat left, right;
				if (dot1 < 0 && dot2 > 0) {
					left = v1;
					right = v2;
				}
				else {
					left = v2;
					right = v1;
				}

				//determine the intersection
				Mat u = right - left;
				double t = left.dot(normal) / u.dot(normal);
				Mat newright = t * u + left;
				if (dot1 < 0 && dot2 > 0) {
					
					projpoints[e1].push_back(Point(cvRound(projected.at<double>(0,e1))+(cols/2), (rows/2) - cvRound(projected.at<double>(1, e1))));
					projpointedges[e1].push_back(e2);
					Mat p2;
					
					projectPoint(newright, p2, a, normal, eyepos, w1, w2, p0);
					
					projpoints[e2].push_back(Point(cvRound(p2.at<double>(0, 0))+cols/2, cvRound(rows/2 - p2.at<double>(1, 0))));
					projpointedges[e2].push_back(e1);

				}
				else {
					projpoints[e2].push_back(Point(cvRound(projected.at<double>(0, e2))+cols/2, cvRound(rows/2 - projected.at<double>(1, e2))));
					projpointedges[e2].push_back(e1);
					Mat p1;
					projectPoint(newright, p1, a, normal, eyepos, w1, w2, p0);
					projpoints[e1].push_back(Point(cvRound(p1.at<double>(0, 0))+cols/2, cvRound(rows/2 - p1.at<double>(1, 0))));
					projpointedges[e1].push_back(e2);
				}
				

			}
		}
		//connect edges when there is a parallel one in the list
		//place the edges into a set
		unordered_set<Edge, Edge::HashFunc> edgeset = unordered_set< Edge, Edge::HashFunc>();
		vector<Edge> todraw = vector<Edge>();
		for (vector<int>& edge : edges) {
			edgeset.insert(Edge(edge[0], edge[1]));
		}
		
		/*for (int v = 0; v < coordinates.cols; v++) {
			line()
		}*/
		//iterate through list of edges
		//if both are len 1:
			//draw edge
		//otherwise:
			//skip
		//if both are greater than one:
			//add to list of edges to draw
		for (vector<int>& ed:edges) {
			if (projpointedges[ed[0]].size() == 1 && projpointedges[ed[1]].size() == 1) {
				line(frame, projpoints[ed[0]][0], projpoints[ed[1]][0],Scalar(0, 0, 0), 3);
			}
			else if (projpointedges[ed[0]].size() > 1 && projpointedges[ed[1]].size() > 1) {
				todraw.push_back(Edge(ed[0], ed[1]));
			}
		}
		//connect all multi-vertices to themselves
		for (int i = 0; i < projpointedges.size(); i++) {
			for (int j = 0; j < projpointedges[i].size(); j++) {
				//draw the vertex
				circle(frame, projpoints[i][j], 3, Scalar(0, 0, 0), -1);
				circle(frame, projpoints[i][j], 5, Scalar(0, 0, 255), 3);
				for (int k = j+1; k < projpointedges[i].size(); k++) {
					line(frame, projpoints[i][j], projpoints[i][j], Scalar(0, 0, 0), 3);
				}
			}
		}
		// go through list of to draw edges:
			//only connect edges that are in the edge set
		for (const Edge& t : todraw) {
			int v1 = t.getv1();
			int v2 = t.getv2();
			for (int i = 0; i < projpoints[v1].size(); i++) {
				for (int j = 0; j < projpoints[v2].size(); j++) {
					if (edgeset.find(Edge(projpointedges[v1][i], projpointedges[v2][j])) != edgeset.end()) {
						line(frame, projpoints[v1][i],projpoints[v2][j], Scalar(0, 0, 0), 3);
					}
				}
			}
		}
		x += M_PI / 180;
		y += M_PI / 180;
		vidFrames.push_back(frame);


	}
	if (writeCoordinates) {
		coordsfile.close();
		coords2dfile.close();
	}




}
//void checkSides(Mat &coordinates, Mat &a, Mat &normal, Mat & )
void part2() {
	//calculations
	//Mat coordinates;
	//vector<vector<int>> edges;
	//if (shape == 0) {
	Mat cubecoordinates = (Mat_<double>(4, 8) << 1, 1, 1, 1, -1, -1, -1, -1,
		1, 1, -1, -1, 1, 1, -1, -1,
		1, -1, 1, -1, 1, -1, 1, -1,
		1, 1, 1, 1, 1, 1, 1, 1);
	vector<vector<int>> cubeedges = { {0,1},{1,5},{1,3},{0,2},{0,4},{4,5},{4,6},{7,5},{6,7},{6,2},{3,2},{3,7} };

	//}
	//else {
	double sqrt3 = sqrt(3);
	Mat tetracoordinates = (Mat_<double>(4, 4) << 1 / sqrt3, 0, -sqrt3 / 6, -sqrt3 / 6,
		0, 0, 0.5, -0.5,
		0, 2 / sqrt(6), 0, 0,
		1, 1, 1, 1);
	vector<vector<int>> tetraedges = { {0,1},{0,3},{0,2},{1,3},{1,2},{2,3} };
	//}
	//making the frames
	vector<Mat> vidFrames;
	//makeVideo("tetrahedron.avi", coordinates, edges, 600, 800);
	makeVideo(150, cubecoordinates, cubeedges, 600, 800, true, vidFrames, 175, 650);
	makeVideo(250, tetracoordinates, tetraedges, 600, 800, false, vidFrames, 110, 650);
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
	part2();
}
