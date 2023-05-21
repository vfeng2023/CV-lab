
#define _USE_MATH_DEFINES // for C++
#include <iostream>
#include <fstream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <cmath>
#include <stdlib.h>
#include <time.h> 
#include <random>


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
		//cout << "eye: " << eye << endl;
		Mat ve = v-eye;
		//cout << "ve: " << ve << endl;
		Mat ae = a-eye;
		//cout << "ae: " << ae << endl;
		double t = ae.dot(normal) / ve.dot(normal);
		/*Mat x;
		ve.convertTo(result, -1, t);
		result += eye;*/
		result = t * ve + eye;
		//cout << result << " result" << endl;
		//cout << "line 122\n";
		/*int pause;
		cin >> pause;*/

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
void getBasis(Mat& coordinates, Mat& a, Mat& normal, Mat& eye, Mat &w1,Mat& w2,Mat &p0, bool writetofile, ofstream &logfile) {
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
	Mat b = p3 - p1;

	//get normalized basis
	//Mat w1;
	//Mat w2;
	gramschmit(aprime, b, w1, w2);
	/*cout << "applied gram schmit\n";
	cout << w1.rows << " " << w1.cols << "size of w1\n";
	cout << vecops::magnitude(w1, 0) << vecops::magnitude(w2, 0) << endl;
	cout << w1.dot(w2) << "dot product w1 w2";*/

	if (writetofile) {
		logfile << "\nVertices I used to create he 2d coordinate system and their projections are: \n";
		for (int i = 0; i < 3; i++) {
			double x = coordinates.at<double>(0, i);
			double y = coordinates.at<double>(1, i);
			double z = coordinates.at<double>(2, i);
			logfile << "v" << (i + 1) << "=" << "(" << x << ", " << y << ", " << z << ")\n";
		}

		logfile << "\npv1 = " << "(" << p1.at<double>(0, 0) << ", " << p1.at<double>(1, 0) << ", " << p1.at<double>(2, 0) << ")\n";
		logfile << "pv2 = " << "(" << p2.at<double>(0, 0) << ", " << p2.at<double>(1, 0) << ", " << p2.at<double>(2, 0) << ")\n";
		logfile << "pv3 = " << "(" << p3.at<double>(0, 0) << ", " << p3.at<double>(1, 0) << ", " << p3.at<double>(2, 0) << ")\n";

		logfile << "\nThe 2 vectors a and b that are in plane are: \n";
		logfile << "a = pv1-pv2 = " << "(" << aprime.at<double>(0, 0) << ", " << aprime.at<double>(1, 0) << ", " << aprime.at<double>(2, 0) << ")\n";
		logfile << "b = pv3-pv1 = " << "(" << b.at<double>(0, 0) << ", " << b.at<double>(1, 0) << ", " << b.at<double>(2, 0) << ")\n";
	}
}
/**
* @brief project the coordinates onto the plane
*/
void projectShape(Mat& coordinates, Mat& proj, Mat &p0, Mat &a, Mat &normal, Mat &eye, Mat &w1, Mat &w2) {
	//place in the matrix projected
	proj = Mat(coordinates.rows, coordinates.cols, CV_64F);
	for (int i = 0; i < coordinates.cols; i++) {
		//identify the projection coordinate
		Mat pv;
		project(coordinates(Range(0,p0.rows),Range(i,i+1)), pv, a, normal, eye);
		//convert to u,v coordinates by solving systems of equations
		//cout << "value of pv: " << pv << endl;
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
		//cout << "solved " <<sol << endl;
		proj.at<double>(0, i) = sol.at<double>(0, 0);
		proj.at<double>(1, i) = sol.at<double>(1, 0);
		/*int pause;
		cin >> pause;*/
	}





}

/**
* @brief check if the vector is to the left of the plane
*/
bool isleft(Mat v, Mat& a, Mat& normal) {
	return normal.dot(v - a) < 0;
}

/**
* @brief return a random number
*/
int randNumber(int a, int b) {
	return rand() % (b - a+1) + a;
}
/**
* Returns true if vertices are all same sign, eye is opposite sign in terms of dot product
*/
bool overlap(vector<Mat> transformationCoords, Mat &a, Mat &eye, Mat &normal) {
	for (int i = 0; i < transformationCoords.size(); i++) {
		Mat& coords = transformationCoords[i];
		//Mat &v = coords(Range(0, 3), Range(i, i + 1))
		bool sign = isleft(coords(Range(0, 3), Range(0, 1)), a, normal);
		for (int c = 1; c < coords.cols; c++) {
			bool mysign = isleft(coords(Range(0, 3), Range(c, c + 1)), a, normal);
			if (mysign != sign) {
				return true;
			}
		}
		bool eyesign = isleft(eye, a, normal);
		if (eyesign == sign) {
			return true;
		}

		
	}
	return false;
}
/**
* @brief generate a, normal, and eye to not overlap with the shape itself
*/
void setRandom(vector<Mat> transformationCoords, Mat& a, Mat& eye, Mat& normal) {
	//coordinates of a between 250 and 450
	//coordinates of eye between 450 and 560
	//normal must be between 0 and 1
	/*std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;
	double a_random_double = unif(re);*/
	do {
		//generate coordinates for a
		std::uniform_real_distribution<double> unif(250, 450);
		random_device rd;
		mt19937 re(rd());
		for (int r = 0; r < a.rows; r++) {
			
			//double a_random_double = unif(re);
			a.at<double>(r, 0) = unif(re);
		}
		cout << "point: " << a << endl;
		uniform_real_distribution<double> eyefinder(450, 650);
		for (int r = 0; r < eye.rows; r++) {
			eye.at<double>(r, 0) = eyefinder(re);
		}
		cout << "eye: " << eye << endl;
		uniform_real_distribution<double> nfinder(0, 1);
		for (int r = 0; r < eye.rows; r++) {
			normal.at<double>(r, 0) = nfinder(re);
		}
		cout << "normal: " << normal << endl;
	} while (overlap(transformationCoords,a, eye,normal));
	

}
/**
* @brief make video
*/
void makeVideo(double scale, Mat& coordinates, vector<vector<int>>& edges, int rows, int cols, bool writeCoordinates, vector<Mat>& vidFrames, Mat &a, Mat &normal, Mat &eyepos, bool genrandom) {
	//vector<Mat> vidFrames;
	//rotate around x and y axis, scale by a factor of 250, create 360 frames (rotate by 1 degree)
	//30 frames per second
	//double scale = 250;
	double x = 0;
	double y = 0;
	double z = 0;

	//define the normal parameters
	//Mat a = (Mat_<double>(3, 1) << 0,plane,0);
	//Mat normal = (Mat_<double>(3, 1) << 0, 1, 0);
	//Mat eyepos = (Mat_<double>(3, 1) << eye, eye, eye);
	//vector<Mat> framevertices;
	int count = 0;
	ofstream coordsfile;
	ofstream coords2dfile;
	ofstream logfile;
	if (writeCoordinates) {
		coordsfile.open("coordinates.txt");
		coordsfile.precision(17);
		coords2dfile.open("coordinates2d.txt");
		coords2dfile.precision(17);
		logfile.open("log.txt");
		logfile.precision(17);
	}
	//make edges random colors
	vector<Scalar> edgeColor = vector<Scalar>();
	for (size_t q = 0; q < edges.size(); q++) {
		int r = rand() % 256;
		int g = rand() % 256;
		int b = rand() % 256;
		edgeColor.push_back(Scalar(b, g, r));
	}

	vector<Mat> transformationCoords = vector<Mat>();
	//get the coordinates matrices
	for (int i = 0; i < 360; i++) {
		Mat tranformationMat;
		findtransform(scale, x, y, z, tranformationMat);
		Mat result = tranformationMat * coordinates;
		transformationCoords.push_back(result);
		x += M_PI / 180;
		y += M_PI / 180;
	}
	//set random values for a, n, and eye
	if (genrandom) {
		setRandom(transformationCoords, a, eyepos, normal);
	}
	cout << "Plane partitions the shape and eye without intersection: " << !overlap(transformationCoords, a, eyepos, normal) << endl;
	//print the values of a, n, and eye
	//write to logfile
	if (writeCoordinates) {
		logfile << "The plane defined by (x-a)*n =0 is: " << endl;
		logfile << "\ta = " << "(" << a.at<double>(0, 0) << ", " << a.at<double>(1, 0) << ", " << a.at<double>(2, 0) << ")\n";
		logfile << "\tn = " << "(" << normal.at<double>(0, 0) << ", " << normal.at<double>(1, 0) << ", " << normal.at<double>(2, 0) << ")\n";

		//print eye
		logfile << "The eye e is: " << endl;
		logfile << "\te = " << "(" << eyepos.at<double>(0, 0) << ", " << eyepos.at<double>(1, 0) << ", " << eyepos.at<double>(2, 0) << ")\n";

		//print info about w1, w2, p0




	}
	//project and move the points onto the plane
	Mat w1, w2, p0;
	getBasis(transformationCoords[0], a, normal, eyepos, w1, w2, p0, writeCoordinates, logfile);
	//cout << p0;
	//write to file
	if (writeCoordinates) {
		stringstream w1str;
		stringstream w2str;
		stringstream p0str;
		logfile << "\nThe w1 and w2 obtained from a and b (these now are perpendicular and of magnitude 1) are:\n";
		w1str << "w1 = " << "(" << w1.at<double>(0, 0) << "," << w1.at<double>(1, 0) << "," << w1.at<double>(2, 0) << ")\n";
		w2str << "w2 = " << "(" << w2.at<double>(0, 0) << "," << w2.at<double>(1, 0) << "," << w2.at<double>(2, 0) << ")\n";
		logfile << w1str.str();
		logfile << w2str.str();
		//write center
		
		logfile << "\nThe center of the cube in first frame and its projection are:\n";
		logfile << "center = " << "(0.0,0.0,0.0)\n";
		p0str << "p0 = " << "(" << p0.at<double>(0, 0) << "," << p0.at<double>(1, 0) << "," << p0.at<double>(2, 0) << ")\n";
		logfile << p0str.str();

		logfile << "\nThe coordinates in the 2d plane x = p0 + u*w1 + v*w2 are:" << endl;
		logfile << "(where p0 is the origin, preferraby the projection of the center of the cube in first frame, w1 and w2 are 2 perpendicular vertices in the plane)\n";
		logfile << "\t" << p0str.str();
		logfile << "\t" << w1str.str();
		logfile << "\t" << w2str.str();



	}


	//vector to save the projected coordinates for first four frames
	vector<vector<vector<double>>> first4proj = vector<vector<vector<double>>>();
	vector<vector<vector<double>>> first43d = vector<vector<vector<double>>>();
	for (int i = 0; i < transformationCoords.size(); i++) {
		Mat& result = transformationCoords[i];
		//cout << "line 71\n";
		Mat projected;
		
		projectShape(result, projected, p0, a, normal,eyepos,w1,w2);
		//create video frame
		//move points to center
		if (count < 4 && writeCoordinates) {
			vector<vector<double>> c3d = vector<vector<double>>();
			vector<vector<double>> c2d = vector<vector<double>>();
			
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
				c2d.push_back({ u,v });
				c3d.push_back({ x,y,z });

			}
			first4proj.push_back(c2d);
			first43d.push_back(c3d);
			count += 1;
		}
		for (int i = 0; i < projected.cols; i++) {
			projected.at<double>(0, i) += double(cols / 2);
			projected.at<double>(1, i) = double(rows / 2) - projected.at<double>(1, i);

		}
		//draw vertices
		Mat frame = Mat(rows, cols, CV_8UC3, Scalar(255, 255, 255));;
		for (int p = 0; p < result.cols; p++) {
			int xval = cvRound(projected.at<double>(0, p));
			int yval = cvRound(projected.at<double>(1, p));

			//make vertices random colors
			//int r = rand() % 256;
			//int g = rand() % 256;
			//int b = rand() % 256;
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
			line(frame, Point(x1, y1), Point(x2, y2), Scalar(0, 0, 0), 3);


		}

		vidFrames.push_back(frame);


	}

	if (writeCoordinates) {
		for (int i = 0; i < 4; i++) {
			logfile << "\nThe frame" << i + 1 << " in 3d has the following edges:\n";
			logfile << "(should be 12 edges, 12 pairs of 3d coordinates (x,y,z))\n";
			for (vector<int>& edge : edges) {
				int e1 = edge[0];
				int e2 = edge[1];
				logfile << "\t(" << first43d[i][e1][0] << "," << first43d[i][e1][1] << ","<< first43d[i][e1][2] << "), ";
				logfile << "(" << first43d[i][e2][0] << "," << first43d[i][e2][1] << "," << first43d[i][e2][2] << ")\n";
			}

			logfile << "\nThe frame" << i + 1 << " in 2d has the following edges :\n";
			logfile << "(should be 12 edges, 12 pairs of 2d coordinates (u,v), if the vertices do not cross the plane these are the projections of the 3d points on the plane)\n";
			for (vector<int>& edge : edges) {
				int e1 = edge[0];
				int e2 = edge[1];
				logfile << "\t(" << first4proj[i][e1][0] << "," << first4proj[i][e1][1]<<"), ";
				logfile << "(" << first4proj[i][e2][0] << "," << first4proj[i][e2][1] << ")\n";
			}
		}


	}
	//close file
	if (writeCoordinates) {
		coordsfile.close();
		coords2dfile.close();
		logfile.close();
	}




}

//void part2() {
//	//calculations
//	//Mat coordinates;
//	//vector<vector<int>> edges;
//	//if (shape == 0) {
//	Mat cubecoordinates = (Mat_<double>(4, 8) << 1, 1, 1, 1, -1, -1, -1, -1,
//		1, 1, -1, -1, 1, 1, -1, -1,
//		1, -1, 1, -1, 1, -1, 1, -1,
//		1, 1, 1, 1, 1, 1, 1, 1);
//	vector<vector<int>> cubeedges = { {0,1},{1,5},{1,3},{0,2},{0,4},{4,5},{4,6},{7,5},{6,7},{6,2},{3,2},{3,7} };
//
//	//}
//	//else {
//	double sqrt3 = sqrt(3);
//	Mat tetracoordinates = (Mat_<double>(4, 4) << 1 / sqrt3, 0, -sqrt3 / 6, -sqrt3 / 6,
//		0, 0, 0.5, -0.5,
//		0, 2 / sqrt(6), 0, 0,
//		1, 1, 1, 1);
//	vector<vector<int>> tetraedges = { {0,1},{0,3},{0,2},{1,3},{1,2},{2,3} };
//	//}
//	//making the frames
//	vector<Mat> vidFrames;
//	//makeVideo("tetrahedron.avi", coordinates, edges, 600, 800);
//	makeVideo(150, cubecoordinates, cubeedges, 600, 800, true, vidFrames, 300, 650, false);
//	makeVideo(250, tetracoordinates, tetraedges, 600, 800, false, vidFrames, 300, 650, false);
//	//write to file
//	VideoWriter outputvid;
//	int codec = VideoWriter::fourcc('M', 'J', 'P', 'G');  // select desired codec (must be available at runtime)
//	double fps = 12;                          // framerate of the created video stream
//	string filename = "rotation.avi";             // name of the output video file
//	outputvid.open(filename, codec, fps, vidFrames[0].size(), true);
//	for (int c = 0; c < vidFrames.size(); c++) {
//		outputvid << vidFrames[c];
//	}
//	outputvid.release();
//}

void convertToVector(string s, Mat& matrix) {
	size_t pos = 0;
	std::string token;
	string delimiter = ",";
	int idx = 0;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		//std::cout << token << std::endl;
		s.erase(0, pos + delimiter.length());
		matrix.at<double>(idx, 0) = stod(token);
		idx += 1;
	}
	matrix.at<double>(idx, 0) = stod(s);
	//std::cout << s << std::endl;
}
void part3(int argc, char* argv[]) {

	//read command line arguments
	Mat a = (Mat_<double>(3, 1) << 10, 301, 23);
	Mat normal = (Mat_<double>(3, 1) << 1, 9, 1);
	Mat eyepos = (Mat_<double>(3, 1) << 50.5, 651, 349);

	bool generaterandom = false;
	for (int i = 1; i < argc; i++) {
		string flag = string(argv[i]);
		if (flag == "-a") {
			string stra = string(argv[++i]);
			//cout << stra;
			convertToVector(stra.substr(1, stra.size() - 2), a);
			generaterandom = false;
		}
		else if (flag == "-n") {
			string strn = string(argv[++i]);
			convertToVector(strn.substr(1, strn.size() - 2), normal);
			generaterandom = false;
		}
		else if (flag == "-e") {
			string streye = string(argv[++i]);
			convertToVector(streye.substr(1, streye.size() - 2), eyepos);
			generaterandom = false;
		}
	}
	//cout << "a" << a << endl;
	//cout << "normal" << normal << endl;
	//cout << "eye" << eyepos << endl;

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
	makeVideo(150, cubecoordinates, cubeedges, 600, 800, true, vidFrames, a, normal, eyepos, generaterandom);
	makeVideo(250, tetracoordinates, tetraedges, 600, 800, false, vidFrames, a, normal, eyepos, generaterandom);
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
int main(int argc, char *argv[]) {
	//part2();
	part3(argc, argv);

}
