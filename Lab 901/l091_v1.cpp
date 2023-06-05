
#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

using namespace std;
using namespace cv;
/**
* @brief Reads video from a specified filename and saves the frames to an vector
* @param filename name of the video file to read
* @param vidFrames reference to vector to save frames
*/
//void readVideo(string filename, vector<Mat> &vidFrames) {
//	//https://stackoverflow.com/questions/13709274/reading-video-from-file-opencv
//	VideoCapture capture(filename);
//	Mat frame;
//	if (!capture.isOpened()) {
//		throw "Error when reading " + filename;
//	}
//	while (true) {
//		capture >> frame;
//		if (frame.empty()) {
//			break;
//		}
//		vidFrames.push_back(frame);
//	}
//	
//}
///**
//* @brief Saves vector of frames as video
//* Precondition - filename must end in .avi
//* @param filename name of the file to save video to
//*/
//void writeVideo(string filename) {
//	//tutorial: https://docs.opencv.org/2.4.1/doc/tutorials/highgui/video-write/video-write.html
//	//VideoWriter outputVideo;
//	//outputVideo.open(filename, ex=-1, )
//}

/**
* @brief draw image on a frame
* @param frame the frame to draw the image on
* @param imagePoints the vector of points
* @param edges list of edges to draw
*/
void drawFigure(Mat& frame, vector<Point2f> imagePoints, vector<vector<int>> edges) {
	for (int i = 0; i < edges.size(); i++) {
		int edge1 = edges[i][0];
		int edge2 = edges[i][1];
		line(frame, imagePoints[edge1], imagePoints[edge2], Scalar(255, 255, 255), 3);
	}
}
/**
* @brief place the cube on every frame. coordinates should be prescaled.
*/
void placeObject(vector<Mat>& vidFrames, vector<Mat>& outputFrames, Mat& coordinates, vector<vector<int>> edges) {
	//identify good features to track
	// 
	////3. for every frame:
	int boardRows = 7;
	int boardCols = 7;
	vector<Point3f> objectPoints = vector<Point3f>();
	//this way center square is assigned origin by starting at -3
	for (int i = -3; i < -3 + boardRows; i++) {
		for (int j = -3; j < -3 + boardCols; j++) {
			Point3d pt = Point3d(i, j, 0);
			objectPoints.push_back(pt);
		}
	}
	int cameracount = 0;
	//intrensic camera properties that do not change frame to frame
	Mat cameraMatrix;
	Mat dist;
	Mat rvecs;
	Mat tvecs;
	//calibrate camera first using first 30 frames
	//list of images used for calibration
	vector<vector<Point3f>> objPointsList;
	vector<vector<Point2f>> planePointsList;
	for (int c = 0; c < 14; c++) {


		//convert to grayscale
		Mat gray;
		cvtColor(vidFrames[c], gray, COLOR_BGR2GRAY);
		//attempt to find chessboard pattern
		vector<Point2f> corner2d;
		bool found = findChessboardCorners(gray, Size(boardCols, boardRows), corner2d, CALIB_CB_ADAPTIVE_THRESH + CALIB_CB_NORMALIZE_IMAGE + CALIB_CB_FAST_CHECK);
		//if chessboard pattern detected: https://docs.opencv.org/4.x/d5/d1f/calib3d_solvePnP.html
		
		//extrensic camera properties that change frame to frame
		Mat dist;
		Mat rvecs;
		Mat tvecs;
		cout << "callibration: " << c << ": " << found << endl;
		if (found) {

			TermCriteria criteria = TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 30, 0.001);
			cornerSubPix(gray, corner2d, Size(11, 11), Size(-1, -1), criteria);

			//use append values to calibratecamera's coordinates
			objPointsList.push_back(objectPoints);
			planePointsList.push_back(corner2d);
			//vector<Point2f> imagePoints;
		}
		//Mat r, t, d;
		//calibrateCamera(objPointsList, planePointsList, vidFrames[0].size(), cameraMatrix, dist, rvecs, tvecs);

	}
	//Mat r, t, d;
	//calibrate the camera
	calibrateCamera(objPointsList, planePointsList, vidFrames[0].size(), cameraMatrix, dist, rvecs, tvecs);
	//values for LK Optical flow tracking
	Mat oldFrame, oldGray;
	vector<Point2f> trackedInFrame;
	vector<Point2f>& featurePoints = planePointsList[0];
	for (int i = 0; i < 30; i++) {

		//convert to grayscale
		Mat gray;
		cvtColor(vidFrames[i], gray, COLOR_BGR2GRAY);
		//attempt to find chessboard pattern
		vector<Point2f> corner2d;
		bool found = findChessboardCorners(gray, Size(boardCols, boardRows), corner2d, CALIB_CB_ADAPTIVE_THRESH + CALIB_CB_NORMALIZE_IMAGE + CALIB_CB_FAST_CHECK);
		//if chessboard pattern detected: https://docs.opencv.org/4.x/d5/d1f/calib3d_solvePnP.html
		//extrensic camera properties that change frame to frame

		cout << i << ": " << found << endl;
		if (found) {
			//apply undistortion transformation and save a copy of the frame (calculations contained in solve PnP
			//4. determine the coordinates of the center four squares, let that be the origin, and the table be z=0
			//5. find the homography matrix of each frame
			//6. transform verticies of cube using homography matrix
			TermCriteria criteria = TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 30, 0.001);
			cornerSubPix(gray, corner2d, Size(11, 11), Size(-1, -1), criteria);

			cout << "frame " << i << endl;
			//use solve pnp to find the projection matrix

			solvePnP(objectPoints, corner2d, cameraMatrix, dist, rvecs, tvecs);
			vector<Point2f> imagePoints;
			projectPoints(coordinates, rvecs, tvecs, cameraMatrix, dist, imagePoints);
			drawFigure(vidFrames[i], imagePoints, edges);
			cameracount += 1;


		}
		else {

		}
		outputFrames.push_back(vidFrames[i]);
		//else if pattern not found:
			//use lukas-kanade optical flow tracking algorithm to identify chessboard points
			//find homography matrix and transform points
			//https://docs.opencv.org/3.4/d4/dee/tutorial_optical_flow.html
	}


}
void part1() {
	//Strategy
	//1. Read video and extract all the frames, storing them in an vector
	string filename = "withChessBoard.MOV";
	vector<Mat> vidFrames = vector<Mat>();
	VideoCapture capture(filename);

	if (!capture.isOpened()) {
		throw "Error when reading " + filename;
	}
	while (true) {
		Mat frame;
		capture >> frame;
		if (frame.empty()) {
			break;
		}
		vidFrames.push_back(frame);
	}
	cout << "length of vidFrames: " << vidFrames.size();
	//coordinates of cube and edgelist
	int coordsize[] = { 1,8 };
	cout << "line 152" << endl;
	//Mat cubecoords = (Mat_<float>(3,8) << 1, 1, 1, 1, -1, -1, -1, -1,
	//									  1, 1, -1, -1, 1, 1, -1, -1,
	//									  1, -1, 1, -1, 1, -1, 1, -1);
	Mat cubecoords = (Mat_<float>(3, 8) << 1, 1, 1, 1, -1, -1, -1, -1,
		1, 1, -1, -1, 1, 1, -1, -1,
		1, -1, 1, -1, 1, -1, 1, -1);
	//translate the cubcoordinates upward
	for (int i = 0; i < 8; i++) {
		cubecoords.at<float>(2, i) += 1;
		//cubecoords.at<float>(2, i) *= 4;
	}
	vector<vector<int>> edges = { {0,1},{1,5},{1,3},{0,2},{0,4},{4,5},{4,6},{7,5},{6,7},{6,2},{3,2},{3,7} };
	vector<Mat> outputFrames;
	placeObject(vidFrames, outputFrames, cubecoords, edges);
	//2. using the first 15-20 frames, identify camera calibration (calculation contained in solve PnP)
	//3. for every frame:
		//if chessboard pattern detected: https://docs.opencv.org/4.x/d5/d1f/calib3d_solvePnP.html
			//apply undistortion transformation and save a copy of the frame (calculations contained in solve PnP
			//4. determine the coordinates of the center four squares, let that be the origin, and the table be z=0
			//5. find the homography matrix of each frame
			//6. transform verticies of cube using homography matrix
			//save the frame and previous vertices
		//else if pattern not found:
			//use lukas-kanade optical flow tracking algorithm to identify chessboard points
			//find homography matrix and transform points
			//https://docs.opencv.org/3.4/d4/dee/tutorial_optical_flow.html


	//save frames to a file
	VideoWriter outputvid;
	//string outname = "vr.avi";
	//Size S = Size((int)capture.get(CAP_PROP_FRAME_WIDTH),    //Acquire input size
	//	(int)capture.get(CAP_PROP_FRAME_HEIGHT));
	//outputvid.open(outname, -1, capture.get(CAP_PROP_FPS), S, true);
	int codec = VideoWriter::fourcc('M', 'J', 'P', 'G');  // select desired codec (must be available at runtime)
	double fps = capture.get(CAP_PROP_FPS);                          // framerate of the created video stream
	string outfile = "vr.avi";             // name of the output video file
	outputvid.open(outfile, codec, fps, outputFrames[0].size(), true);
	for (int i = 0; i < outputFrames.size(); i++) {
		outputvid << outputFrames[i];
	}
	outputvid.release();
}
int main() {
	part1();
}
