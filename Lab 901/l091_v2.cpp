
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/video.hpp>

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
		line(frame, imagePoints[edge1], imagePoints[edge2], Scalar(255, 128, 128), 3);
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
			/*vector<Point3f> selectedObj;
			vector<Point2f> selectedCorners;
			for (int i = 0; i < corner2d.size(); i++) {
				if (i % 2 == 0) {
					selectedObj.push_back(objectPoints[i]);
					selectedCorners.push_back(corner2d[i]);
				}*/
			//}
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
	//undistort the first points in objectPointsList
	// find the different between them
	// 
	//values for LK Optical flow tracking
	Mat oldFrame = vidFrames[0];
	Mat oldGray;
	cvtColor(oldFrame, oldGray, COLOR_BGR2GRAY);
	vector<Point2f> trackedInFrame;
	vector<Point2f> prevfeaturePoints = planePointsList[0];
	vector<Point3f> prevObjPoints = objectPoints;
	//use feature points from LK optical tracking for determination of points to use in soldpnp
	for (int i = 1; i < vidFrames.size(); i++) {

		//convert to grayscale
		Mat gray;
		cvtColor(vidFrames[i], gray, COLOR_BGR2GRAY);
		//attempt to find chessboard pattern
		vector<Point2f> corner2d;
		//bool found = findChessboardCorners(gray, Size(boardCols, boardRows), corner2d, CALIB_CB_ADAPTIVE_THRESH + CALIB_CB_NORMALIZE_IMAGE + CALIB_CB_FAST_CHECK);
		//if chessboard pattern detected: https://docs.opencv.org/4.x/d5/d1f/calib3d_solvePnP.html
		//extrensic camera properties that change frame to frame
		if (i%20==0) { //reasses the corners every 20 frames
			//use the corners detected from findChessBoard to calibrate lkoptical flow
			bool found = findChessboardCorners(gray, Size(boardCols, boardRows), corner2d, CALIB_CB_ADAPTIVE_THRESH + CALIB_CB_NORMALIZE_IMAGE + CALIB_CB_FAST_CHECK);
			if (found) {
				TermCriteria criteria = TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 10, 0.01);
				cornerSubPix(gray, corner2d, Size(11, 11), Size(-1, -1), criteria);
				solvePnP(objectPoints, corner2d, cameraMatrix, dist, rvecs, tvecs);
				//project onto image matrix
				vector<Point2f> imagePoints;
				projectPoints(coordinates, rvecs, tvecs, cameraMatrix, dist, imagePoints);
				drawFigure(vidFrames[i], imagePoints, edges);
				oldFrame = gray;
				prevObjPoints = objectPoints; //updates prevObjPoints with valid points from this iteration
				prevfeaturePoints = corner2d; //update prevfeaturePoints with good 2dcorners from this iteration
				continue;
			}
			else {
				corner2d.clear();
			}
			


		}
		
			TermCriteria criteria = TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 10, 0.01);
			//cornerSubPix(gray, corner2d, Size(11, 11), Size(-1, -1), criteria);
			vector<uchar> status;
			vector<float> err;
			vector<Point2f> cornerCandidates;
			calcOpticalFlowPyrLK(oldGray, gray, prevfeaturePoints, cornerCandidates, status, err, Size(15, 15), 2, criteria);
			//object Points and good corner points
			vector<Point3f> goodObjects;
			int count = 0;
			for (uint i = 0; i < prevfeaturePoints.size(); i++) {
				if (status[i] == 1){
					goodObjects.push_back(prevObjPoints[i]); // tracks the good object points
					corner2d.push_back(cornerCandidates[i]); //tracks the valid 2d corners
					count++;
				}
			}
			cout << "frame " << i << "length of good objects and corner2d: " << corner2d.size() << " " << goodObjects.size() << endl;
			//use solve pnp to find the projection matrix

			solvePnP(goodObjects, corner2d, cameraMatrix, dist, rvecs, tvecs);
			//project onto image matrix
			vector<Point2f> imagePoints;
			projectPoints(coordinates, rvecs, tvecs, cameraMatrix, dist, imagePoints);
			drawFigure(vidFrames[i], imagePoints, edges);
			prevObjPoints = goodObjects; //updates prevObjPoints with valid points from this iteration
			prevfeaturePoints = corner2d; //update prevfeaturePoints with good 2dcorners from this iteration
			oldGray = gray;
			//cameracount += 1;
		
		
		//cameracount += 1;


		/*}
		else {

		}*/
		if (CV_VERSION_MAJOR < 4 || CV_VERSION_MAJOR==4 && CV_VERSION_MINOR < 5) {
			flip(vidFrames[i], vidFrames[i], -1);
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
		2, 0, 2, 0, 2, 0, 2, 0);
	//translate the cubcoordinates upward
	//for (int i = 0; i < 8; i++) {
	//	cubecoords.at<float>(2, i) += 1;
	//	//cubecoords.at<float>(2, i) *= 4;
	//}
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
