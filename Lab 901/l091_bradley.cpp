#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

void augmentedCube()
{
    Mat gray;
    Mat guessframe;
    Mat frame;

    VideoCapture vidIn("withChessBoard.MOV");
    int fps = vidIn.get(CAP_PROP_FPS);
    Size boardsize(7, 7);
    Size framesize;
    vector<Point2f> corners;
    vector<Point3f> vertices;
    vertices.push_back(Point3f{ -1,1,0 });
    vertices.push_back(Point3f{ 1,1,0 });
    vertices.push_back(Point3f{ -1,-1,0 });
    vertices.push_back(Point3f{ 1,-1,0 });
    vector<vector<Point2f>> previouscorners;
    vector<vector<Point2f>> allcorners;
    vector<vector<Point3f>> previous3dcoords;
    vector<uchar> status;
    vector<float> error;
    vidIn >> frame;
    VideoWriter video("vr.avi", VideoWriter::fourcc('M', 'J', 'P', 'G'), fps, frame.size());

    int framecache = 0;
    while (!frame.empty())
    {
        framesize = frame.size();
        flip(frame, frame, -1);
        cvtColor(frame, gray, COLOR_BGR2GRAY);
        vector<Point2f> framecorners;
        vector<Point2f> framecubecorners;
        if (findChessboardCorners(frame, boardsize, framecorners))
        {
            framecubecorners.push_back(framecorners[16]);
            framecubecorners.push_back(framecorners[18]);
            framecubecorners.push_back(framecorners[30]);
            framecubecorners.push_back(framecorners[32]);
        }
        else
        {
            TermCriteria criteria = TermCriteria((TermCriteria::COUNT)+(TermCriteria::EPS), 10, 0.01);
            calcOpticalFlowPyrLK(guessframe, gray, corners, framecubecorners, status, error, Size(15, 15), 2, criteria);
        }
        gray.copyTo(guessframe);
        allcorners.push_back(framecubecorners);
        corners = framecubecorners;
        if (framecache == 20)
        {
            previouscorners.push_back(framecubecorners);
            previous3dcoords.push_back(vertices);
            framecache = 0;
        }
        framecache += 1;
        vidIn >> frame;
    }

    Mat cameraMatrix;
    Mat distCoeffs;
    Mat rvecs;
    Mat tvecs;
    Mat stdDeviationsIntrinsics;
    Mat stdDeviationsExtrinsics;
    Mat perViewErrors;
    calibrateCamera(previous3dcoords, previouscorners, framesize, cameraMatrix, distCoeffs, rvecs, tvecs, stdDeviationsIntrinsics, stdDeviationsExtrinsics, perViewErrors);

    vidIn = VideoCapture("withChessBoard.MOV");
    vidIn >> frame;
    int index = 0;
    vector<Point3f> realvertices;
    realvertices.push_back(Point3f{ -1,1,0 });
    realvertices.push_back(Point3f{ 1,1,0 });
    realvertices.push_back(Point3f{ -1,-1,0 });
    realvertices.push_back(Point3f{ 1,-1,0 });
    realvertices.push_back(Point3f{ -1,1,2 });
    realvertices.push_back(Point3f{ 1,1,2 });
    realvertices.push_back(Point3f{ -1,-1,2 });
    realvertices.push_back(Point3f{ 1,-1,2 });
    while (!frame.empty())
    {
        flip(frame, frame, -1);
        solvePnP(vertices, allcorners[index], cameraMatrix, distCoeffs, rvecs, tvecs);
        vector<Point2f> projectedPoints;
        projectPoints(realvertices, rvecs, tvecs, cameraMatrix, distCoeffs, projectedPoints);
        line(frame, projectedPoints[0], projectedPoints[1], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[0], projectedPoints[2], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[1], projectedPoints[3], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[2], projectedPoints[3], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[4], projectedPoints[5], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[4], projectedPoints[6], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[5], projectedPoints[7], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[6], projectedPoints[7], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[0], projectedPoints[4], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[1], projectedPoints[5], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[2], projectedPoints[6], Scalar(147, 20, 255), 5);
        line(frame, projectedPoints[3], projectedPoints[7], Scalar(147, 20, 255), 5);
        video.write(frame);
        index += 1;
        vidIn >> frame;
    }
    video.release();
    vidIn.release();
}

int main(int argc, char** argv)
{
    augmentedCube();
}