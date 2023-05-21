
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <sstream>
#include <float.h>
#include <list>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <iterator>
#include <stack>
#include <unordered_set>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;

 // #include <numbers>

#define _USE_MATH_DEFINES

using namespace std::chrono;
using namespace std;

const int PRECISION = 23;
const double ERROR = 1e-7;
int CIRCLE_RADIUS = 180;
// typedef std::numeric_limits< double > dbl;
// ofstream results;

typedef std::numeric_limits<double> dbl;


void cannyEdge(Mat &image, Mat &output, int low, int high) {
    //applies the canny edge detector to the grayscale coins image
    //returns the output array
    Mat blurredimage;
    //erode(image, image, getStructuringElement(MORPH_ELLIPSE, Size(4,4)),Point(-1,-1),4);
    //dilate(image, image, getStructuringElement(MORPH_ELLIPSE, Size(4, 4)), Point(-1,-1),4);
    //GaussianBlur(image, image, Size(3, 3), 1);
    //stackBlur(image, image,3);
    
    //bilateralFilter(image, blurredimage, 5, 20, 20);
    dilate(image, image, getStructuringElement(MORPH_ELLIPSE, Size(3, 3)),Point(-1,-1),3);
    erode(image, image, getStructuringElement(MORPH_CROSS, Size(4, 4)),Point(-1,-1),3);
    image.convertTo(image, -1, 1.2);
    //imwrite("blurred.jpg", blurredimage);
    Canny(image, output, low, high);
}
void findCircles(Mat& grayImage, vector<Vec3f> &circles, double highThreshold, double minRad, double maxRad, double thresholdCount, double minDist) {
    /*Mat blurredimage;
    bilateralFilter(grayImage, blurredimage, 5, 20, 20);
    blur(blurredimage, blurredimage, Size(3, 3));*/
    HoughCircles(grayImage,circles,HOUGH_GRADIENT,1,minDist,highThreshold,thresholdCount,minRad,maxRad);
}


void identifyCoins(Mat &colorImg, bool &ishard,vector<Vec3f> &circleCandidates, vector<int> &identified){
    //isolate the pennies in the image by color
    //obtain the average radii of the penny, and then scale coins accordingly
    //identify reference radii
    //for c in circles:
        //classify each circle
    //https://stackoverflow.com/questions/38673027/how-to-identify-gold-color-in-image
    //blur(img, img, Size(50, 50) );
    //find the average penny radius to use as a reference size
    Mat& img = colorImg;
    GaussianBlur(img, img, Size(7, 7), 1.4);
    Mat adjusted;
    img.convertTo(adjusted, -1, 1.5);
    cvtColor(adjusted, adjusted, COLOR_BGR2HSV);
    Mat saturation;
    extractChannel(adjusted, saturation, 1);
    Mat binarized;
    threshold(saturation, binarized, 50, 255, THRESH_BINARY);
    //dilate(binarized, binarized, getStructuringElement(MORPH_DILATE, Size(4, 4)));
    morphologyEx(binarized, binarized, MORPH_CLOSE, getStructuringElement(MORPH_ELLIPSE, Size(4, 4)));
    //floodFill(binarized, cv::Point(0, 0), Scalar(255));
    //isolate the pennies,
    Mat isolated;
    bitwise_and(saturation, saturation, isolated, binarized);
    imwrite("binary.jpg", isolated);
    vector<Vec3f> circles;
    HoughCircles(isolated, circles, HOUGH_GRADIENT, 1, 80, 100, 25, 75, 90);
    //compute the average radius of pennies
    double avgPenny = 0;
    for (size_t i = 0; i < circles.size(); i++) {
        avgPenny += circles[i][2];
    }
    avgPenny /= circles.size();
    cout << "Average penny size: " << avgPenny;
    if (abs(avgPenny) < 1) {
        avgPenny = 87;
    }
    double tolerance = 0.05;
    if (!ishard) {
        double quarterratio = 1.25;
        double nickelratio = 1.13;
        double pennyratio = 1.03;
        //compute ratio
        /*
       * Representation of coins:
       0 -- penny
       1 - nickel
       2 - dime
       3- quarter
       4 - half-dollar
       5 - dollar
       */
        //if computed ratio falls within 0.05 of reference ratio, classify as said coin

        for (size_t idx = 0; idx < circleCandidates.size(); idx++) {
            double radius = circleCandidates[idx][2];
            double ratio = radius / avgPenny;


            /* Get white area (aka high saturation area) */
            //get mask of region containing circle
            //count number of white pixels
            //if number of white pixels exceeds 25%, count as nickel
            int x = cvRound(circleCandidates[idx][0]);
            int y = cvRound(circleCandidates[idx][1]);
            int length = cvRound(radius * 2);
            int halflen = cvRound(radius);
            const double pi = 3.1415926535;
            int area = pi * radius * radius;
            double totalArea = 0;
            if (x < binarized.cols - length / 2 && y < binarized.rows - length / 2 && x > halflen && y > halflen) {
                Rect dim(x - halflen, y - halflen, length, length);
                Mat roi = binarized(dim);
                std::vector<std::vector<cv::Point> > contours;
                std::vector<cv::Vec4i> hierarchy;
                //floodFill(roi, Point(0, 0), 255);
                //invert(roi, roi);
                //bitwise_not(roi, roi);

                findContours(roi, contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
                
                
                bool penny = false;

                //Mat colorImg;
                //cvtColor(roi, colorImg, COLOR_GRAY2BGR);
                //drawContours(colorImg, contours, -1, Scalar(0,0,255),3);
                //imshow("owo", colorImg);
                
                for (int i = 0; i < contours.size(); i++)
                {
                    /*if (contourArea(contours[i]) > thresholdarea) {
                        identified.push_back(0);
                        penny = true;
                        break;
                    }*/
                    totalArea += contourArea(contours[i]);
                }
                



            }
            
            //classify coins
            if (ratio < pennyratio) {
                identified.push_back(0);
            }
            else if (abs(ratio - pennyratio) < tolerance) {
                if (totalArea < 0.25*area) {
                    identified.push_back(1);
                }
                else {
                    identified.push_back(0);
                }
                
            }else if(abs(ratio-nickelratio) < tolerance){
              
                if (totalArea > 0.25*area) {
                    identified.push_back(0);
                }
                else {
                    identified.push_back(1);
                }
            }
            else if (abs(ratio - quarterratio) < tolerance) {
                if (totalArea > 0.5*area) {
                    identified.push_back(0);
                }
                else {
                    identified.push_back(3);
                }
                
            }
            else {
                identified.push_back(3);
            }
        }
    }
}

void summary(vector<int>& coins, string filename)
{
    ofstream results;
    results.open(filename);
    vector<int> counts = vector<int>(6, 0);
    for (const int& coinType : coins)
    {
        // switch (coinType)
        // {
        // case(0):
        //     penny++;
        //     break;
        // case(1):
        //     dime++;
        //     break;
        // case(2):
        //     nickel++;
        //     break;
        // case(3):
        //     quarter++;

        // default:
        //     break;
        // }
        if (coinType >= 0)
        {
            counts[coinType] += 1;
        }
    }

    results << "penny - " << counts[0] << endl;
    cout << "\npenny - " << counts[0] << endl;

    results << "nickel - " << counts[1] << endl;
    cout << "nickel - " << counts[1] << endl;

    results << "dime - " << counts[2] << endl;
    cout << "dime - " << counts[2] << endl;

    results << "quarter - " << counts[3] << endl;
    cout << "quarter - " << counts[3] << endl;

    results << "half dollar - " << counts[4] << endl;
    cout << "half dollar - " << counts[4] << endl;

    results << "dollar - " << counts[5] << endl;
    cout << "dollar - " << counts[5] << endl;
    // results.precision(2);
    cout << std::fixed << std::showpoint;
    cout << std::setprecision(2);
    results << std::fixed << std::showpoint;
    results << std::setprecision(2);
    int totalAmt = counts[0] + 5 * counts[1] + 10 * counts[2] + 25 * counts[3] + 50 * counts[4] + 100 * counts[5];
    double totalAmtDollar = totalAmt / 100.0;

    results << "Total: $" << totalAmtDollar << endl;
    cout << "Total: $" << totalAmtDollar << endl;
    results.close();
}
/**
 * @brief Part 1 of Lab 7
 *
 * @param argc
 * @param argv
 */
void part3(int argc, char** argv)
{
    /*
    Optimal Parameters for Hard image:
    Q1:
    lt = 100, ht = 250
    TC =
    https://stackoverflow.com/questions/26855264/identifying-different-coin-values-from-an-image-using-matlab
    */
    string imageName = "coinsHarder.jpg";
    
    int high = 200;
    int low = high/2;
    string fg = "imageg.jpg";
    string ff = "imagef.jpg";
    string fci = "imageCircles.jpg";
    string fco = "imageCoins.jpg";
    string fr = "results.txt";
    int thresholdCount = 25;
    int tCircle = 20;
    bool override = false;
    bool tCircleOverride = true;
    //bool thresholdOverride = false;
    bool hardimage = false;
    int minCenterDist = 50;
    int minRad = 80;
    int maxRad = 120;


    // read the command line arguments
    for (int i = 1; i < argc; i++)
    {
        string flag = string(argv[i]);
        if (flag == "-f")
        {
            // cout << "flag read" << endl;
            imageName = argv[++i];
        }
        //else if (flag == "-lt")
        //{
        //    stringstream ss(argv[++i]);
        //    ss >> low;
        //    //thresholdOverride = true;
        //}
        else if (flag == "-ht")
        {
            stringstream ss(argv[++i]);
            ss >> high;
            low = high / 2;
            //thresholdOverride = true;
        }
        else if (flag == "-ff")
        {
            ff = argv[++i];
        }
        else if (flag == "-fg")
        {
            fg = argv[++i];
        }
        else if (flag == "-TC")
        {
            stringstream ss(argv[i + 1]);
            ss >> thresholdCount;
            override = true;
        }
        else if (flag == "-TCircle")
        {
            stringstream ss(argv[++i]);
            ss >> tCircle;
            tCircleOverride = true;
        }
        else if (flag == "-fCi")
        {
            stringstream ss(argv[++i]);
            ss >> fci;
        }
        else if (flag == "-fCo")
        {
            stringstream ss(argv[++i]);
            ss >> fco;
        }
        else if (flag == "-fR")
        {
            stringstream ss(argv[++i]);
            ss >> fr;
        }
        else if (flag == "-mincentdist")
        {
            stringstream ss(argv[++i]);
            ss >> minCenterDist;
        }
        else if (flag == "-minR")
        {
            stringstream ss(argv[++i]);
            ss >> minRad;
        }
        else if (flag == "-maxR")
        {
            stringstream ss(argv[++i]);
            ss >> maxRad;
        }
        else if (flag == "-hard")
        {   
            hardimage = true;
        }

        // i++;
    }
    Mat image;
    image = imread(imageName, IMREAD_COLOR);
    //create grayscale image
   /* Mat highContrast;
    image.convertTo(highContrast, -1, 1.5);*/
    Mat imageg;
    cvtColor(image, imageg, COLOR_BGR2GRAY);
    ////apply glare elimination
    //cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
    //clahe->setClipLimit(4);
    //clahe->setTilesGridSize(Size(2,2));
    ////cv::Mat dst;
    //clahe->apply(imageg, imageg);
    
    Mat edges;
    //apply canny edge
    cannyEdge(imageg, edges, low, high);
    imwrite(fg, imageg);
    imwrite(ff, edges);
    //identify circles
    vector<Vec3f> circles;
    findCircles(imageg,circles,high,minRad,maxRad,thresholdCount,minCenterDist);
    //classify coins
    vector<int> coinTypes;
    identifyCoins(image,hardimage,circles,coinTypes);
    //draw circles on the original image
    cout << "Len of Circles: " << circles.size() << endl;
    for (size_t i = 0; i < circles.size(); i++)
    {
        Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
        int radius = cvRound(circles[i][2]);
        // draw the circle center
        circle(image, center, 3, Scalar(0,0, 255), -1, 8, 0);
        // draw the circle outline
        circle(image, center, radius, Scalar(0, 0, 255), 3, 8, 0);
    }
    imwrite(fci, image);
    
    cout << "Len of identified coins:  " << coinTypes.size() << endl;
    for (size_t i = 0; i < coinTypes.size(); i++) {
        Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
        int radius = cvRound(circles[i][2]);
        if (coinTypes[i] == 0) {
            circle(image, center, radius, Scalar(0,0,255),3,8);
        }
        else if (coinTypes[i] == 1) {
            circle(image, center, radius, Scalar(0, 255, 255), 3, 8);
        }
        else if (coinTypes[i] == 2) {
            circle(image, center, radius, Scalar(255, 0, 0), 3, 8);
        }
        else if (coinTypes[i] == 3) {
            circle(image, center, radius, Scalar(128, 0, 128), 3, 8);
        }
        else if (coinTypes[i] == 4) {
            circle(image, center, radius, Scalar(0, 255, 0), 3, 8);
        }
        else if (coinTypes[i] == 5) {
            circle(image, center, radius, Scalar(128, 0, 255), 3, 8);
        }
        else {
            circle(image, center, radius, Scalar(255, 255, 255), 3, 8);
        }
        //cout << coinTypes[i] << " " << endl;
    }
    imwrite(fco, image);
    //summary();
    
    
}
int main(int argc, char** argv)
{
    // part1();
    auto start = high_resolution_clock::now();
    // part1(argc, argv);
    // part2(argc, argv);
    part3(argc, argv);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    std::cout << "Time taken " << duration.count() << endl;
    waitKey(0);

    ////looking at hsv image
    ////Mat img;
    ////img = imread("coinsHardest.jpg");
    ////Mat imginverted;
    ////imginverted = ~img;
    ////Mat hsvimg;
    ////cvtColor(imginverted, hsvimg, COLOR_BGR2HSV);
    ////imwrite("coinsHardestinverted.jpg", imginverted);
    ////imwrite("invertedhsv.jpg", hsvimg);
    ////waitKey(0);
    ////destroyAllWindows();

   /* testing the enhance for coin detection*/
    //Mat img;
    //img = imread("coinsEasy.jpg");
    ////blur(img, img, Size(50, 50) );
    //GaussianBlur(img, img, Size(7, 7), 1.4);
    //Mat adjusted;
    //img.convertTo(adjusted, -1, 1.5);
    //cvtColor(adjusted, adjusted, COLOR_BGR2HSV);
    //Mat saturation;
    //extractChannel(adjusted, saturation, 1);
    //Mat binarized;
    //threshold(saturation, binarized, 80,255, THRESH_BINARY);
    ////dilate(binarized, binarized, getStructuringElement(MORPH_DILATE, Size(4, 4)));
    //morphologyEx(binarized, binarized, MORPH_CLOSE, getStructuringElement(MORPH_ELLIPSE,Size(4,4)));
    ////floodFill(binarized, cv::Point(0, 0), Scalar(255));
    ////isolate the pennies,
    //Mat isolated;
    //bitwise_and(saturation, saturation, isolated, binarized);
    //vector<Vec3f> circles;
    //HoughCircles(isolated, circles, HOUGH_GRADIENT,1, 80, 100, 25, 75, 97);
    //cout << circles.size() << endl;
    //int avgPenny = 0;
    //cvtColor(isolated, isolated, COLOR_GRAY2RGB);
    //for (size_t i = 0; i < circles.size(); i++)
    //{
    //    Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
    //    int radius = cvRound(circles[i][2]);
    //    // draw the circle center
    //    circle(img, center, 3, Scalar(0, 255, 0), -1, 8, 0);
    //    // draw the circle outline
    //    circle(img, center, radius, Scalar(0, 0, 255), 3, 8, 0);
    //    avgPenny += radius;
    //}
    //avgPenny /= circles.size();
    //cout << avgPenny << endl;
    //
    //namedWindow("Display Image", WINDOW_NORMAL);
    //imshow("Display Image",img);
    //waitKey(0);
    //glare elimination
    //Mat bgr_image = imread("coinsHarder.jpg");
    //Mat lab_image;
    //cvtColor(bgr_image, lab_image, COLOR_BGR2Lab);

    //// Extract the L channel
    //std::vector<cv::Mat> lab_planes(3);
    //split(lab_image, lab_planes);  // now we have the L image in lab_planes[0]

    //// apply the CLAHE algorithm to the L channel
    //Ptr<cv::CLAHE> clahe = cv::createCLAHE();
    //clahe->setClipLimit(4);
    //cv::Mat dst;
    //clahe->apply(lab_planes[0], dst);

    //// Merge the the color planes back into an Lab image
    //dst.copyTo(lab_planes[0]);
    //cv::merge(lab_planes, lab_image);

    //// convert back to RGB
    //Mat image_clahe;
    //cvtColor(lab_image, image_clahe, COLOR_Lab2BGR);

    //// display the results  (you might also want to see lab_planes[0] before and after).
    //namedWindow("image original", WINDOW_NORMAL);
    //namedWindow("image CLAHE", WINDOW_NORMAL);
    //cv::imshow("image original", bgr_image);
    //cv::imshow("image CLAHE", image_clahe);
    //cv::waitKey(0);
    //destroyAllWindows();
    //Mat image;
    //image = imread("coinsHarder.jpg", IMREAD_GRAYSCALE);
    ///*cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
    //clahe->setClipLimit(4);
    //clahe->setTilesGridSize(Size(2,2));*/
    //cv::Mat dst;
    ////clahe->apply(image, dst);
    ////equalizeHist(image, dst);
    //imwrite("glaredfreed.jpg", dst);
  
}
