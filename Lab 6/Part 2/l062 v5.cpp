/**
 * @file l061.cpp
 * @author Vivian Feng
 * @brief Lab 6, Part 1 Coin Dectection
 * @version 5.3
 * @date 2023-1-20
 *
 * @copyright Copyright (c) 2023
 *
 */
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
// #include <numbers>

#define _USE_MATH_DEFINES

using namespace std::chrono;
using namespace std;

const int PRECISION = 23;
const double ERROR = 1e-7;
int CIRCLE_RADIUS = 120;
// typedef std::numeric_limits< double > dbl;
// ofstream results;

typedef std::numeric_limits<double> dbl;

class IntegerPoint
{

private:
    int x;
    int y;
    // vector<vector<int>> grid; //should be 800x800 vector
public:
    // IntegerPoint(){
    //     x = 0;
    //     y = 0;
    // }
    IntegerPoint(int xp, int yp)
    {
        x = xp;
        y = yp;
        // grid = gr;
    }
    IntegerPoint()
    {
        x = 0;
        y = 0;
    }
    void operator=(const IntegerPoint &rhs)
    {
        x = rhs.x;
        y = rhs.y;
    }
    string toString(bool fullPrecision = false)
    {
        if (!fullPrecision)
        {
            return "(" + to_string(x) + "," + to_string(y) + ")";
        }
        else
        {
            stringstream stream;
            stream.precision(dbl::max_digits10);
            stream << fixed << "(" << x << "," << y << ")";
            return stream.str();
        }
    }
    int getX() const
    {
        return x;
    }
    int getY() const
    {
        return y;
    }
    void setX(int xp)
    {
        x = xp;
    }
    void setY(int yp)
    {
        y = yp;
    }
};
class DoublePoint
{
private:
    double x;
    double y;

public:
    DoublePoint()
    {
        x = 0;
        y = 0;
    }
    DoublePoint(double xp, double yp)
    {
        x = xp;
        y = yp;
    }
    DoublePoint(const DoublePoint &dp)
    {
        x = dp.x;
        y = dp.y;
    }

    void operator=(const DoublePoint &rhs)
    {
        x = rhs.x;
        y = rhs.y;
    }
    void setX(double xp)
    {
        x = xp;
    }
    void setY(double yp)
    {
        y = yp;
    }
    double getX()
    {
        return x;
    }
    double getY()
    {
        return y;
    }
    // void set(double xp, double yp){
    //     x = xp;
    //     y = yp;
    // }
    ~DoublePoint() {}
    string toString(bool fullPrecision = false)
    {
        if (!fullPrecision)
        {
            return ("(" + to_string(x)) + "," + to_string(y) + ")";
        }
        else
        {
            stringstream stream;
            stream.precision(17);
            stream << fixed << "(" << x << "," << y << ")";
            return stream.str();
        }
    }
    IntegerPoint scaleToInt(int height, int width)
    {
        IntegerPoint intpoint((int)(round(width * x)), (int)(round(height * y)));
        return intpoint;
    }
    // DoublePoint operator/ (const DoublePoint &lhs,double num){
    //     return DoublePoint(lhs/num,)
    // }
    // void printPoint(std::ofstream outlog){
    //     outlog.precision(dbl::max_digits10);
    //     outlog << "(" << x <<","<<y<<")";
    // }
};

/**
 * @brief Pixel class necessary to restore the original image
 *
 */
class Pixel
{
private:
    vector<int> radii;
    int red;
    int green;
    int blue;
    int gray;
    bool isCenter;
    

public:
    Pixel()
    {
        red = 0;
        blue = 0;
        green = 0;
        gray = 0;
        isCenter = false;
    }

    Pixel(int red, int green, int blue)
    {
        this->red = red;
        this->blue = blue;
        this->green = green;
        this->gray = (red + blue + green) / 3;
        isCenter = false;
    }
    bool getCenter() const{
        return isCenter;
    }
    void setCenter(bool val){
        isCenter = val;
    }
    int getGray() const
    {
        return gray;
    }
    void setGray(int val){
        gray = val;
    }
    int getRed()
    {
        return red;
    }
    int getGreen()
    {
        return green;
    }
    int getBlue()
    {
        return blue;
    }
    vector<int> getRadii(){
        return radii;
    }
    void addRadii(int r){
        radii.push_back(r);
    }
    Pixel(const Pixel &pix)
    {
        this->red = pix.red;
        this->blue = pix.blue;
        this->green = pix.green;
        this->gray = pix.gray;
    }
};

class Grid
{
private:
    int width;
    int height;
    int maxIntensity;
    vector<vector<Pixel>> grid;

public:
    Grid(int width, int height, int maxIntensity)
    {
        grid = vector<vector<Pixel>>(height, vector<Pixel>(width, Pixel()));
        // grid = new int *[mySize];
        // size = mySize;
        this->height = height;
        this->width = width;
        this->maxIntensity = maxIntensity;
    }
    Grid(const Grid &gr)
    {
        grid = gr.grid;
        width = gr.width;
        height = gr.height;
    }
    ~Grid()
    {
        // std::cout << "Grid deleted"<<endl;
    }
    void setPixel(int r, int c, int val)
    {
        if (r >= 0 && r < height && c >= 0 && c < width)
        {
            grid[r][c].setGray(val);
        }
    }
    void setPixel(int r, int c, int red, int green, int blue)
    {
        if (r >= 0 && r < height && c >= 0 && c < width)
        {
            grid[r][c] = Pixel(red, green, blue);
            // return true;
        }
        // return false;
    }
    void drawPixel(int x, int y, int red, int green, int blue)
    {
        setPixel(y, x, red, green, blue);
        // return outcome;
    }
    int getPixel(int r, int c)
    {
        return grid[r][c].getGray();
    }
    int getHeight()
    {
        return height;
    }
    int getWidth()
    {
        return width;
    }
    void markCenter(int r,int c,bool val){
        grid[r][c].setCenter(val);
    }
    bool getCenter(int r, int c){
        return grid[r][c].getCenter();
    }
    void toPPM(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "P3"
               << " " << width << " " << height << " " << maxIntensity << endl;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                myfile << grid[i][j].getGray() << " " << grid[i][j].getGray() << " " << grid[i][j].getGray() << " ";
            }
            myfile << endl;
        }
        myfile.close();
    }
    void toColorPPM(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "P3"
               << " " << width << " " << height << " " << maxIntensity << endl;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                myfile << grid[i][j].getRed() << " " << grid[i][j].getGreen() << " " << grid[i][j].getBlue() << " ";
            }
            myfile << endl;
        }
        myfile.close();
    }
    void addRadii(int r, int c,int radius){
        grid[r][c].addRadii(radius);
    }
    int getMaxIntensity()
    {
        return maxIntensity;
    }
    void drawCircle(IntegerPoint intCenter, int radius, int red, int green, int blue)
    {
        // int size = gr.getSize();
        // IntegerPoint intCenter = center.scaleToInt(size);

        // midpoint circle algorithm

        int r = radius;
        // if (customRad > 0)
        // {
        //     r = customRad;
        // }
        int x, y, y2, y2_new, ty;
        // vector<IntegerPoint> lastPts = vector<IntegerPoint>();
        //  xmax = (int)(round(radius*size * 0.70710678)); // maximum x at radius/sqrt(2) (go to 46 degrees )
        y = r;
        y2 = y * y;
        ty = (2 * y) - 1;
        // ty = 3-2*y;
        y2_new = y2;
        // int yfinal = y;
        for (x = 0; x <= y; x++)
        {
            if ((y2 - y2_new) >= ty)
            {
                y2 -= ty;
                y -= 1;
                ty -= 2;
            }
            // shift center from given algorithm
            drawPixel(x + intCenter.getX(), y + intCenter.getY(), red, green, blue);
            // cout << x+intCenter.getX()<<" "<<y+intCenter.getY()<<endl;
            drawPixel(x + intCenter.getX(), -y + intCenter.getY(), red, green, blue);
            drawPixel(-x + intCenter.getX(), y + intCenter.getY(), red, green, blue);
            drawPixel(-x + intCenter.getX(), -y + intCenter.getY(), red, green, blue);
            drawPixel(y + intCenter.getX(), x + intCenter.getY(), red, green, blue);
            drawPixel(y + intCenter.getX(), -x + intCenter.getY(), red, green, blue);
            drawPixel(-y + intCenter.getX(), x + intCenter.getY(), red, green, blue);
            drawPixel(-y + intCenter.getX(), -x + intCenter.getY(), red, green, blue);
            y2_new -= (2 * x) - 3;
            // yfinal = y;
            // cout << gr.toString()<<endl;
        }
    }
};


/**
 * @brief Returns a pointer to a Grid Object on the heap containing the image's intensity matrix
 *
 * @param filename
 * @return Grid*
 */
Grid *readImage(string filename)
{
    // open image file
    /*
    std::ifstream file("filename");
    std::vector<int> array;
    int number;
    while(file >> number) {
        array.push_back(number);
    }

    */
    ifstream imgFile(filename);
    string enc;
    imgFile >> enc;
    int width, height;
    imgFile >> width >> height;
    int intensity;
    imgFile >> intensity;
    Grid *gr = new Grid(width, height, intensity);
    for (int r = 0; r < height; r++)
    {
        for (int c = 0; c < width; c++)
        {
            int red, blue, green;
            imgFile >> red >> blue >> green;
            // int gray = (red + blue + green) / 3;
            gr->setPixel(r, c, red, green, blue);
        }
    }
    CIRCLE_RADIUS = max(120,min(width / 10, height / 10));
    imgFile.close();
    return gr;
}

// PART 2 starts here
/**
 * @brief Hashable object
 *
 */
class Coordinate
{
private:
    int row;
    int col;

public:
    Coordinate()
    {
        row = 0;
        col = 0;
    }
    Coordinate(int row, int col)
    {
        this->row = row;
        this->col = col;
    }
    bool operator==(const Coordinate &other) const
    {
        return this->row == other.row && this->col == other.col;
    }
    int getX() const
    {
        return col;
    }
    int getY() const
    {
        return row;
    }
    class HashFunc
    {
    public:
        size_t operator()(const Coordinate &point) const
        {
            size_t xHash = std::hash<int>()(point.getX());
            size_t yHash = std::hash<int>()(point.getY()) << 1;
            return xHash ^ yHash;
        }
    };
};
Grid *sobelOperatorPt2(Grid *gr, int low, int high, string of, string f2, string ff, vector<vector<double>> &angles);
void areaFill(Grid &gr, int row, int col, unordered_set<Coordinate, Coordinate::HashFunc> &visited);
int roundAngle(double angle);

// void convertToGrid(vector<vector<int>> *v);
void gaussBlur(Grid &gr, double sigma, int k){
    vector<vector<double>> filter = vector<vector<double>>(2*k+1, vector<double>(2*k+1,0));
    //fill filter
    double multiplier = 1.0/(2*M_PI*sigma*sigma);
    int size = 2*k+1;
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            double expr = -((i-(k+1))*(i-(k+1)) + (j-(k+1))*(j-(k+1)))/(2.0*sigma*sigma);
            filter[i][j] = exp(expr)*multiplier;
        }
    }




    for(int r=2; r < gr.getHeight()-2; r++){
        for(int c=2; c < gr.getWidth()-2; c++){
            // int pix = gr.getPixel(r,c);
            int startrow = r-k;
            int startcol = c-k;
            double val = 0;
            for (size_t i = 0; i < filter.size(); i++)
            {
                for (size_t j = 0; j < filter.size(); j++)
                {
                    double temp = gr.getPixel(startrow + i, startcol + j) * filter[i][j];
                    val += temp;
                }
            }
            gr.setPixel(r,c,(int)val);
        }
    }
    cout << "Applied Gaussian Blur\n";
}


/**
 * @brief Applies the sobel operator and double thereshold and hysterisis and nonmax suppression
 *
 * @param gr pointer to grid
 * @param low low threshold
 * @param high high threshold
 * @param of name for image1.ppm
 * @param f2 name for image2.ppm
 * @param ff name of final image
 */
Grid *sobelOperatorPt2(Grid *gr, int low, int high, string of, string f2, string ff, vector<vector<double>> &angles)
{
    Grid &myGrid = *gr;
    
    // vector<vector<int>> gx = vector<vector<int>>(myGrid.getHeight(), vector<int>(myGrid.getWidth(), 0));
    // vector<vector<int>> gx = vector<vector<int>>(myGrid.getHeight(), vector<int>(myGrid.getWidth(), 0));
    vector<vector<int>> xdirsobel =
        {{-1, 0, 1},
         {-2, 0, 2},
         {-1, 0, 1}};
    vector<vector<int>> ydirsobel =
        {{-1, -2, -1},
         {0, 0, 0},
         {1, 2, 1}};

    vector<vector<int>> magG = vector<vector<int>>(myGrid.getHeight(), vector<int>(myGrid.getWidth(), 0));
    // fill gx
    for (int r = 1; r < myGrid.getHeight() - 1; r++)
    {
        for (int c = 1; c < myGrid.getWidth() - 1; c++)
        {
            int val = 0;
            int gxstartrow = r - 1;
            int gxstartcol = c - 1;
            for (size_t i = 0; i < xdirsobel.size(); i++)
            {
                for (size_t j = 0; j < xdirsobel.size(); j++)
                {
                    int temp = myGrid.getPixel(gxstartrow + i, gxstartcol + j) * xdirsobel[i][j];
                    val += temp;
                }
            }
            int gx = val;
            int gy = 0;
            int gystartrow = r - 1;
            int gystartcol = c - 1;
            for (size_t i = 0; i < ydirsobel.size(); i++)
            {
                for (size_t j = 0; j < ydirsobel.size(); j++)
                {
                    int temp = myGrid.getPixel(gystartrow + i, gystartcol + j) * ydirsobel[i][j];
                    gy += temp;
                }
            }
            magG[r][c] = gx * gx + gy * gy;
            angles[r][c] = atan2(gy, gx) * 180 / M_PI;
            // cout << angles[r][c] << endl;
        }
    }
    // grid filled
    int lowsq = low * low;
    int highsq = high * high;
    // apply double threshold, returning pointer to vector on the heap
    Grid *finalImage = new Grid(myGrid.getWidth(), myGrid.getHeight(), 1);
    // Grid &edges = *finalImage;
    Grid &ed = *finalImage;
    for (int r = 0; r < myGrid.getHeight(); r++)
    {
        for (int c = 0; c < myGrid.getWidth(); c++)
        {
            if (magG[r][c] < lowsq)
            {
                ed.setPixel(r, c, 0);
            }
            else if (magG[r][c] < highsq)
            {
                ed.setPixel(r, c, 2);
            }
            else
            {
                ed.setPixel(r, c, 1);
            }
        }
    }

    std::cout << "Grid filled" << endl;
    // apply area fill to this
    unordered_set<Coordinate, Coordinate::HashFunc> visited = unordered_set<Coordinate, Coordinate::HashFunc>();
    for (int r = 0; r < ed.getHeight(); r++)
    {
        for (int c = 0; c < ed.getWidth(); c++)
        {
            Coordinate point(r, c);
            if (ed.getPixel(r, c) == 1)
            {
                areaFill(ed, r, c, visited);
            }
        }
    }
    std::cout << "Finished areaFill" << endl;
    // place in final grid
    // grid will be filled with 0, 1, 2, 3. This is where part 2 formerly was used to create image1.ppm
    Grid &newGrid = ed;
    for (int i = 0; i < newGrid.getHeight(); i++)
    {
        for (int j = 0; j < newGrid.getWidth(); j++)
        {
            if (ed.getPixel(i, j) == 0 || ed.getPixel(i, j) == 2)
            {
                newGrid.setPixel(i, j, 0);
            }
            else
            {
                newGrid.setPixel(i, j, 1);
            }
        }
    }
    //  newGrid.toPPM(of);
    //  cout << "created image1.ppm\n";
    // apply nonmax suppression
    // compute angles
    // vector<vector<int>> angles = vector<vector<int>>(myGrid.getHeight(), vector<int>(myGrid.getWidth(), 0));
    // computeTangents(gx, gy, angles);
    // produce image2
    Grid img2 = Grid(myGrid.getWidth(), myGrid.getHeight(), 1);
    for (size_t r = 1; r < angles.size() - 1; r++)
    {
        for (size_t c = 1; c < angles[0].size() - 1; c++)
        {
            int p1 = magG[r][c];
            int p0 = 0;
            int p2 = 0;
            int roundedAngle = roundAngle(angles[r][c]);
            if (abs(roundedAngle) == 0 || abs(roundedAngle) == 180)
            {
                p0 = magG[r][c - 1];
                p2 = magG[r][c + 1];
            }
            else if (roundedAngle == 90 || roundedAngle == -90)
            {
                p0 = magG[r - 1][c];
                p2 = magG[r + 1][c];
            }
            else if (roundedAngle == -45 || roundedAngle == 135)
            {
                p0 = magG[r - 1][c + 1];
                p2 = magG[r + 1][c - 1];
            }
            else if (roundedAngle == 45 || roundedAngle == -135)
            {
                p0 = magG[r - 1][c - 1];
                p2 = magG[r + 1][c + 1];
            }
            if (p1 >= p0 && p1 >= p2)
            {
                img2.setPixel(r, c, 1);
            }

            // modify edges final Image
            if (img2.getPixel(r, c) == 1 && (visited.find(Coordinate(r, c)) != visited.end()))
            {
                finalImage->setPixel(r, c, 1);
            }
            else
            {
                finalImage->setPixel(r, c, 0);
                // cout << "Changed";
            }
        }
    }
    // img2.toPPM(f2);

    // create final image by doing and of the values in img2 and newGrid

    // for (int r = 0; r < myGrid.getHeight(); r++)
    // {
    //     for (int c = 0; c < myGrid.getWidth(); c++)
    //     {
    //         if (img2.getPixel(r, c) == 1 && (visited.find(Coordinate(r, c)) != visited.end()))
    //         {
    //             finalImage->setPixel(r, c, 1);
    //         }else{
    //             finalImage->setPixel(r,c,0);
    //             //cout << "Changed";
    //         }
    //     }
    // }
    // newGrid.toPPM(of);
    // finalImage->toPPM(ff);
    return finalImage;
}

/**
 * @brief accordingly changes weak edges in the vector gr to 0 and 2. Marks visited squares with 3
 *
 * @param gr a reference to the vector
 * @param row starting row
 * @param col starting column
 */
void areaFill(Grid &gr, int row, int col, unordered_set<Coordinate, Coordinate::HashFunc> &visited)
{
    // std::cout << row << " "<< endl;
    Coordinate curr(row, col);
    if (row < 0 || col < 0 || row >= gr.getHeight() || col >= gr.getWidth() || visited.find(curr) != visited.end() || gr.getPixel(row, col) == 0)
    {
        return;
    }
    // for(int i=0; i < gr.size(); i++){
    //     for(int j = 0; j < gr[0].size(); j++){
    //         std::cout << gr[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << "end iteration";
    if (gr.getPixel(row, col) == 2)
    {
        gr.setPixel(row, col, 1);
    }
    else if (gr.getPixel(row, col) == 1)
    {
        // Coordinate point(row,col);
        visited.insert(curr);
    }
    areaFill(gr, row + 1, col, visited);
    areaFill(gr, row + 1, col + 1, visited);
    areaFill(gr, row, col + 1, visited);
    areaFill(gr, row - 1, col + 1, visited);
    areaFill(gr, row - 1, col, visited);
    areaFill(gr, row - 1, col - 1, visited);
    areaFill(gr, row, col - 1, visited);
    areaFill(gr, row + 1, col - 1, visited);
}

// LAB 6 PART 1 STARTS HERE
/**
 * @brief Rounds an angle
 *
 */
int roundAngle(double angle)
{
    if (-180 <= angle && angle <= -157.5)
    {
        return -180;
    }
    else if (-157.5 < angle && angle <= -112.5)
    {
        return -135;
    }
    else if (-112.5 < angle && angle <= -67.5)
    {
        return -90;
    }
    else if (-67.5 < angle && angle <= 22.5)
    {
        return -45;
    }
    else if (-22.5 < angle && angle <= 22.5)
    {
        return 0;
    }
    else if (22.5 < angle && angle <= 67.5)
    {
        return 45;
    }
    else if (67.5 < angle && angle <= 112.5)
    {
        return 90;
    }
    else if (112.5 < angle && angle <= 167.5)
    {
        return 135;
    }
    else
    {
        return 180;
    }
}

// Class VoteMatrix
/**
 * @brief Matrix for votes. X is columns, Y is rows!!
 *
 */
class VoteMatrix
{
private:
    vector<vector<int>> grid;
    int height;
    int width;
    int maxIntensity;

public:
    VoteMatrix(int height, int width)
    {
        grid = vector<vector<int>>(height, vector<int>(width, 0));
        this->height = height;
        this->width = width;
        maxIntensity = 0;
    }
    void drawPixel(IntegerPoint &p1)
    {
        if (p1.getX() < width && p1.getY() < height && p1.getX() >= 0 && p1.getY() >= 0)
        {
            grid[p1.getY()][p1.getX()] += 1;
            // grid[p1.getY()][p1.getX()]= 1;
            maxIntensity = max(maxIntensity, grid[p1.getY()][p1.getX()]);
        }
    }
    // void drawPixel(int )

    void drawPixel(int x, int y)
    {
        if (x < width && y < height && x >= 0 && y >= 0)
        {
            grid[y][x] += 1;
            // grid[y][x]=1;
            maxIntensity = max(grid[y][x], maxIntensity);
        }
    }
    int getPixel(int r, int c)
    {
        return grid[r][c];
    }
    int getHeight() const
    {
        return height;
    }
    int getWidth() const
    {
        return width;
    }
    int getMaxIntensity()
    {
        return maxIntensity;
    }
    void voteAlongLine(int x, int y, double angle)
    {
        double dblx = x * 1.0 / width;
        double dbly = y * 1.0 / height;
        double x0 = 0; // y = 0
        double x1 = 0; // y=1
        double y0 = 0; // x = 0
        double y1 = 0; // x = 1;
        double slope = tan(angle * M_PI / 180);
        // std::cout << dblx << " " << dbly << " " << slope << endl;
        // std::cout << slope << endl;
        if (abs(abs(angle) - 90) < ERROR)
        { // if vertical line
            x0 = dblx;
            x1 = dblx;
            IntegerPoint i0 = DoublePoint(x0, 0).scaleToInt(height, width);
            // std::cout << "i0" << i0.toString() << "\n";
            IntegerPoint i1 = DoublePoint(x1, 1).scaleToInt(height, width);
            // std::cout << "i0: " << i0.toString() << "\n";
            // std::cout << "i1 " << i1.toString() << endl;
            drawLine(i0, i1);
        }
        else if (abs(angle - 180) < ERROR || abs(angle) < ERROR || abs(angle - (-180)) < ERROR)
        { // if horizontal line
            y0 = dbly;
            y1 = dbly;
            IntegerPoint i0 = DoublePoint(0, y0).scaleToInt(height, width);
            // std::cout << "i0" << i0.toString() << "\n";
            IntegerPoint i1 = DoublePoint(1, y1).scaleToInt(height, width);
            // std::cout << "i0: " << i0.toString() << "\n";
            // std::cout << "i1 " << i1.toString() << endl;
            drawLine(i0, i1);
        }
        else
        {

            double magn = sqrt(1 + slope * slope);
            double unitdx = 1 / magn;
            double unitdy = slope / magn;
            // total change in X: scale unit dx by radius/width
            // total change in Y: scal unitdx by radius/height
            double scaledRadiusx = CIRCLE_RADIUS * 1.0 / width;
            double scaledRadiusy = CIRCLE_RADIUS * 1.0 / height;
            double deltax = unitdx * scaledRadiusx;
            double deltay = unitdy * scaledRadiusy;
            DoublePoint start = DoublePoint(dblx + deltax, dbly + deltay);
            DoublePoint end = DoublePoint(dblx - deltax, dbly - deltay);
            IntegerPoint startInt = start.scaleToInt(height, width);
            IntegerPoint endInt = end.scaleToInt(height, width);
            drawLine(startInt, endInt);
            // drawLine(scaledX,0,0,scaledY);
            // IntegerPoint i1 = d1.scaleToInt(size);
            // IntegerPoint i2 = d2.scaleToInt(size);
            // drawLine(i1, i2);
        }
    }

    void drawLine(int x1, int y1, int x2, int y2)
    {
        IntegerPoint p1(x1, y1);
        IntegerPoint p2(x2, y2);
        drawLine(p1, p2);
    }

    /*
        draws a line given two pixel coordinates on the ppm
    */
    void drawLine(IntegerPoint &p1, IntegerPoint &p2)
    {
        // determine case
        int dx = p2.getX() - p1.getX();
        int dy = p2.getY() - p1.getY();
        if (abs(dx) >= abs(dy))
        {
            if (dx > 0 && dy > 0)
            {
                drawDAxPlus(p1.getX(), p1.getY(), p2.getX(), p2.getY());
            }
            else if (dx < 0 && dy < 0)
            {
                drawDAxPlus(p2.getX(), p2.getY(), p1.getX(), p1.getY());
            }
            else if (dy == 0)
            {
                drawHoriz(p1.getX(), p1.getY(), p2.getX(), p2.getY());
            }
            else if (dx > 0 && dy < 0)
            {
                drawDAxMinus(p1.getX(), p1.getY(), p2.getX(), p2.getY());
            }
            else if (dx < 0 && dy > 0)
            {
                drawDAxMinus(p2.getX(), p2.getY(), p1.getX(), p1.getY());
            }
        }
        else
        {
            if (dx > 0 && dy > 0)
            {
                drawDAyPlus(p1.getX(), p1.getY(), p2.getX(), p2.getY());
            }
            else if (dx < 0 && dy < 0)
            {
                drawDAyPlus(p2.getX(), p2.getY(), p1.getX(), p1.getY());
            }
            else if (dx == 0)
            {
                drawVert(p1.getX(), p1.getY(), p2.getX(), p2.getY());
            }
            else if (dx < 0 && dy > 0)
            {
                drawDAyMinus(p1.getX(), p1.getY(), p2.getX(), p2.getY());
            }
            else if (dx > 0 && dy < 0)
            {
                drawDAyMinus(p2.getX(), p2.getY(), p1.getX(), p1.getY());
            }
        }
        // draw appropriately
    }
    /*
        Precondition: x1 < x2, y1 < y2, and dx > dy
        xmult and ymult are transformations applied to drawing other cases
    */
    void drawDAxPlus(int x1, int y1, int x2, int y2)
    {
        int dx = x2 - x1;
        int dy = y2 - y1;
        int j = y1;
        int eps = dy - dx;
        for (int i = x1; i <= x2; i++)
        {
            drawPixel(i, j);
            if (eps >= 0)
            {
                j += 1;
                eps -= dx;
            }
            eps += dy;
        }
    }
    /*
        Precondition dy > dx, dy > 0, dx > 0
    */
    void drawDAyPlus(int x1, int y1, int x2, int y2)
    {
        int dx = x2 - x1;
        int dy = y2 - y1;
        int j = x1;
        int eps = dx - dy;
        for (int i = y1; i <= y2; i++)
        {
            drawPixel(j, i);
            if (eps >= 0)
            {
                j += 1;
                eps -= dy;
            }
            eps += dx;
        }
    }
    /*
        Precondition: x1 < x2, y1 > y2
    */
    void drawDAxMinus(int x1, int y1, int x2, int y2)
    {
        int dx = x2 - x1;
        int dy = y2 - y1;
        int eps = abs(dy) - abs(dx);
        int absdy = abs(dy);
        int j = y1;
        for (int i = x1; i <= x2; i++)
        {
            drawPixel(i, j);
            if (eps >= 0)
            {
                j -= 1;
                eps -= dx;
            }
            eps += absdy;
        }
    }
    /*
Precondition: x1 > x2, y1 < y2
*/
    void drawDAyMinus(int x1, int y1, int x2, int y2)
    {
        int dx = x2 - x1;
        int dy = y2 - y1;
        int eps = abs(dx) - abs(dy);
        int absdx = abs(dx);
        int j = x1;
        for (int i = y1; i <= y2; i++)
        {
            drawPixel(j, i);
            if (eps >= 0)
            {
                j -= 1;
                eps -= dy;
            }
            eps += absdx;
        }
    }
    void drawHoriz(int x1, int y1, int x2, int y2)
    {
        for (int i = min(x1, x2); i <= max(x1, x2); i++)
        {
            drawPixel(i, y1);
        }
    }
    void drawVert(int x1, int y1, int x2, int y2)
    {
        for (int i = min(y2, y1); i <= max(y2, y1); i++)
        {
            drawPixel(x1, i);
        }
    }
    /*
        String representation of the grid. Useful for small representations
    */
    string toString()
    {
        string toRet = "";
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                toRet += to_string(grid[i][j]) + " ";
                // std::cout << grid[i][j];
            }
            toRet += "\n";
        }
        return toRet;
    }
    /**
     * @brief Keeps top percent of votes
     *
     * @param kernelsize -- the size of half the window
     * @param percent
     * @param threshold
     */
    void cleanUp(int kernelsize, double percent, int threshold, unordered_set<Coordinate, Coordinate::HashFunc> &centers,bool override=false)
    {
        int filtersize = 2 * kernelsize + 1;
        cout << "Filter size: " << filtersize;
        // vector<Coordinate> centers =
        for (size_t i = kernelsize; i < grid.size() - kernelsize; i += kernelsize)
        {
            for (size_t j = kernelsize; j < grid[0].size() - kernelsize; j += kernelsize)
            {
                int startr = i - kernelsize;
                int startc = j - kernelsize;
                int q1threshold = maxIntensity * 0.25;
                int q2threshold = maxIntensity * 0.25;
                int q3threshold = maxIntensity * 0.35;
                int q4threshold = maxIntensity *0.35;
                vector<int> localCounts = vector<int>();
                for (int y = startr; y < startr + filtersize; y++)
                {
                    for (int x = startc; x < startc + filtersize; x++)
                    {
                        // cout << grid[y][x];
                        if(override){
                            if (getPixel(y, x) >= threshold)
                            {
                                // Coordinate coord(y,x);
                                // cout << "HEllo";
                                localCounts.push_back(grid[y][x]);
                            }
                        }else{
                            

                            /*
                            """"""""""
                            *I  | II *
                            *III| IV *
                            * *******
                            * 
                            */
                            if(y < grid.size()/2){
                                if(x < grid[0].size()/2){
                                    if(getPixel(y,x) >= q1threshold){
                                        localCounts.push_back(grid[y][x]);
                                    }
                                }else{
                                    if(getPixel(y,x) >= q2threshold){
                                        localCounts.push_back(grid[y][x]);
                                    }
                                }
                                
                            }else{
                                if(x < grid[0].size()/2){
                                    if(getPixel(y,x) >= q3threshold){
                                        localCounts.push_back(grid[y][x]);
                                    }
                                }else{
                                    if(getPixel(y,x) >= q4threshold){
                                        localCounts.push_back(grid[y][x]);
                                    }
                                }
                            }
                        }
                        
                    }
                }
                std::sort(localCounts.begin(), localCounts.end(), std::greater<int>());
                if (localCounts.size() == 0)
                    continue;
                int localThreshold = localCounts[int(percent * localCounts.size())];
                // int localThreshold = localCounts[0];
                // int localThreshold=threshold;
                for (int r = startr; r < startr + filtersize; r++)
                {
                    for (int c = startc; c < startc + filtersize; c++)
                    {
                        if (grid[r][c] >= localThreshold)
                        {
                            centers.insert(Coordinate(r, c));
                        }
                        else
                        {
                            grid[r][c] = 0;
                        }
                    }
                }
            }
        }
    }
    void top10(unordered_set<Coordinate, Coordinate::HashFunc> &centers, int filtersize){
        int totalsize = filtersize*2+1;
        for(int r = filtersize; r < grid.size()-filtersize; r+=filtersize){
            for(int c=filtersize; c < grid[0].size()-filtersize; c+=filtersize){
                int localMax = 0;
                for(int y=r-filtersize; y < r+filtersize; y++){
                    for(int x=c-filtersize; x < c+filtersize; x++){
                        localMax = max(localMax,grid[y][x]);
                    }
                }

                int thresh = localMax *0.9;
                for(int y=r-filtersize; y < r+filtersize; y++){
                    for(int x=c-filtersize; x < c+filtersize; x++){
                        if(grid[y][x] > thresh){
                            
                        }else{
                            Coordinate coord(y,x);
                            if(centers.find(coord)!=centers.end()){
                                centers.erase(coord);
                            }
                        }
                    }
                }
            }
        }
    }
    void toPPM(int maxIntensity, string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "P3"
               << " " << grid[0].size() << " " << grid.size() << " " << maxIntensity << endl;
        for (size_t i = 0; i < grid.size(); i++)
        {
            for (size_t j = 0; j < grid[0].size(); j++)
            {
                myfile << grid[i][j] << " " << grid[i][j] << " " << grid[i][j] << " ";
            }
            myfile << endl;
        }
        myfile.close();
    }
};

void vote(Grid &edges, VoteMatrix &votes, vector<vector<double>> &angles)
{
    int count = 0;
    for (int r = 0; r < edges.getHeight(); r++)
    {
        for (int c = 0; c < edges.getWidth(); c++)
        {
            if (edges.getPixel(r, c) == 1)
            {
                votes.voteAlongLine(c, r, angles[r][c]);
                count += 1;
            }
            // if(count==1){

            //     std::cout << "Row, COl" << r << " " << c << " angle: " << angles[r][c]<<endl;
            //     return;
            // }
        }
    }
}


class Circle
{
private:
    Coordinate center;
    int radius;
    double matchFraction;
    int edgeCount;

public:
    Circle()
    {
        center = Coordinate(0, 0);
        radius = 0;
        matchFraction = 0;
        edgeCount = 0;
    }
    Circle(const Coordinate &cent, int rad)
    {
        center = cent;
        radius = rad;
        matchFraction = 0;
    }
    Circle(const Coordinate &cent, int rad,double frac, int edge)
    {
        center = cent;
        radius = rad;
        matchFraction = frac;
        edgeCount = edge;
    }
    Circle(const Circle &c)
    {
        center = c.center;
        radius = c.radius;
        matchFraction = c.matchFraction;
    }
    Coordinate getCenter() const{
        return center;
    }
    int getRadius() const{
        return radius;
    }
    void setMatchFraction(double f){
        matchFraction = f;
    }
    int getEdgeCount() const{
        return edgeCount;
    }
    double getMatchFraction() const{
        //cout << matchFraction << endl;
        return matchFraction;
    }
    int squareddistanceToCenter(const Circle &other){
        int dx = other.center.getX() - center.getX();
        int dy = other.center.getY() - this->center.getY();
        return dx * dx + dy*dy;
    }
    static bool compareMatchFraction(const Circle &a, const Circle &b){
        return a.getMatchFraction() > b.getMatchFraction();
    }
};

void medianBlur(Grid &image,int filtersize){
    int totalsize = 2*filtersize-1;
    for(int i=filtersize; i < image.getHeight()-filtersize; i++){
        for(int j=filtersize; j < image.getWidth()-filtersize; j++){
            vector<int> vals = vector<int>();
            for(int r=i-filtersize; r < i+filtersize; r++){
                for(int c=j-filtersize; c < j+filtersize; c++){
                    vals.push_back(image.getPixel(r,c));
                }
                
            }
            std::sort(vals.begin(),vals.end());
            int median = vals[vals.size()/2];
            image.setPixel(i,j,median);
        }
    }
}
bool match(Grid &edges, int x, int y){
    if(x >= 0 && x < edges.getWidth() && y >= 0 && y < edges.getHeight()){
        return edges.getPixel(y,x);
    }
    return false;
}
void findCircles(Grid &image,unordered_set<Coordinate, Coordinate::HashFunc> centers, Grid &edges,vector<Circle> &circs)
{
    int minRad = 60;//min(60, edges.getHeight() / 50);
    int maxRad = max(120, edges.getHeight() / 25);
    cout << "Min radius, max radius " << minRad << " " << maxRad << endl;
    int size = maxRad - minRad + 1;
    for (const Coordinate &cent : centers)
    {
        vector<int> counts = vector<int>(size, 0);
        vector<int> edgcnts = vector<int>(size,0);
        for (int radius = minRad; radius <= maxRad; radius++)
        {
            int countsIdx = radius - minRad;
            int x, y, y2, y2_new, ty;
            int totalCount = 0;
            int xCC = cent.getX();
            int yCC = cent.getY();
            y = radius;
            y2 = y * y;
            ty = 2 * y - 1;
            y2_new = y2;
            
            for (x = 0; x <= y; x++)
            {
                if ((y2 - y2_new) >= ty)
                {
                    y2 -= ty;
                    y -= 1;
                    ty -= 2;
                }
                // shift center from given algorithm
                if (match(edges,x + xCC, y + yCC))
                {
                    counts[countsIdx] += 1;
                }
                if (match(edges,x + xCC, -y + yCC))
                {
                    counts[countsIdx] += 1;
                }
                if (match(edges,-x + xCC, y + yCC))
                {
                    counts[countsIdx] += 1;
                }
                if (match(edges,-x + xCC, -y + yCC))
                {
                    counts[countsIdx] += 1;
                }
                if (match(edges,y + xCC, x + yCC) )
                {
                    counts[countsIdx] += 1;
                }
                if (match(edges,y + xCC, -x + yCC))
                {
                    counts[countsIdx] += 1;
                }
                if (match(edges,-y + xCC, x + yCC))
                {
                    counts[countsIdx] += 1;
                }
                if (match(edges,-y + xCC, -x + yCC))
                {
                    counts[countsIdx] += 1;
                }
                // cout << x+intCenter.getX()<<" "<<y+intCenter.getY()<<endl;
                totalCount += 8;
                y2_new -= (2 * x) - 3;
                // yfinal = y;
                // cout << gr.toString()<<endl;
            }
            // int percentageThreshold = 0.20;
            // double decimalthreshold = percentageThreshold/100.0;
            // int edgeCountThreshold = decimalthreshold * totalCount;
            // if(counts[countsIdx] > edgeCountThreshold){
            //     circs.push_back(Circle(cent,radius));
            //     // cout << "in here";
            // }
            edgcnts[countsIdx] = totalCount;
        }
        int percentageThreshold = 20;
        double decimalthreshold = percentageThreshold/100.0;
        
        int binhalf = 1;
        int binsize = 2*binhalf+1;
        for(int r=minRad+1; r < maxRad; r+=1 ){
            int idx = r-minRad;
            int bincount = counts[idx] + counts[idx-1];
            // for(int i=idx-binhalf; i <= idx+binhalf; i++){
            //     bincount += counts[i];
            // }
            if(bincount > decimalthreshold * edgcnts[idx]){
                double matchFrac=bincount*1.0/edgcnts[idx];
                circs.push_back(Circle(cent,r,matchFrac,bincount));
                image.markCenter(cent.getY(), cent.getX(),true);
            }
        }
        
    }
}
void eliminateFalseCircles(Grid &gr, vector<Circle> &circles,int halfkernelsize,vector<Circle> &finalCircles){
    unordered_set<int> idxToRemove = unordered_set<int>();
    //for circle in circles: 
    int distanceThreshold = 70;
    int distsqred = distanceThreshold * distanceThreshold;
    for(size_t i=0; i < circles.size(); i++){
        //if circle is not to be removed
        if(idxToRemove.find(i)==idxToRemove.end()){
            //for circle in other circles:
                //if circle comes within radius of other center
                    //mark the one with the lower match fraction for deletion
                    //when finished -- continue
            for(size_t j=i+1; j < circles.size();j++){
                Circle &c1 = circles[i];
                Circle &c2 = circles[j];
                int dist = c1.squareddistanceToCenter(c2);
                if(dist < distsqred){
                        // int maxEdgeCount = max(c1.getEdgeCount(),c2.getEdgeCount());
                        // int minEdgeCount = min(c1.getEdgeCount(),c2.getEdgeCount());
                        // if(minEdgeCount < 0.1*maxEdgeCount){
                        //     if(c1.getEdgeCount() < c2.getEdgeCount()){
                        //         idxToRemove.insert(j);
                        //     }else{
                        //         idxToRemove.insert(i);
                        //     }
                        // }else{
                        if(c1.getMatchFraction() < c2.getMatchFraction()){
                            idxToRemove.insert(i);
                        }
                        else{
                            idxToRemove.insert(j);
                        }
                        //}
                    
                    // if(c1.getMatchFraction() < c2.getMatchFraction()){
                    //     //if c1 < 10 of max of edgeCount:
                    //         //do i2
                    //     //else:
                    //         //do i1;
                    //     int maxEdgeCount = max(c1.getEdgeCount(),c2.getEdgeCount());
                    //     int minEdgeCount = min(c1.getEdgeCount(),c2.getEdgeCount());
                    //     if(minEdgeCount < 0.1*maxEdgeCount){
                    //         if(c1.getEdgeCount() < c2.getEdgeCount()){
                    //             idxToRemove.insert(j);
                    //         }else{
                    //             idxToRemove.insert(i);
                    //         }
                    //     }else{
                    //         idxToRemove.insert(i);
                    //     }
                        
                    // }else{
                    //     idxToRemove.insert(j);
                    // }
                }


            }
        }
    }
    for(int k=0; k < circles.size(); k++){
        if(idxToRemove.find(k)==idxToRemove.end()){
            finalCircles.push_back(circles[k]);
        }
    }
    
        
        

}
/**
 * @brief Part 2 of lab 1
 *
 * @param argc
 * @param argv
 */
void part2(const int &argc, char **&argv)
{
    string imageName = "image.ppm";
    int low = 50;
    int high = 75;
    string fg = "imageg.ppm";
    string of = "image1.ppm";
    string f2 = "image2.ppm";
    string ff = "imagef.ppm";
    string fv = "imagev.ppm";
    string fcc = "imageCC.ppm";
    int thresholdCount = 25;
    bool override = false;
    // string imageName = "coinsPiece.ppm";
    // int low = 50;
    // int high = 75;
    // string fg = "imageg.ppm";
    // string of = "image1.ppm";
    // string f2 = "image2.ppm";
    // string ff = "imagef.ppm";
    // string fv = "imagev.ppm";
    // string fcc = "imageCC.ppm";
    // int thresholdCount = 25;
    

    // read the command line arguments
    for (int i = 1; i < argc; i++)
    {
        string flag = string(argv[i]);
        if (flag == "-f")
        {
            // cout << "flag read" << endl;
            imageName = argv[i + 1];
        }
        else if (flag == "-lt")
        {
            stringstream ss(argv[i + 1]);
            ss >> low;
        }
        else if (flag == "-ht")
        {
            stringstream ss(argv[i + 1]);
            ss >> high;
        }
        else if (flag == "-f1")
        {
            of = argv[i + 1];
        }
        else if (flag == "-f2")
        {
            f2 = argv[i + 1];
        }
        else if (flag == "-ff")
        {
            ff = argv[i + 1];
        }
        else if (flag == "-fg")
        {
            fg = argv[i + 1];
        }
        else if (flag == "-TC")
        {
            stringstream ss(argv[i + 1]);
            ss >> thresholdCount;
            override= true;
        }
        else if (flag == "-fv")
        {
            stringstream ss(argv[i + 1]);
            ss >> fv;
        }
        else if (flag == "-fcc")
        {
            stringstream ss(argv[i + 1]);
            ss >> fcc;
        }
        i++;
    }
    // print the values of the flags
    std::cout << "-f " << imageName << "-lt " << low << " -ht " << high << "-of " << of << "-f2 " << f2 << "-ff " << ff << "-fg " << fg << " -TC " << thresholdCount << "-fv " << fv << "-fcc " << fcc << endl;
    // read in the image
    Grid *gr = readImage(imageName);
    
    cout << "Image read\n";
    //medianBlur(*gr,1);
    gaussBlur(*gr,1.4,2);
    //gr->toPPM(fg);
    cout << "Gray image created\n";
    // create the version with edge detection
    vector<vector<double>> angles = vector<vector<double>>(gr->getHeight(), vector<double>(gr->getWidth(), 0));
    Grid *edgebinary = sobelOperatorPt2(gr, low, high, of, f2, ff, angles);
    std::cout << "Edge dection finished\n";
    edgebinary->toPPM(ff);

    // apply voting and threshold, write to original image
    // vector<vector<int>> votes = vector<vector<int>>(gr->getHeight(), vector<int>(gr->getWidth(), 0));
    VoteMatrix votes = VoteMatrix(edgebinary->getHeight(), edgebinary->getWidth());
    vote(*edgebinary, votes, angles);
    // votes.voteAlongLine(120,120,45);
    //  for(int i=-180; i<=180; i++){
    //      votes.voteAlongLine(120,120,i);
    //  }
    // std::cout << "Voting finished\n";
    int maxIntensity = votes.getMaxIntensity();
    // for (int i = 0; i < votes.getHeight(); i++)
    // {
    //     for (int j = 0; j < votes.getWidth(); j++)
    //     {
    //         maxIntensity = max(votes.getPixel(i,j), maxIntensity);
    //     }
    // }
    cout << "Voting finished, max intensity = " << maxIntensity << endl;
    votes.toPPM(maxIntensity, fv);
    cout << "Vote image created\n";
    unordered_set<Coordinate, Coordinate::HashFunc> centers = unordered_set<Coordinate, Coordinate::HashFunc>();

    votes.cleanUp(5, 0.1, thresholdCount, centers,override);
    votes.top10(centers,7);
    cout << "\nlen of centers " << centers.size() << endl;
    for (const Coordinate &center : centers)
    {
        for (int i = 0; i < 4; i++)
        {
            gr->drawCircle(IntegerPoint(center.getX(), center.getY()), i, gr->getMaxIntensity(), 0, 0);
        }
    }
    //gr->toColorPPM(fcc);
    vector<Circle> circles = vector<Circle>();
    findCircles(*gr,centers,*edgebinary,circles);
    cout << "len of circles " << circles.size() << endl;
    vector<Circle> finalCircles = vector<Circle>();
    eliminateFalseCircles(*gr,circles,20,finalCircles);
    for (const Circle &c : finalCircles)
    {
        int radius = c.getRadius();
        Coordinate center = c.getCenter();
        for (int i = radius-1; i <= radius+1; i++)
        {
            gr->drawCircle(IntegerPoint(center.getX(), center.getY()), i, 0, gr->getMaxIntensity(), 0);
            
        }
        
        //gr->addRadii(center.getY(),center.getX(),radius);
    }
    
    cout<< "len of finalCircles" << finalCircles.size() << endl;
    gr->toColorPPM("imageCircles.ppm");
    delete gr;
    delete edgebinary;
}
int main(int argc, char **argv)
{
    // part1();
    auto start = high_resolution_clock::now();
    //part1(argc, argv);
    part2(argc, argv);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    std::cout << "Time taken " << duration.count() << endl;
}