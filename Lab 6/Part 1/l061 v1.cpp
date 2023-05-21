/**
 * @file l052.cpp
 * @author Vivian Feng
 * @brief Lab 5, Part 3 Canny Edge detection Commplete
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

#define _USE_MATH_DEFINES

using namespace std::chrono;
using namespace std;

const int PRECISION = 23;
// typedef std::numeric_limits< double > dbl;
// ofstream results;

typedef std::numeric_limits< double > dbl;

class IntegerPoint
{

private:
    int x;
    int y;
    // vector<vector<int>> grid; //should be 800x800 vector
public:
    IntegerPoint(int xp, int yp)
    {
        x = xp;
        y = yp;
        // grid = gr;
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
    int getX()
    {
        return x;
    }
    int getY()
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
    DoublePoint(const DoublePoint &dp){
        x = dp.x;
        y = dp.y;
    }

    void operator=(const DoublePoint &rhs){
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
    int red;
    int green;
    int blue;
    int gray;

public:
    Pixel()
    {
        red = 0;
        blue = 0;
        green = 0;
        gray = 0;
    }

    Pixel(int red, int green, int blue)
    {
        this->red = red;
        this->blue = blue;
        this->green = green;
        this->gray = (red + blue + green) / 3;
    }
    int getGray() const
    {
        return gray;
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
        // cout << "Grid deleted"<<endl;
    }
    void setPixel(int r, int c, int val)
    {
        grid[r][c] = Pixel(val, val, val);
    }
    void setPixel(int r, int c, int red, int green, int blue)
    {
        grid[r][c] = Pixel(red, green, blue);
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
    static void toPPM(vector<vector<int>> matrix, int maxIntensity,string filename){
        ofstream myfile;
        myfile.open(filename);
        myfile << "P3"
               << " " << matrix[0].size() << " " << matrix.size() << " " << maxIntensity << endl;
        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = 0; j < matrix[0].size(); j++)
            {
                myfile << matrix[i][j] << " " << matrix[i][j] << " " << matrix[i][j] << " ";
            }
            myfile << endl;
        }
        myfile.close();
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
        return row;
    }
    int getY() const
    {
        return col;
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
Grid *sobelOperatorPt2(Grid *gr, int low, int high, string of, string f2, string ff, vector<vector<int>> &gx, vector<vector<int>> &gy, vector<vector<double>> &angles);
void areaFill(vector<vector<int>> &gr, int row, int col, unordered_set<Coordinate, Coordinate::HashFunc> &visited);
int roundAngle(double angle);
void computeTangents(vector<vector<int>> &gx, vector<vector<int>> &gy, vector<vector<double>> &result);
// void convertToGrid(vector<vector<int>> *v);

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

Grid *sobelOperatorPt2(Grid *gr, int low, int high, string of, string f2, string ff, vector<vector<int>> &gx, vector<vector<int>> &gy, vector<vector<double>> &angles)
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
    // fill gx
    for (size_t r = 1; r < gx.size() - 1; r++)
    {
        for (size_t c = 1; c < gx[0].size() - 1; c++)
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
            gx[r][c] = val;
        }
    }
    // fill gy
    for (size_t r = 1; r < gy.size() - 1; r++)
    {
        for (size_t c = 1; c < gy[0].size() - 1; c++)
        {
            int val = 0;
            int gystartrow = r - 1;
            int gystartcol = c - 1;
            for (size_t i = 0; i < ydirsobel.size(); i++)
            {
                for (size_t j = 0; j < ydirsobel.size(); j++)
                {
                    int temp = myGrid.getPixel(gystartrow + i, gystartcol + j) * ydirsobel[i][j];
                    val += temp;
                }
            }
            gy[r][c] = val;
        }
    }
    // compute magG
    vector<vector<int>> magG = vector<vector<int>>(myGrid.getHeight(), vector<int>(myGrid.getWidth(), 0));
    for (size_t r = 0; r < magG.size(); r++)
    {
        for (size_t c = 0; c < magG[0].size(); c++)
        {
            magG[r][c] = gx[r][c] * gx[r][c] + gy[r][c] * gy[r][c];
        }
    }
    // grid filled
    int lowsq = low * low;
    int highsq = high * high;
    // apply double threshold, returning pointer to vector on the heap
    vector<vector<int>> edges = vector<vector<int>>(myGrid.getHeight(), vector<int>(myGrid.getWidth(), 0));
    vector<vector<int>> &ed = edges;
    for (size_t r = 0; r < magG.size(); r++)
    {
        for (size_t c = 0; c < magG[0].size(); c++)
        {
            if (magG[r][c] < lowsq)
            {
                ed[r][c] = 0;
            }
            else if (magG[r][c] < highsq)
            {
                ed[r][c] = 1;
            }
            else
            {
                ed[r][c] = 2;
            }
        }
    }

    cout << "Grid filled" << endl;
    // apply area fill to this
    unordered_set<Coordinate, Coordinate::HashFunc> visited = unordered_set<Coordinate, Coordinate::HashFunc>();
    for (size_t r = 0; r < ed.size(); r++)
    {
        for (size_t c = 0; c < ed.size(); c++)
        {
            Coordinate point(r, c);
            if (ed[r][c] == 2)
            {
                areaFill(ed, r, c, visited);
            }
        }
    }
    cout << "Finished areaFill" << endl;

    // apply nonmax suppression
    // compute angles
    // vector<vector<int>> angles = vector<vector<int>>(myGrid.getHeight(), vector<int>(myGrid.getWidth(), 0));
    computeTangents(gx, gy, angles);
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
        }
    }
    // img2.toPPM(f2);
    // place in final grid
    // grid will be filled with 0, 1, 2, 3. This is where part 2 formerly was used to create image1.ppm
    //  Grid newGrid = Grid(myGrid.getWidth(),myGrid.getHeight(),1);
    //  for(int i=0; i < newGrid.getHeight(); i++){
    //      for(int j=0; j < newGrid.getWidth(); j++){
    //          if(ed[i][j] == 0|| ed[i][j] == 1){
    //              newGrid.setPixel(i,j,0);
    //          }else{
    //              newGrid.setPixel(i,j,1);
    //          }
    //      }
    //  }
    //  newGrid.toPPM(of);

    // create final image by doing and of the values in img2 and newGrid
    Grid *finalImage = new Grid(myGrid.getWidth(), myGrid.getHeight(), 1);
    for (int r = 0; r < myGrid.getHeight(); r++)
    {
        for (int c = 0; c < myGrid.getWidth(); c++)
        {
            if (img2.getPixel(r, c) == 1 && (visited.find(Coordinate(r, c)) != visited.end() || ed[r][c] == 2))
            {
                finalImage->setPixel(r, c, 1);
            }
        }
    }
    // newGrid.toPPM(of);
    //finalImage->toPPM(ff);
    return finalImage;
}

/**
 * @brief accordingly changes weak edges in the vector gr to 0 and 2. Marks visited squares with 3
 *
 * @param gr a reference to the vector
 * @param row starting row
 * @param col starting column
 */
void areaFill(vector<vector<int>> &gr, int row, int col, unordered_set<Coordinate, Coordinate::HashFunc> &visited)
{
    // cout << row << " "<< endl;
    Coordinate curr(row, col);
    if (row < 0 || col < 0 || row >= (int)gr.size() || col >= (int)gr[0].size() || visited.find(curr) != visited.end() || gr[row][col] == 0)
    {
        return;
    }
    // for(int i=0; i < gr.size(); i++){
    //     for(int j = 0; j < gr[0].size(); j++){
    //         cout << gr[i][j] << " ";
    //     }
    //     cout << "\n";
    // }
    // cout << "end iteration";
    if (gr[row][col] == 1)
    {
        gr[row][col] = 2;
    }
    else if (gr[row][col] == 2)
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
/**
 * @brief COmpute the angles of the image pixels and puts them in result
 *
 * @param gx
 * @param gy
 * @param result
 */
void computeTangents(vector<vector<int>> &gx, vector<vector<int>> &gy, vector<vector<double>> &result)
{
    for (size_t r = 0; r < gx.size(); r++)
    {
        for (size_t c = 0; c < gx[0].size(); c++)
        {
            double angle = atan2(gy[r][c], gx[r][c]) * 180 / M_PI;
            // int scaled = (int)((angle + 180)*8/360+0.4);
            // result[r][c] = scaled*45-180;
            // if angle <=22.5 and
            // double absangle = angle;
            result[r][c] = angle;
        }
    }
}

//FUNCTION DECLARATIONS
void vote(Grid &edges, vector<vector<int>> &voteMatrix,vector<vector<double>> angles);
void voteLine(int r,int c, vector<vector<int>> &voteMatrix,vector<vector<double>> angle);
void drawPixel(vector<vector<int>> &matrix, int r,int c);
void drawLine(vector<vector<int>> &voteMatrix, IntegerPoint &p1, IntegerPoint &p2);
void drawDAxPlus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2);
void drawDAyPlus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2);
void drawDAxMinus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2);
void drawDAyMinus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2);
void drawHoriz(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2);
void drawVert(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2);


/**
 * @brief Votes along the lines of pixels for each edge
 * 
 * @param gr 
 * @param voteMatrix 
 */
void vote(Grid &edges, vector<vector<int>> &voteMatrix,vector<vector<double>> angles){
    for(int r=0; r < edges.getHeight();r++){
        for(int c=0; c < edges.getWidth();c++){
            if(edges.getPixel(r,c)==1)
                voteLine(r,c,voteMatrix,angles);
        }
    }

}
void voteLine(int r,int c, vector<vector<int>> &voteMatrix,vector<vector<double>> angle){
    //find the start and end of the line originating from r,c in direction of angle
    double y = r*1.0/voteMatrix.size();
    double x = c*1.0/voteMatrix[0].size();
    //draw using bresheram's algorithm
    if(abs(abs(angle[r][c])-90) < 1e-7){
        double y0 = 0;
        double y1 = 1;
        IntegerPoint i0 = DoublePoint(x,0).scaleToInt(voteMatrix.size(),voteMatrix[0].size());
        IntegerPoint i1 = DoublePoint(x,1).scaleToInt(voteMatrix.size(),voteMatrix[0].size());
        drawLine(voteMatrix,i0,i1);

    }else{
        double slope = tan(angle[r][c]);
        double x0 = 1/slope* (-y) + x;
        double x1 = 1/slope*(1-x)+x;
        IntegerPoint i0 = DoublePoint(x0,0).scaleToInt(voteMatrix.size(),voteMatrix[0].size());
        IntegerPoint i1 = DoublePoint(x1,1).scaleToInt(voteMatrix.size(),voteMatrix[0].size());
        drawLine(voteMatrix, i0,i1);

    }

    
}    
void drawPixel(vector<vector<int>> &matrix, int r,int c){
    if(r >= 0 && r < matrix.size() && c >=0 && c < matrix[0].size()){
        matrix[r][c] += 1;
        //cout << "P";
    }
}
/*
    draws a line given two pixel coordinates on the ppm
*/
void drawLine(vector<vector<int>> &voteMatrix, IntegerPoint &p1, IntegerPoint &p2)
{
    // determine case
    int dx = p2.getX() - p1.getX();
    int dy = p2.getY() - p1.getY();
    if (abs(dx) >= abs(dy))
    {
        if (dx > 0 && dy > 0)
        {
            drawDAxPlus(voteMatrix,p1.getX(), p1.getY(), p2.getX(), p2.getY());
        }
        else if (dx < 0 && dy < 0)
        {
            drawDAxPlus(voteMatrix, p2.getX(), p2.getY(), p1.getX(), p1.getY());
        }
        else if (dy == 0)
        {
            drawHoriz(voteMatrix,p1.getX(), p1.getY(), p2.getX(), p2.getY());
        }
        else if (dx > 0 && dy < 0)
        {
            drawDAxMinus(voteMatrix,p1.getX(), p1.getY(), p2.getX(), p2.getY());
        }
        else if (dx < 0 && dy > 0)
        {
            drawDAxMinus(voteMatrix, p2.getX(), p2.getY(), p1.getX(), p1.getY());
        }
    }
    else
    {
        if (dx > 0 && dy > 0)
        {
            drawDAyPlus(voteMatrix, p1.getX(), p1.getY(), p2.getX(), p2.getY());
        }
        else if (dx < 0 && dy < 0)
        {
            drawDAyPlus(voteMatrix, p2.getX(), p2.getY(), p1.getX(), p1.getY());
        }
        else if (dx == 0)
        {
            drawVert(voteMatrix, p1.getX(), p1.getY(), p2.getX(), p2.getY());
        }
        else if (dx < 0 && dy > 0)
        {
            drawDAyMinus(voteMatrix, p1.getX(), p1.getY(), p2.getX(), p2.getY());
        }
        else if (dx > 0 && dy < 0)
        {
            drawDAyMinus(voteMatrix, p2.getX(), p2.getY(), p1.getX(), p1.getY());
        }
    }
    // draw appropriately
}
/*
    Precondition: x1 < x2, y1 < y2, and dx > dy
    xmult and ymult are transformations applied to drawing other cases
*/
void drawDAxPlus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    int j = y1;
    int eps = dy - dx;
    for (int i = x1; i <= x2; i++)
    {
        drawPixel(voteMatrix, i, j);
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
void drawDAyPlus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    int j = x1;
    int eps = dx - dy;
    for (int i = y1; i <= y2; i++)
    {
        drawPixel(voteMatrix,j, i);
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
void drawDAxMinus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    int eps = abs(dy) - abs(dx);
    int absdy = abs(dy);
    int j = y1;
    for (int i = x1; i <= x2; i++)
    {
        drawPixel(voteMatrix,i, j);
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
void drawDAyMinus(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    int eps = abs(dx) - abs(dy);
    int absdx = abs(dx);
    int j = x1;
    for (int i = y1; i <= y2; i++)
    {
        drawPixel(voteMatrix,j, i);
        if (eps >= 0)
        {
            j -= 1;
            eps -= dy;
        }
        eps += absdx;
    }
}
void drawHoriz(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2)
{
    for (int i = min(x1, x2); i <= max(x1, x2); i++)
    {
        drawPixel(voteMatrix,i, y1);
    }
}
void drawVert(vector<vector<int>> &voteMatrix,int x1, int y1, int x2, int y2)
{
    for (int i = min(y2, y1); i <= max(y2, y1); i++)
    {
        drawPixel(voteMatrix,x1, i);
    }
}
void part1(const int &argc, char **&argv)
{
    string imageName = "imgs.ppm";
    int low = 150;
    int high = 175;
    string fg = "imageg.ppm";
    string of = "image1.ppm";
    string f2 = "image2.ppm";
    string ff = "imagef.ppm";
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
        else if (flag == "-of")
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
        i++;
    }
    // print the values of the flags
    cout << imageName << " " << low << " " << high << " " << fg << " " << of << " " << f2 << " " << ff << endl;
    // read in the image
    Grid *gr = readImage(imageName);
    // gr->toPPM(fg);
    // create the version with edge detection
    vector<vector<double>> angles = vector<vector<double>>(gr->getHeight(), vector<double>(gr->getWidth(), 0));
    vector<vector<int>> gx = vector<vector<int>>(gr->getHeight(), vector<int>(gr->getWidth(), 0));
    vector<vector<int>> gy = vector<vector<int>>(gr->getHeight(), vector<int>(gr->getWidth(), 0));
    Grid* edgebinary = sobelOperatorPt2(gr, low, high, of, f2, ff,gx,gy,angles);
    edgebinary->toPPM("imagef.ppm");
    cout << "Edge dection finished\n";
    // apply voting and threshold, write to original image
    vector<vector<int>> votes = vector<vector<int>>(gr->getHeight(), vector<int>(gr->getWidth(), 0));
    vote(*edgebinary, votes,angles);
    cout << "Voting finished\n";
    int maxIntensity = 0;
    for(int i=0; i < votes.size(); i++){
        for(int j=0; j < votes[0].size(); j++){
            maxIntensity = max(votes[i][j],maxIntensity);
        }
    }
    Grid::toPPM(votes,maxIntensity,"imagev.ppm");
    cout << "Votes image created "<<maxIntensity <<endl;
    //Grid::toPPM()
    // create the voting image
    // TODO: standalone PPM, color PPM, port over breshnam
    delete gr;
    delete edgebinary;
}
int main(int argc, char **argv)
{
    // part1();
    auto start = high_resolution_clock::now();
    part1(argc, argv);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "Time taken " << duration.count() << endl;
}