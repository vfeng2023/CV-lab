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

using namespace std::chrono;
using namespace std;

const int PRECISION = 23;
// typedef std::numeric_limits< double > dbl;
// ofstream results;
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
            stream.precision(PRECISION);
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
    DoublePoint(const DoublePoint &dp)
    {
        x = dp.x;
        y = dp.y;
    }
    DoublePoint(double xp, double yp)
    {
        x = xp;
        y = yp;
    }
    DoublePoint(DoublePoint &&dp)
    {
        x = dp.x;
        y = dp.y;
    }

    void operator=(const DoublePoint &rhs)
    {
        x = rhs.x;
        y = rhs.y;
    }
    void operator=(DoublePoint &&rhs)
    {
        x = rhs.x;
        y = rhs.y;
    }
    // bool operator==(DoublePoint &rhs){
    //     return x == rhs.x && y==rhs.y;
    // }
    // bool operator!=(DoublePoint)
    void setX(double xp)
    {
        x = xp;
    }
    void setY(double yp)
    {
        y = yp;
    }
    double getX() const
    {
        return x;
    }
    double getY() const
    {
        return y;
    }
    void set(double xp, double yp)
    {
        x = xp;
        y = yp;
    }
    void set(DoublePoint &dp)
    {
        x = dp.x;
        y = dp.y;
    }
    DoublePoint operator-(DoublePoint const &p2)
    {
        double retX = x - p2.x;
        double retY = y - p2.y;
        return DoublePoint(retX, retY);
    }
    // DoublePoint operator-(DoublePoint p1,DoublePoint const &p2){

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
            stream.precision(PRECISION);
            stream << fixed << "(" << x << "," << y << ")";
            return stream.str();
        }
    }

    IntegerPoint scaleToInt(int size)
    {
        IntegerPoint intpoint((int)(round(size * x)), (int)(round(size * y)));
        return intpoint;
    }
    bool static comparePoint(DoublePoint dp1, DoublePoint dp2)
    {
        return dp1.getX() < dp2.getX();
    }
    bool static comparePointY(DoublePoint dp1, DoublePoint dp2)
    {
        return dp1.getY() < dp2.getY();
    }
};

// ostream& operator<<(ostream& os,DoublePoint& dt){
//     os << dt.toString(true);
//     return os;
// }

class Grid
{
private:
    int size;
    int **grid;

public:
    Grid(int mySize)
    {
        // grid = vector<vector<int>>(mySize, vector<int>(mySize, 1));
        grid = new int *[mySize];
        size = mySize;
        for (int i = 0; i < mySize; i++)
        {
            grid[i] = new int[mySize];
        }
        for (int r = 0; r < mySize; r++)
        {
            for (int c = 0; c < mySize; c++)
            {
                grid[r][c] = 1;
            }
        }
    }
    Grid(const Grid &gr)
    {
        grid = new int *[gr.size];
        size = gr.size;
        for (int i = 0; i < size; i++)
        {
            grid[i] = new int[size];
        }
        for (int r = 0; r < size; r++)
        {
            for (int c = 0; c < size; c++)
            {
                grid[r][c] = gr.grid[r][c];
            }
        }
    }
    ~Grid()
    {
        for (int i = 0; i < size; i++)
        {
            delete[] grid[i];
        }
        delete[] grid;
        // cout << "Grid deleted"<<endl;
    }
    void drawPixel(IntegerPoint &p1)
    {
        if (p1.getX() < size && p1.getY() < size && p1.getX() >= 0 && p1.getY() >= 0)
            grid[size - 1 - p1.getY()][p1.getX()] = 0;
    }
    /*
        Choices in color - 0 = black, 1 = white, 2 = red, 3 = green, 4 = blue
    */
    void drawPixel(IntegerPoint &p1, int color)
    {
        if (p1.getX() < size && p1.getY() < size && p1.getX() >= 0 && p1.getY() >= 0)
            grid[size - 1 - p1.getY()][p1.getX()] = color;
    }
    void drawPixel(int x, int y, int color = 0)
    {
        if (x < size && y < size && x >= 0 && y >= 0)
            grid[size - 1 - y][x] = color;
    }
    // void drawPixel(int x,int y,int c){
    //     if (x < size && y < size && x >=0 && y >= 0)
    //         grid[x][y] = 2;
    // }
    int getSize()
    {
        return size;
    }
    void drawLine(int x1, int y1, int x2, int y2)
    {
        IntegerPoint p1(x1, y1);
        IntegerPoint p2(x2, y2);
        drawLine(p1, p2);
    }
    /*
        Given two unscaled points, draws line on grid that goes through these two points
    */
    void drawLine(DoublePoint &d1, DoublePoint &d2)
    {
        double dx = d2.getX() - d1.getX();
        double dy = d2.getY() - d1.getY();
        if (abs(dx) < 1e-8)
        {
            DoublePoint start(d1.getX(), 0);
            DoublePoint end(d1.getX(), 1);
            IntegerPoint i1 = start.scaleToInt(size);
            IntegerPoint i2 = end.scaleToInt(size);
            drawLine(i1, i2);
            cout << "\nslope: vertical line";
            return;
        }
        double slope = dy / dx;
        if (abs(slope) < 1e-8)
        {
            DoublePoint start(0, d1.getY());
            DoublePoint end(1, d1.getY());
            IntegerPoint i1 = start.scaleToInt(size);
            IntegerPoint i2 = end.scaleToInt(size);
            drawLine(i1, i2);
            cout << "\nslope: horizontal line" << endl;
            return;
        }
        // line equation: slope * (x-d1.getX()) + d1.getY()
        //  double yintercept = slope *(-d1.getX())+d1.getY();
        //  double xintercept = -d1.getY()/slope + d1.getX();
        //  int scaledY = (int)(round(yintercept*size));
        //  int scaledX = (int)(round(xintercept*size));

        // cout <<"y intercept: "<<yintercept<<endl;
        // cout << "xintercept: "<<xintercept<<endl;
        cout << "\nslope " << slope << endl;
        // find intersection with all sides of the bounding box: x = 0, x=1, y=0, y=1
        // take the points which have a positive coordinate
        double yx0 = slope * (-1.0 * d1.getX()) + d1.getY();
        double yx1 = slope * (1.0 - d1.getX()) + d1.getY();
        double xy0 = -d1.getY() / slope + d1.getX();
        double xy1 = (1.0 - d1.getY()) / slope + d1.getX();
        bool chosestart = false;

        DoublePoint start(0.0, yx0);
        DoublePoint end(1, yx1);
        if (yx0 >= 0 && yx0 <= 1)
        {
            if (!chosestart)
            {
                start.setX(0);
                start.setY(yx0);
                chosestart = true;
            }
            else
            {
                end.setX(0);
                end.setY(yx0);
            }
        }
        if (yx1 >= 0 && yx1 <= 1)
        {
            if (!chosestart)
            {
                start.setX(1);
                start.setY(yx1);
                chosestart = true;
            }
            else
            {
                end.setX(1);
                end.setY(yx1);
            }
        }
        if (xy0 >= 0 && xy0 <= 1)
        {
            if (!chosestart)
            {

                start.setX(xy0);
                start.setY(0);
                chosestart = true;
            }
            else
            {
                end.setX(xy0);
                end.setY(0);
            }
        }
        if (xy1 >= 0 && xy1 <= 1)
        {
            if (!chosestart)
            {
                start.setX(xy1);
                start.setY(1.0);
                chosestart = true;
            }
            else
            {
                end.setX(xy1);
                end.setY(1.0);
            }
        }
        cout << "\n"
             << start.toString() << endl;
        cout << end.toString() << endl;
        // DoublePoint start(0.0, y0);
        // DoublePoint end(1, y1);
        IntegerPoint startInt = start.scaleToInt(size);
        IntegerPoint endInt = end.scaleToInt(size);
        drawLine(startInt, endInt);
        // drawLine(scaledX,0,0,scaledY);
        IntegerPoint i1 = d1.scaleToInt(size);
        IntegerPoint i2 = d2.scaleToInt(size);
        drawLine(i1, i2);
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
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                toRet += to_string(grid[i][j]) + " ";
                // cout << grid[i][j];
            }
            toRet += "\n";
        }
        return toRet;
    }
    /*
        Outputs ppm file contents. Change the output to be triangle.ppm like this: ./a.out > outputfile.txt
    */
    void toPPM(string filename = "test.ppm")
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "P3"
               << " " << size << " " << size << " " << 1 << endl;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if (grid[i][j] == 1)
                {
                    myfile << "1 1 1 ";
                }
                else if (grid[i][j] == 0)
                {
                    myfile << "0 0 0 ";
                }
                else if (grid[i][j] == 2)
                {
                    myfile << "1 0 0 ";
                }
                else if (grid[i][j] == 3)
                {
                    myfile << "0 1 0 ";
                }
                else if (grid[i][j] == 4)
                {
                    myfile << "0 0 1 ";
                }
            }
            myfile << endl;
        }
        myfile.close();
    }
};
/*
    takes two DoublePoints. Given a third point, returns a line which is perpendicular
    can draw the line on the grid using grid's drawLine method
*/
class Line
{
private:
    DoublePoint p1;
    DoublePoint p2;
    double slope;
    double ycept;
    bool isDef;

public:
    Line() {}
    Line(double x1, double y1, double x2, double y2)
    {
        p1 = DoublePoint(x1, y1);
        p2 = DoublePoint(x2, y2);
        double dx = x2 - x1;
        double dy = y2 - y1;
        if (abs(dx) < 1e-7)
        {
            isDef = false;
        }
        else
        {
            slope = dy / dx;
            ycept = slope * (-x1) + y1;
            isDef = true;
        }
    }
    Line(const Line &line)
    {
        p1 = line.p1;
        p2 = line.p2;
        slope = line.slope;
        ycept = line.ycept;
        isDef = line.isDef;
    }
    Line(DoublePoint &u, DoublePoint &v)
    {
        // p1 = u;
        // p2 = v;
        p1.setX(u.getX());
        p1.setY(u.getY());
        p2.setX(v.getX());
        p2.setY(v.getY());
        double dx = p2.getX() - p1.getX();
        double dy = p2.getY() - p1.getY();
        if (abs(dx) < 1e-7)
        {
            isDef = false;
        }
        else
        {
            slope = dy / dx;
            ycept = slope * (-p1.getX()) + p1.getY();
            isDef = true;
        }
    }
    // Line(Line &line){
    //     p1 = D;
    //     p2 = DoublePoint(line.p2);
    // }
    /*
        Returns a line that is perpendicular to this line
    */
    double getSlope()
    {
        if (isDef)
        {
            return slope;
        }
        else
        {
            return 10000;
        }
    }
    double getIntercept()
    {
        return ycept;
    }
    Line findPerpFromPoint(DoublePoint &perp)
    {
        // slope of current line
        double dx = p2.getX() - p1.getX();
        double dy = p2.getY() - p1.getY();
        if (abs(dx) < 1e-7)
        { // vertical line
            return Line(p1.getX(), perp.getY(), perp.getX(), perp.getY());
        }
        else if (abs(dy) < 1e-7)
        {
            return Line(perp.getX(), perp.getY(), perp.getX(), p1.getY());
        }
        // otherwise, find slope of line, take negative reciprocal return owo
        double mySlope = dy / dx;
        double perpSlope = -1 / mySlope;                         // slope
        double ycept = perpSlope * (-perp.getX()) + perp.getY(); // yinterceipt
        return Line(perp.getX(), perp.getY(), 0.0, ycept);
    }
    /*
        Finds the value of x at a particular point on the line
    */
    double evaluateAt(double x)
    {
        if (isDef)
        {
            return slope * x + ycept;
        }
        else
        {
            return 999999999;
        }
    }
    /*
     Returns a DoublePoint that is the intersection of the line with another line. If there is no intersection, returns (0.0,0.0)
    */
    DoublePoint findIntersect(Line &other)
    { // returns o
        if (!other.isDef && !isDef)
        {
            return DoublePoint(0, 0);
        }
        else if (!isDef)
        {
            double yp = other.evaluateAt(p1.getX());
            return DoublePoint(p1.getX(), yp);
        }
        else if (!other.isDef)
        {
            double yp = evaluateAt(other.p1.getX());
            return DoublePoint(other.p1.getX(), yp);
        }
        else if (abs(other.slope - slope) < 1e-7)
        {
            return DoublePoint(0, 0);
        }
        double x = (other.ycept - ycept) / (slope - other.slope);
        return DoublePoint(x, evaluateAt(x));
    }
    double calcDist()
    {
        double d2 = (p1.getX() - p2.getX()) * (p1.getX() - p2.getX()) + (p1.getY() - p2.getY()) * (p1.getY() - p2.getY());
        return sqrt(d2);
    }
    double static calcDist(const DoublePoint &d1, const DoublePoint &d2)
    {
        double dist2 = (d1.getX() - d2.getX()) * (d1.getX() - d2.getX()) + (d1.getY() - d2.getY()) * (d1.getY() - d2.getY());
        return sqrt(dist2);
    }
    double static calcDist2(DoublePoint d1, DoublePoint d2)
    {
        double dist2 = (d1.getX() - d2.getX()) * (d1.getX() - d2.getX()) + (d1.getY() - d2.getY()) * (d1.getY() - d2.getY());
        return sqrt(dist2);
    }
    void draw(Grid &gr)
    {
        gr.drawLine(p1, p2);
    }
};
class Circle
{
private:
    DoublePoint center;
    double radius;

public:
    Circle() {}
    Circle(Circle &c)
    {
        center = c.center;
        radius = c.radius;
    }
    Circle(double x, double y, double rad)
    {
        center = DoublePoint(x, y);
        radius = rad;
    }
    ~Circle()
    {
        // delete center;
        //  cout << "Circle deleted"<<endl;
    }
    DoublePoint getCenter()
    {
        return center;
    }
    double getRadius()
    {
        return radius;
    }
    string toString()
    {
        return "Center: " + center.toString() + "; radius=" + to_string(radius);
    }
    /*
        Draws the circle on a grid;
    */
    void drawCircle(Grid &gr, int color = 0, int customRad = -1)
    {
        int size = gr.getSize();
        IntegerPoint intCenter = center.scaleToInt(size);

        // midpoint circle algorithm

        int r = (int)(round(radius * size));
        if (customRad > 0)
        {
            r = customRad;
        }
        int x, y, y2, y2_new, ty;
        vector<IntegerPoint> lastPts = vector<IntegerPoint>();
        // xmax = (int)(round(radius*size * 0.70710678)); // maximum x at radius/sqrt(2) (go to 46 degrees )
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
            gr.drawPixel(x + intCenter.getX(), y + intCenter.getY(), color);
            // cout << x+intCenter.getX()<<" "<<y+intCenter.getY()<<endl;
            gr.drawPixel(x + intCenter.getX(), -y + intCenter.getY(), color);
            gr.drawPixel(-x + intCenter.getX(), y + intCenter.getY(), color);
            gr.drawPixel(-x + intCenter.getX(), -y + intCenter.getY(), color);
            gr.drawPixel(y + intCenter.getX(), x + intCenter.getY(), color);
            gr.drawPixel(y + intCenter.getX(), -x + intCenter.getY(), color);
            gr.drawPixel(-y + intCenter.getX(), x + intCenter.getY(), color);
            gr.drawPixel(-y + intCenter.getX(), -x + intCenter.getY(), color);
            y2_new -= (2 * x) - 3;
            // yfinal = y;
            // cout << gr.toString()<<endl;
        }
    }
};

vector<DoublePoint> splitpt2(string filename)
{
    fstream myFile;
    myFile.open(filename);
    vector<DoublePoint> toRet;
    string line;
    while (getline(myFile, line))
    {
        stringstream ss(line);
        vector<string> coords;
        string buffer;
        while (ss >> buffer)
        {
            if (!buffer.empty())
                coords.push_back(buffer);
        }
        double x = stod(coords[0]);
        double y = stod(coords[1]);
        DoublePoint dp(x, y);
        toRet.push_back(dp);
    }
    // cout << toRet[0].toString();
    myFile.close();
    return toRet;
}
void generatePoints(int numPoints)
{

    srand(time(NULL));
    ofstream pointsFile;
    pointsFile.open("points.txt");
    pointsFile.precision(PRECISION);
    for (int i = 0; i < numPoints; i++)
    {
        double x = rand() * 1.0 / RAND_MAX;
        double y = rand() * 1.0 / RAND_MAX;
        pointsFile << fixed << x << "  " << y << endl;
    }
    pointsFile.close();
    //}
}
/*
    draws a bolded circle
*/
void highLightPoint(Grid &gr, DoublePoint &p, int rad, int color = 2)
{
    Circle c(p.getX(), p.getY(), rad * 1.0 / gr.getSize());
    c.drawCircle(gr, color, rad);
    c.drawCircle(gr, color, rad - 1);
}

/*
    Part 1 is due Monday 12/05

 a) name your file l041.cpp (lower case L followed by the digits 0,4,1)

 b) in the main you should have only a call to part1 method (with arguments or not, is your choice)

 c) generate 60 random points in the unit square (save them in a data structure of your chosing) and save them in points.txt in the same format as 3.4 (each line has the x and y coordinate of a point separated by 2 spaces, use at least 20 digits precision)

 d) apply the QuickHull algorithm to find a convex hull that encloses all 60 points and the vertices are points from the set of 60 points you generated

 e) create a ppm of size 400x400 named quickhull.ppm in which you display a circle of radius 3 for each scaled point you generated, as well as connect the vertices you obtained for the convex hull (so in the image I should see 60 points and a convex hull with vertices circles of radus 3 that are conected forming a convex hull)

  f) complete the attached form Project 4 Part 1 QuickHull.docx and turn in this document as well as the cpp file

*/
/**
 https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
*/
/**
 * Given two points a and b which lie on a line, and a point c, tells if point c is to right of oriented line AB
 */
bool isRight(DoublePoint &a, DoublePoint &b, DoublePoint &c)
{
    return ((b.getX() - a.getX()) * (c.getY() - a.getY()) - (b.getY() - a.getY()) * (c.getX() - a.getX())) < 0;
}
double calcDistfromLine(DoublePoint &l1, DoublePoint &l2, DoublePoint &p)
{
    DoublePoint u = p - l1;
    DoublePoint v = l2 - l1;
    double magv2 = v.getX() * v.getX() + v.getY() * v.getY();
    double magu2 = u.getX() * u.getX() + u.getY() * u.getY();
    double compSquared = pow((u.getX() * v.getX()) + (u.getY() * v.getY()), 2) / magv2;
    double d2 = magu2 - compSquared;

    return sqrt(d2);
}
/*
    points -- contains the og points
    subset -- contains the indices of the poins to the right of ab
    idx1 and idx2 -- the indices of A & B in the convex hull array poly
*/

void findHull(vector<DoublePoint> &points, vector<int> &subset, list<int> &poly, list<int>::iterator &idx1, list<int>::iterator &idx2)
{
    // if no points return
    int aidx = *idx1;
    int bidx = *idx2;
    // cout << "\n" << aidx << " "<< bidx;
    // cout << " subset: ";
    // for(int j:subset){
    //     cout << j << " ";
    // }
    // cout << "\n";
    if (subset.size() == 0)
    {
        return;
    }

    // find furthest point from ab
    int smallestIdx = subset[0];
    double maxDist = calcDistfromLine(points[aidx], points[bidx], points[smallestIdx]);
    for (int idx : subset)
    {
        double currDist = calcDistfromLine(points[aidx], points[bidx], points[idx]);
        if (currDist > maxDist)
        {
            smallestIdx = idx;
            maxDist = currDist;
        }
    }
    // insert between ab
    poly.insert(idx2, smallestIdx);
    auto newpos = prev(idx2, 1);
    // call findhull again
    vector<int> s1;
    vector<int> s2;
    for (int pointsIdx : subset)
    {
        if (pointsIdx == aidx || pointsIdx == bidx || pointsIdx == smallestIdx)
            continue;
        if (isRight(points[aidx], points[smallestIdx], points[pointsIdx]))
        {
            s1.push_back(pointsIdx);
        }
        else if (isRight(points[smallestIdx], points[bidx], points[pointsIdx]))
        {
            s2.push_back(pointsIdx);
        }
    }

    // cout << *idx1 << " "<< *newpos;
    findHull(points, s1, poly, idx1, newpos);
    findHull(points, s2, poly, newpos, idx2);
}
/*
    Precondition: points is nonzero, hull is empty
*/
list<int> quickHull(vector<DoublePoint> &points)
{
    // find left and right most points
    int aidx = 0; // leftmost point
    for (size_t i = 0; i < points.size(); i++)
    {
        if (points[i].getX() < points[aidx].getX())
        {
            aidx = i;
        }
        else if (abs(points[i].getX() - points[aidx].getX()) < 1e-15)
        {
            if (points[i].getY() < points[aidx].getY())
            {
                aidx = i;
            }
        }
    }
    // cout << "Aidx" << aidx;
    int bidx = 0; // leftmost point
    for (size_t i = 0; i < points.size(); i++)
    {
        if (points[i].getX() > points[bidx].getX())
        {
            bidx = i;
        }
        else if (abs(points[i].getX() - points[bidx].getX()) < 1e-15)
        {
            if (points[i].getY() < points[bidx].getY())
            {
                bidx = i;
            }
        }
    }
    // cout << "bidx" << bidx;
    //  if equal, use min y as tiebreaker
    list<int> polygon = {aidx, bidx, aidx};
    // define s1 and s2
    vector<int> s1;
    vector<int> s2;

    for (int j = 0; j < (int)points.size(); j++)
    {
        if (j == aidx || j == bidx)
            continue;
        if (isRight(points[aidx], points[bidx], points[j]))
        {
            s1.push_back(j);
        }
        else if (isRight(points[bidx], points[aidx], points[j]))
        {
            s2.push_back(j);
        }
    }
    // find hull ab, ba using the indices of the vectors in points
    // TODO: work out the indexing for this thingy
    list<int>::iterator apos = polygon.begin();
    list<int>::iterator bpos = next(polygon.begin(), 1);
    list<int>::iterator apos2 = next(bpos, 1);
    // cout << "in findHull\n";
    // cout << "s1: ";
    // for(int i:s1){
    //     cout << i<<" ";
    // }
    // cout << "\n";
    // cout << "s2: ";
    // for(int j:s2){
    //     cout << j << " ";
    // }
    // cout << "\n";
    findHull(points, s1, polygon, apos, bpos);
    findHull(points, s2, polygon, bpos, apos2);
    return polygon;
}
void part1()
{

    //  c) generate 60 random points in the unit square (save them in a data structure of your chosing) and save them in points.txt in the same format as 3.4 (each line has the x and y coordinate of a point separated by 2 spaces, use at least 20 digits precision)
    generatePoints(60);
    vector<DoublePoint> points = splitpt2("points.txt");
    //  d) apply the QuickHull algorithm to find a convex hull that encloses all 60 points and the vertices are points from the set of 60 points you generated
    // vector<DoublePoint> hull;
    // vector<int> hull;
    list<int> hull = quickHull(points);
    // quickHull()
    // quickHull(points);
    //  e) create a ppm of size 400x400 named quickhull.ppm in which you display a circle of radius 3 for each scaled point you generated, as well as connect the vertices you obtained for the convex hull (so in the image I should see 60 points and a convex hull with vertices circles of radus 3 that are conected forming a convex hull)
    Grid gr(400);
    for (int i = 0; i < (int)points.size(); i++)
    {
        highLightPoint(gr, points[i], 3, 0);
        // cout << points[i].toString() << "\n";
    }

    highLightPoint(gr, points[*hull.begin()], 4, 2);
    // cout << "\n";
    //  for(auto it=hull.begin();it!=hull.end();it++){
    //      //cout << points[*it].toString() << "\n";
    //  }

    for (auto j = next(hull.begin(), 1); j != hull.end(); j++)
    {
        int idx1 = *j;
        int idx2 = *std::prev(j, 1);
        highLightPoint(gr, points[idx1], 4, 2);
        IntegerPoint p1 = points[idx1].scaleToInt(gr.getSize());
        IntegerPoint p2 = points[idx2].scaleToInt(gr.getSize());
        gr.drawLine(p1, p2);
    }
    gr.toPPM("quickhull.ppm");
}
/**
 * @brief the comparator function for the Graham scan
 *
 */

class GrahamScanComparer
{

private:
    DoublePoint basePoint;

public:
    GrahamScanComparer(DoublePoint &dp)
    {
        basePoint = dp;
    }
    static double calcCos(const DoublePoint &u, const DoublePoint &v)
    {
        double magu = sqrt(u.getX() * u.getX() + u.getY() * u.getY());
        double magv = sqrt(v.getX() * v.getX() + v.getY() * v.getY());
        double dot = u.getX() * v.getX() + u.getY() * v.getY();
        return dot / (magu * magv);
    }
    double calcCos(const DoublePoint &u)
    {
        DoublePoint xaxis(1, 0);
        DoublePoint lxp(u.getX() - basePoint.getX(), u.getY() - basePoint.getY());
        return calcCos(lxp, xaxis);
    }
    bool operator()(const DoublePoint &lx, const DoublePoint &rx) const
    {
        // calculate the cosine of each

        DoublePoint xaxis(1, 0);
        DoublePoint lxp(lx.getX() - basePoint.getX(), lx.getY() - basePoint.getY());
        DoublePoint rxp(rx.getX() - basePoint.getX(), rx.getY() - basePoint.getY());
        double lcos = calcCos(lxp, xaxis);
        double rcos = calcCos(rxp, xaxis);
        if (abs(lxp.getX()) < 1e-15 && abs(lxp.getY()) < 1e-15)
        {
            lcos = 1.0;
        }
        if (abs(rxp.getX()) < 1e-15 && abs(rxp.getY()) < 1e-15)
        {
            rcos = 1.0;
        }
        // smaller cosine = bigger angle
        if (lcos > rcos)
        {
            return true;
        }
        else if (abs(lcos - rcos) < 1e-15)
        {
            return Line::calcDist(lx, basePoint) > Line::calcDist(rx, basePoint);
        }
        return false;
    }
};
/**
 * @brief Runs the graham scan algorithm
 *
 */
stack<DoublePoint> grahamScan(vector<DoublePoint> &points)
{
    // let points be the list of points
    // let stack = empty_stack()

    // find the lowest y-coordinate and leftmost point, called P0
    DoublePoint p0 = points[0];
    size_t idx = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
        if (points[i].getY() < p0.getY())
        {
            p0 = points[i];
            idx = i;
        }
        else if (abs(points[i].getY() - p0.getY()) < 1e-15)
        {
            if (points[i].getX() < p0.getX())
            {
                p0 = points[i];
                idx = i;
            }
        }
    }
    vector<DoublePoint> sortedPoints;
    for(size_t i=0; i < points.size();i++){
        if(i!=idx){
            sortedPoints.push_back(points[i]);
        }
    }
    // stk.push_back(p0);
    //  sort points by polar angle with P0, if several points have the same polar angle then only keep the farthest
    std::sort(sortedPoints.begin(), sortedPoints.end(), GrahamScanComparer(p0));
    //cout << "P0: " << points[0].toString() << "\n";
    // for point in points:
    //     # pop the last point from the stack if we turn clockwise to reach this point
    //     while count stack > 1 and notALeftTurn(A,B, point):
    //         pop stack
    //     push point to stack
    // end
    // cout << "Processing order: ";
    stack<DoublePoint> stk;
    stk.push(p0);
    stk.push(sortedPoints[0]);
    GrahamScanComparer comp(p0);
    for (size_t j = 0; j < sortedPoints.size(); j++)
    {
        // DoublePoint &top = stk.back();
        // DoublePoint &nextToTop = *std::prev(stk.end(),2);
        double jcos = comp.calcCos(sortedPoints[j]);
        double prevCos = comp.calcCos(sortedPoints[j - 1]);
        if (abs(jcos - prevCos) < 1e-15)
        {
            continue;
        }
        DoublePoint B = stk.top();
        stk.pop();
        DoublePoint A = stk.top();
        while (stk.size() > 1 && isRight(A, B, sortedPoints[j]))
        {
            B = stk.top();
            stk.pop();
            A = stk.top();
        }
        stk.push(B);
        stk.push(sortedPoints[j]);
        // cout << points[j].toString()<<endl;
    }
    stk.push(p0);
    // cout << "Done";
    return stk;
}
void part2()
{
    //generatePoints(60);
    vector<DoublePoint> points = splitpt2("points.txt");
    stack<DoublePoint> hull = grahamScan(points);

    Grid gr(400);
    for (int i = 0; i < (int)points.size(); i++)
    {
        highLightPoint(gr, points[i], 3, 0);
        // cout << points[i].toString() << "\n";
    }

    highLightPoint(gr, hull.top(), 4, 3);
    // B = top
    // A = next
    // while len of stack > 0
    // draw line between A & B
    // push A back on stack
    DoublePoint p2 = hull.top();
    hull.pop();
    DoublePoint p1 = hull.top();
    // cout << "line 1211\n";
    while (hull.size() > 1)
    {
        highLightPoint(gr, p1, 4, 3);
        IntegerPoint intp1 = p1.scaleToInt(gr.getSize());
        IntegerPoint intp2 = p2.scaleToInt(gr.getSize());
        gr.drawLine(intp1, intp2);
        p2 = p1;
        hull.pop();
        p1 = hull.top();
    }
    IntegerPoint intp1 = p1.scaleToInt(gr.getSize());
    IntegerPoint intp2 = p2.scaleToInt(gr.getSize());
    gr.drawLine(intp1, intp2);

    gr.toPPM("grahamscan.ppm");
}

int main()
{
    // part1();
    part2();
}