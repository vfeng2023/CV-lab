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
    double getX()
    {
        return x;
    }
    double getY()
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
            grid[p1.getY()][p1.getX()] = 0;
    }
    /*
        Choices in color - 0 = black, 1 = white, 2 = red, 3 = green, 4 = blue
    */
    void drawPixel(IntegerPoint &p1, int color)
    {
        if (p1.getX() < size && p1.getY() < size && p1.getX() >= 0 && p1.getY() >= 0)
            grid[p1.getY()][p1.getX()] = color;
    }
    void drawPixel(int x, int y, int color = 0)
    {
        if (x < size && y < size && x >= 0 && y >= 0)
            grid[y][x] = color;
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
    double static calcDist(DoublePoint &d1, DoublePoint &d2)
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

list<DoublePoint> split(string &filename)
{
    fstream myFile;
    myFile.open(filename);
    list<DoublePoint> toRet;
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
    return toRet;
}
vector<DoublePoint> splitpt2(string &filename)
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
/*
    calculates the things
*/
void highLightPoint(Grid &gr, DoublePoint &p, int rad, int color = 2)
{
    Circle c(p.getX(), p.getY(), rad * 1.0 / gr.getSize());
    c.drawCircle(gr, color, rad);
    c.drawCircle(gr, color, rad - 1);
}

void part0(int numPoints)
{
    // int numPoints = 60;
    // cout << "Do you want to generate random points(Yes/No, y for Yes, n for No): ";
    // char ans;
    // cin >> ans;
    srand(time(NULL));
    //if (ans == 'y' || ans == 'Y')
    //{
        ofstream pointsFile;
        pointsFile.open("points.txt");
        pointsFile.precision(PRECISION);

        // cout << "\nWhat size of data do you wish to generate? ";
        // int numPoints;
        //cin >> numPoints;
        for (int i = 0; i < numPoints; i++)
        {
            double x = rand() * 1.0 / RAND_MAX;
            double y = rand() * 1.0 / RAND_MAX;
            pointsFile << fixed << x << "  " << y << endl;
        }
        pointsFile.close();
    //}
}
string part1()
{
    // read file to a list
    string filename = "points.txt";
    list<DoublePoint> myPoints = split(filename);
    double minDist = DBL_MAX;
    DoublePoint mindp1;
    DoublePoint mindp2;
    // apply brute force
    // for it = l.begin, it!=l.end, it++
    // for it2 = it; ...
    // calc distance
    // if less than mindistance, save
    auto start = high_resolution_clock::now();
    for (list<DoublePoint>::iterator first = myPoints.begin(); first != myPoints.end(); first++)
    {
        for (list<DoublePoint>::iterator second = next(first); second != myPoints.end(); second++)
        {
            DoublePoint &dp1 = *first;
            DoublePoint &dp2 = *second;
            double dist = Line::calcDist(dp1, dp2);
            if (dist < minDist)
            {
                mindp1 = dp1;
                mindp2 = dp2;
                minDist = dist;
            }
        }
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    // draw points in ppm
    // Grid gr(800);
    // for (list<DoublePoint>::iterator it = myPoints.begin(); it != myPoints.end(); it++)
    // {
    //     highLightPoint(gr, *it, 3, 0);
    // }
    // highLightPoint(gr, mindp1, 3, 2);
    // highLightPoint(gr, mindp2, 3, 2);
    stringstream stream;
    stream.precision(PRECISION);
    stream << "Brute force method (Part 1)" << endl;
    stream << "Time taken: " << duration.count() << " microseconds" << endl;
    stream << "Point 1: " << mindp1.toString(true) << endl;
    stream << "Point 2: " << mindp2.toString(true) << endl;
    stream << "Minimum Distance: " << minDist << endl;
    // gr.toPPM("points.ppm");
    return stream.str();
    // die
}
class MergeReturnObj
{
private:
    DoublePoint dp1;
    DoublePoint dp2;
    double distance;

public:
    MergeReturnObj()
    {
        distance = 0;
    }
    MergeReturnObj(double d, DoublePoint p1, DoublePoint p2)
    {
        dp1 = p1;
        dp2 = p2;
        distance = d;
    }
    MergeReturnObj(const MergeReturnObj &mro)
    {
        dp1 = mro.dp1;
        dp2 = mro.dp2;
        distance = mro.distance;
    }
    double getDistance()
    {
        return distance;
    }
    DoublePoint getPoint1()
    {
        return DoublePoint(dp1);
    }
    DoublePoint getPoint2()
    {
        return DoublePoint(dp2);
    }
};
MergeReturnObj baseCaseMerge(vector<DoublePoint> &points, int start, int end)
{
    double minDist = DBL_MAX;
    DoublePoint mindp1;
    DoublePoint mindp2;
    for (int i = start; i <= end; i++)
    {
        for (int j = i + 1; j <= end; j++)
        {
            double dist = Line::calcDist(points[i], points[j]);
            if (dist < minDist)
            {
                minDist = dist;
                mindp1 = points[i];
                mindp2 = points[j];
            }
        }
    }
    return MergeReturnObj(minDist, mindp1, mindp2);
}
MergeReturnObj zipperMerge(vector<DoublePoint> &points, int start, int mid, int end)
{
    double minDist = DBL_MAX;
    int mindp1 = 0;
    int mindp2 = 0;
    // cout << "line 960";
    // cout << "start and end: "<<start << " "<< end<<endl;
    for (int i = start; i <= end; i++)
    {
        for (int j = i + 1; j <= end; j++)
        {
            // if (i != j)
            //{
            DoublePoint &p1 = points[i];
            DoublePoint &p2 = points[j];
            double dist = Line::calcDist(p1, p2);
            if (dist < minDist)
            {
                mindp1 = i;
                mindp2 = j;
                minDist = dist;
            }
            //}
        }
    }
    return MergeReturnObj(minDist, points[mindp1], points[mindp2]);
}
MergeReturnObj merge(vector<DoublePoint> &points, int start, int end)
{
    if (end - start + 1 <= 3)
    {
        return baseCaseMerge(points, start, end);
    }
    else
    {
        // static int count = 0;
        // int currCount = count;

        // cout << "merge function start and end: " << start <<" "<<end<< "level"<<count<<endl;
        int midpoint = (start + end) / 2;
        MergeReturnObj m1 = merge(points, start, midpoint - 1);
        MergeReturnObj m2 = merge(points, midpoint, end);
        double smallest = m2.getDistance();
        MergeReturnObj *smallestObj = &m2;
        if (m1.getDistance() < m2.getDistance())
        {
            smallest = m1.getDistance();
            smallestObj = &m1;
        }
        // create strip
        double lbound = points[midpoint].getX() - smallest;
        double rbound = points[midpoint].getX() + smallest;
        // brute force on the strip
        int l = midpoint;
        // while(l > start && abs(points[l].getX()-points[midpoint].getX()) <= smallest){
        while (l > start && points[l].getX() > lbound)
        {
            l--;
        }
        int r = midpoint;
        // while( r < end && abs(points[r].getX()-points[midpoint].getX()) <= smallest){
        while (r < end && points[r].getX() < rbound)
        {
            r++;
        }
        // cout << "l="<<l<<" r="<<r<< " count "<< currCount<<endl;
        MergeReturnObj strip = zipperMerge(points, l, midpoint, r);
        // count += 1;
        if (strip.getDistance() < smallest)
        {
            return strip;
        }
        else
        {
            return *smallestObj;
        }
        // return points in the strip that have the smallest distance
    }
}
string part2()
{
    string filename = "points.txt";
    vector<DoublePoint> myPoints = splitpt2(filename);
    auto start = high_resolution_clock::now();
    std::sort(myPoints.begin(), myPoints.end(), DoublePoint::comparePoint);
    // for(DoublePoint i:myPoints)
    //     cout << "Point"<<i << endl;
    MergeReturnObj mergeRet = merge(myPoints, 0, myPoints.size() - 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    stringstream stream;
    stream.precision(PRECISION);
    stream << "\n\nMerge sort method (Part 2) " << endl;
    stream << "Time taken: " << duration.count() << " microseconds " << endl;
    // double distance = mergeRet.getDistance();
    DoublePoint dp1 = mergeRet.getPoint1();
    DoublePoint dp2 = mergeRet.getPoint2();

    stream << "Point 1: " << dp1.toString(true) << endl;
    stream << "Point 2: " << dp2.toString(true) << endl;
    stream << "Minimum Distance: " << mergeRet.getDistance() << endl;
    return stream.str();
}

// PART 3
MergeReturnObj zipperMerge3(vector<DoublePoint> &vectorStrip)
{
    double minDist = DBL_MAX;
    DoublePoint mindp1;
    DoublePoint mindp2;
    // cout << "line 960";
    // cout << "start and end: "<<start << " "<< end<<endl;
    // vector<DoublePoint> vectorStrip;
    // for (int k = start; k <= end; k++)
    // {
    //     vectorStrip.push_back(points[k]);
    // }
    std::sort(vectorStrip.begin(), vectorStrip.end(), DoublePoint::comparePointY);
    for (size_t i = 0; i < vectorStrip.size(); i++)
    {
        size_t upper = min(i + 15, vectorStrip.size());
        // cout << (upper==i+15) << endl;
        for (size_t j = i + 1; j < upper; j++)
        {
            DoublePoint p1 = vectorStrip[i];
            DoublePoint p2 = vectorStrip[j];
            double dist = Line::calcDist(p1, p2);
            if (dist < minDist)
            {
                mindp1 = p1;
                mindp2 = p2;
                minDist = dist;
            }
        }
    }
    // cout << "1144\n";
    // cout << "vectorStrip.size()=" << vectorStrip.size();
    // cout << "\nmindp1 = "<<mindp1<<endl;
    // cout << "\nmindp2 = "<< mindp2 << endl;
    return MergeReturnObj(minDist, mindp1, mindp2);
}
MergeReturnObj merge3(vector<DoublePoint> &points, int start, int end)
{
    if (end - start + 1 <= 3)
    {
        return baseCaseMerge(points, start, end);
    }
    else
    {
        // static int count = 0;
        // int currCount = count;

        // cout << "merge function start and end: " << start <<" "<<end<< "level"<<count<<endl;
        int midpoint = (start + end) / 2;
        MergeReturnObj m1 = merge3(points, start, midpoint - 1);
        MergeReturnObj m2 = merge3(points, midpoint, end);
        double smallest = m2.getDistance();
        MergeReturnObj smallestObj(m2);
        if (m1.getDistance() < m2.getDistance())
        {
            smallest = m1.getDistance();
            smallestObj = m1;
        }
        // create strip
        // cout << "Line 1169\n";
        vector<DoublePoint> vectorStrip;
        for (int k = start; k <= end; k++)
        {
            if (abs(points[k].getX() - points[midpoint].getX()) < smallest)
                vectorStrip.push_back(points[k]);
        }
        // cout << "1176\n";
        // MergeReturnObj strip = zipperMerge3(points, l, midpoint, r);
        MergeReturnObj strip = zipperMerge3(vectorStrip);
        // count += 1;
        if (strip.getDistance() < smallest)
        {
            return strip;
        }
        else
        {
            return smallestObj;
        }
        // return points in the strip that have the smallest distance
    }
}
double part3()
{
    string filename = "points.txt";
    vector<DoublePoint> myPoints = splitpt2(filename);
    // cout << "Line 1193\n";
    auto start = high_resolution_clock::now();
    std::sort(myPoints.begin(), myPoints.end(), DoublePoint::comparePoint);
    // for(DoublePoint i:myPoints)
    //     cout << "Point"<<i << endl;
    MergeReturnObj mergeRet = merge3(myPoints, 0, myPoints.size() - 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    stringstream stream;
    stream.precision(PRECISION);
    stream << "\n\nMerge sort method (Part 3) " << endl;
    stream << "Time taken: " << duration.count() << " microseconds " << endl;
    // double distance = mergeRet.getDistance();
    DoublePoint dp1 = mergeRet.getPoint1();
    DoublePoint dp2 = mergeRet.getPoint2();

    stream << "Point 1: " << dp1.toString(true) << endl;
    stream << "Point 2: " << dp2.toString(true) << endl;
    stream << "Minimum Distance: " << mergeRet.getDistance() << endl;
    chrono::duration<double> sec = stop-start;
    return sec.count();
}

// PART 4
class Box
{
private:
    unsigned long long x;
    unsigned long long y;

public:
    Box(unsigned long long int xp, unsigned long long int yp)
    {
        x = xp;
        y = yp;
    }
    Box()
    {
        x = 0;
        y = 0;
    }
    Box(const Box &b)
    {
        x = b.x;
        y = b.y;
    }
    Box(Box &&b)
    {
        x = b.x;
        y = b.y;
    }
    unsigned long long getX() const
    {
        return x;
    }
    unsigned long long getY() const
    {
        return y;
    }
    // string to_string(){
    //     return "(" + to_string(x) + ","+to_string(y)+")";
    // }
    bool operator==(const Box &b) const
    {
        return x == b.x && y == b.y;
    }
};
class BoxHash
{
public:
    size_t operator()(const Box &b) const
    {
        // unsigned long long x = b.getX();
        // unsigned long long y = b.getY();
        return hash<unsigned long long>()(b.getX()) ^ (hash<unsigned long long>()(b.getY()) << 1);
    }
};
void shuffle(vector<DoublePoint> &points)
{
    int size = (int)points.size();
    srand(time(NULL));
    for (int i = 0; i < size - 1; i++)
    {
        int swapidx = rand() % (size - i) + i;
        std::iter_swap(points.begin() + i, points.begin() + swapidx);
        // cout << points[i].toString()<<endl;
    }
}
Box getBox(DoublePoint &dp, double minDist)
{
    double delta = minDist / 2;
    double eps = 1e-9;
    unsigned long long bx = (unsigned long long)(dp.getX() / delta+eps);
    unsigned long long by = (unsigned long long)(dp.getY() / delta+eps);
    return Box(bx, by);
}
ostream &operator<<(ostream &stream, const Box &b)
{
    stream << "[" << b.getX() << "," << b.getY() << "]";
    return stream;
}
MergeReturnObj randomizedAlgorithm(vector<DoublePoint> &points)
{
    double mindist = Line::calcDist(points[0], points[1]);
    // cout << "min dist: " << mindist << endl;
    unordered_map<Box, size_t, BoxHash> pointsdictionary = unordered_map<Box, size_t, BoxHash>(); // maps Box to index
    // make dictionary
    Box b0 = getBox(points[0], mindist);
    Box b1 = getBox(points[1], mindist);
    pointsdictionary[b0] = 0;
    pointsdictionary[b1] = 1;
    // for i=1..n:
    size_t idx1 = 0;
    size_t idx2 = 1;
    for (size_t i = 2; i < points.size(); i++)
    {
        // cout << pointsdictionary<< endl;

        // determine subsquare containing p1
        Box myBox = getBox(points[i], mindist);
        unsigned long long startRow = myBox.getX() - 2;
        if (myBox.getX() < 2)
        {
            startRow = 0;
        }
        unsigned long long startCol = myBox.getY() - 2;
        if (myBox.getY() < 2)
        {
            startCol = 0;
        }
        // looking up the close subsquares
        double currMin = mindist;
        int count = 0;
        // double nextMin = mindist;
        // look up the 25 subsquares close to pi
        for (unsigned long long r = startRow; r <= myBox.getX() + 2; r++)
        {
            for (unsigned long long c = startCol; c <= myBox.getY() + 2; c++)
            {
                Box key(r, c);
                count += 1;
                // compute the distance form points[i] to any of these subsquares
                if (pointsdictionary.find(key) != pointsdictionary.end())
                {
                    size_t p = pointsdictionary[key];
                    double dist = Line::calcDist(points[i], points[p]);
                    // find and store the min
                    if (dist < mindist)
                    {
                        mindist = dist;
                        idx1 = i;
                        idx2 = p;
                    }
                }
            }
        }
        // look up the 25 close subsquares
        // if there is a point pj w/ delta less than mindist
        // remake dictionary
        // if (count > 25)
        //     cout << "count : " << count << endl;
        if (mindist < currMin)
        {
            // cout << "Remake dictionary\n";
            pointsdictionary.clear();
            for (size_t k = 0; k <= i; k++)
            {
                Box newBox = getBox(points[k], mindist);
                pointsdictionary[newBox] = k;
            }
        }
        else
        {
            Box keyBox = getBox(points[i], mindist);
            pointsdictionary[keyBox] = i;
        }
        // add subsquare into dictionary
        //  cout << "Iteration "<< i<<endl;
        //  for(auto it=pointsdictionary.begin();it!=pointsdictionary.end();it++){
        //      //cout << "("<< it->first<
        //      cout << it->first << "=>" << points[it->second].toString()<<endl;
        //  }
        //  cout << "\n\n";
    }
    return MergeReturnObj(mindist, points[idx1], points[idx2]);
}
double part4()
{
    // read points (vector DoublePoint)
    string filename = "points.txt";
    vector<DoublePoint> points = splitpt2(filename);
    auto start = high_resolution_clock::now();
    // shuffle points
    shuffle(points);
    // apply algorithm
    MergeReturnObj mergeRet = randomizedAlgorithm(points);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    stringstream stream;
    stream.precision(PRECISION);
    stream << "\n\nRandomized Algorithm (Part 4) " << endl;
    stream << "Time taken: " << duration.count() << " microseconds " << endl;
    // double distance = mergeRet.getDistance();
    DoublePoint dp1 = mergeRet.getPoint1();
    DoublePoint dp2 = mergeRet.getPoint2();

    stream << "Point 1: " << dp1.toString(true) << endl;
    stream << "Point 2: " << dp2.toString(true) << endl;
    stream << "Minimum Distance: " << mergeRet.getDistance() << endl;
    chrono::duration<double> sec = stop-start;
    return sec.count();
}
int main()
{   
    
    
    ofstream results;
    results.open("results.txt");
    int npoints;
    cout << "Number of points: ";
    cin >> npoints;
    double p3avg = 0;
    double p4avg = 0;
    int numruns = 5;
    for(int i=0; i < numruns; i++){
        part0(npoints);
        cout << "Generated\n";
        double p3 = part3();
        double p4 = part4();
        p3avg += p3;
        p4avg += p4;
    }
    p3avg /= numruns;
    p4avg /=numruns;
    cout.precision(PRECISION);
    results.precision(PRECISION);
    cout << "Part 3 average runtime (Seconds): "<<p3avg << "\n";
    results << "Part 3 average runtime (Seconds): "<<p3avg << "\n";
    
    cout << "Part 4 average runtime(seconds): " << p4avg << endl;
    results << "Part 4 average runtime(seconds): " << p4avg << endl;
    
    results.close();
}