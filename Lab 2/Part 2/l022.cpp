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

using namespace std;

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
    IntegerPoint scaleToInt(int size)
    {
        IntegerPoint intpoint((int)(round(size * x)), (int)(round(size * y)));
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
    Grid(const Grid &gr){
        grid = new int*[gr.size];
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
        if(abs(slope) < 1e-8){
            DoublePoint start(0,d1.getY());
            DoublePoint end(1,d1.getY());
            IntegerPoint i1 = start.scaleToInt(size);
            IntegerPoint i2 = end.scaleToInt(size);
            drawLine(i1, i2);
            cout << "\nslope: horizontal line"<<endl;
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
        //find intersection with all sides of the bounding box: x = 0, x=1, y=0, y=1
        //take the points which have a positive coordinate
        double yx0 = slope * (-1.0 * d1.getX()) + d1.getY();
        double yx1 = slope * (1.0 - d1.getX()) + d1.getY();
        double xy0 = -d1.getY()/slope + d1.getX();
        double xy1 = (1.0 - d1.getY())/slope + d1.getX();
        bool chosestart = false;

        DoublePoint start(0.0, yx0);
        DoublePoint end(1, yx1);
        if(yx0 >= 0 && yx0 <=1){
            if(!chosestart){
                start.setX(0);
                start.setY(yx0);
                chosestart = true;
            }else{
                end.setX(0);
                end.setY(yx0);
            }
            
        }
        if(yx1 >=0 && yx1 <=1){
            if(!chosestart){
                start.setX(1);
                start.setY(yx1);
                chosestart=true;
            }else{
                end.setX(1);
                end.setY(yx1);
            }
        }
        if(xy0 >=0 && xy0 <=1){
            if(!chosestart){
                
                start.setX(xy0);
                start.setY(0);
                chosestart =true;
            }else{
                end.setX(xy0);
                end.setY(0);
               
            }
        }
        if(xy1 >=0 && xy1 <=1){
            if(!chosestart){
                start.setX(xy1);
                start.setY(1.0);
                chosestart=true;
            }else{
                end.setX(xy1);
                end.setY(1.0);
            }
        }
        cout << "\n"<<start.toString()<<endl;
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
    Line(const Line &line){
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
    void draw(Grid &gr)
    {
        gr.drawLine(p1, p2);
    }
};
class Circle
{
private:
    DoublePoint *center;
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
        center = new DoublePoint(x, y);
        radius = rad;
    }
    ~Circle()
    {
        delete center;
        // cout << "Circle deleted"<<endl;
    }
    DoublePoint &getCenter()
    {
        return *center;
    }
    double getRadius()
    {
        return radius;
    }
    string toString()
    {
        return "Center: " + center->toString() + "; radius=" + to_string(radius);
    }
    /*
        Draws the circle on a grid;
    */
    void drawCircle(Grid &gr,int color = 0,int customRad=-1)
    {
        int size = gr.getSize();
        IntegerPoint intCenter = center->scaleToInt(size);

        // midpoint circle algorithm

        int r = (int)(round(radius * size));
        if(customRad > 0){
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
/*
    Part 1 generates three random points, then while the fourth point in the triangle, generates new points
    Uses Point-in-Polygon algorithm. Starts with ray, and looks at pairs of points to see if an horizontal ray emanating form the point intersects these segments
*/
/*
    Finds the cross product of two points (aka determinant)
*/
// double crossProduct(DoublePoint p2, DoublePoint v1, DoublePoint v2){
//     double crossProd = (p1.getX()-v1.getX())*(p2.getY() - p1.getY())* p2.getX();
//     cout << "p1: "<<p1.toString()<< " ";
//     cout << " p2 "<<p2.toString() <<" " << crossProd << endl;
//     return crossProd;
// }
/*
    Uses barycentric coords to check. u and v >= 0
    // Compute vectors
    v0 = C - A
    v1 = B - A
    v2 = P - A

    // Compute dot products
    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1)
    dot02 = dot(v0, v2)
    dot11 = dot(v1, v1)
    dot12 = dot(v1, v2)

    // Compute barycentric coordinates
    invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom

    // Check if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v < 1)
    https://blackpawn.com/texts/pointinpoly/default.html
*/
double dotProduct(DoublePoint &v1, DoublePoint &v2)
{
    return v1.getX() * v2.getX() + v1.getY() * v2.getY();
}
bool sameSideChecker(DoublePoint p, DoublePoint a, DoublePoint b, DoublePoint c)
{
    DoublePoint v0(c.getX() - a.getX(), c.getY() - a.getY());
    DoublePoint v1(b.getX() - a.getX(), b.getY() - a.getY());
    DoublePoint v2(p.getX() - a.getX(), p.getY() - a.getY());
    double dot00, dot01, dot02, dot11, dot12;
    dot00 = dotProduct(v0, v0);
    dot01 = dotProduct(v0, v1);
    dot02 = dotProduct(v0, v2);
    dot11 = dotProduct(v1, v1);
    dot12 = dotProduct(v1, v2);
    double denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * denom;
    double v = (dot00 * dot12 - dot01 * dot02) * denom;
    return u >= 0 && v >= 0 && (u + v) < 1;
}
/*
    - Uses the determinant algorithm to check containment
    for each point, determine the sign of c1, c2, c3 are the same
    - return true if the triangle contains the point at i
*/
bool checkIntersect(const int &i, DoublePoint (&dplist)[4])
{
    vector<DoublePoint> tri = vector<DoublePoint>();
    for (int j = 0; j < 4; j++)
    {
        if (j != i)
            tri.push_back(dplist[j]);
    }
    return sameSideChecker(dplist[i], tri[0], tri[1], tri[2]);
    // math from https://blackpawn.com/texts/pointinpoly/default.html
}
void part1()
{
    srand(time(NULL));
    ofstream myLog;
    myLog.open("log.txt");
    DoublePoint dplist[4];
    for (int i = 0; i < 4; i++)
    {
        double x = rand() * 1.0 / RAND_MAX;
        double y = (double)(rand()) / RAND_MAX;
        dplist[i] = DoublePoint(x, y);
    }
    for (int i = 0; i < 3; i++)
    {
        std::cout << dplist[i].toString();
        myLog << dplist[i].toString();
        if (i < 2)
        {
            std::cout << " , ";
            myLog << " , ";
        }
    }
    // dplist[0].getX() = 0;
    // dplist[0].getY() = .5;
    // dplist[1].getX()=.5;
    // dplist[1].getY()=1;
    // dplist[2].getX() = 1;
    // dplist[2].getY() = .4;
    // DoublePoint p4(rand()/RAND_MAX,rand()/RAND_MAX);
    bool inTriangle = true;
    /*
        PIP algorithm - count the number of times horizontal ray of p4 crosses the sides.
        If Even number of times, then point is outside triangle
    */
    // while point overlaps triangle
    // 1. find intersection of x and the segment
    // 2. if POI is greater than x, within bounds of segment, then increment count
    // if count is even, break
    cout << "\n";
    myLog << "\n";
    while (inTriangle)
    {
        // int crosscount = 0;
        dplist[3].setX((double)(rand()) / RAND_MAX);
        dplist[3].setY((double)(rand()) / RAND_MAX);
        cout << "testing point " << dplist[3].toString() << endl;
        myLog << "testing point " << dplist[3].toString() << endl;
        bool intersect = false;
        for (int i = 0; i < 4; i++)
        {
            if (checkIntersect(i, dplist))
            {
                intersect = true;
                break;
            }
        }
        if (!intersect)
        {
            break;
        }
        // for(int i=0; i < 4; i++){
        //     std::cout << dplist[i].toString()<<",";
        // }
        // std::cout << "\n";
        // string s;

        // cin >> s;
    }
    cout << "\n"
         << "Done. Points: " << endl;
    myLog << "\n"
          << "Done. Points: " << endl;
    // cout.precision(15);
    // myLog.precision(15);
    ofstream myPoints;
    myPoints.open("points.txt");
    //myPoints.open("points_mrjurg.txt");
    myPoints.precision(17);
    for (int i = 0; i < 4; i++)
    {
        std::cout << "(" << dplist[i].getX() << "," << dplist[i].getY() << ")";
        myLog << "(" << dplist[i].getX() << "," << dplist[i].getY() << ")";
        myPoints << fixed << "(" << dplist[i].getX() << "," << dplist[i].getY() << ")";
        if (i < 3)
        {
            std::cout << " , ";
            myLog << " , ";
            myPoints << " , ";
        }
    }
    cout << '\n';
    myLog << "\n";
    myPoints.close();
    myLog.close();
}
vector<DoublePoint> split(string &myLine)
{
    // int start = 0;
    // vector<string> toRet;
    // size_t delimLen = delim.size();

    // while(myLine.find(delim)!=)
    string toNum = "";
    // get the list of doubles
    // cout << myLine<<endl;
    vector<double> ptvals;

    for (char ch : myLine)
    {
        if (ch == '.' || (ch - '0' < 10 && ch - '0' >= 0))
        {
            toNum += ch;
        }
        else
        {
            if (toNum.length() > 0)
            {
                // std::cout << toNum<<","<<endl;
                ptvals.push_back(stod(toNum));
            }

            // ptvals.push_back(stod(toNum));
            toNum = "";
        }
    }
    if (toNum.length() > 0)
    {
        // ptvals.push_back(stod(toNum));
    }
    vector<DoublePoint> myPoints;
    for (unsigned long int i = 0; i < ptvals.size(); i += 2)
    {
        DoublePoint pt(ptvals[i], ptvals[i + 1]);
        myPoints.push_back(pt);
    }
    return myPoints;
}
/*
    calculates the things
*/
void highLightPoint(Grid &gr, DoublePoint &p, int rad, int color = 2)
{
    Circle c(p.getX(), p.getY(), rad*1.0/gr.getSize());
    c.drawCircle(gr, color,rad);
}
/*
    Returns array containing the vertices of the square given four points random points.
    returned vector is of length 4
*/
pair<vector<DoublePoint>,vector<Line>> calculateLines(DoublePoint &a, DoublePoint &b, DoublePoint &c, DoublePoint &d, int dir = 1)
{
    // 1. pick two points that will end up on opposite sides, then connect them
    Line lineac(a, c);
    // 2. select other point and draw perpendicular to 1
    Line perp = lineac.findPerpFromPoint(b);
    // Below checks visually for peperndicularity
    // Grid gr(800);
    // lineac.draw(gr);
    // gr.drawLine(a,c);
    // perp.draw(gr);
    // highLightPoint(gr,a,3);
    // highLightPoint(gr,b,3);
    // highLightPoint(gr,c,3);
    // highLightPoint(gr,d,3);
    // gr.toPPM();
    // 3. Pick E on line from step 2 such that BE = AC in length
    // find unit vector taking step in direction of x with slope of perp
    // return E such that E ix on the line be
    cout << "Perpendicular slope: " << perp.getSlope() << endl;
    cout << "Line ac slope: " << lineac.getSlope() << endl;
    double magn = sqrt(1 + (perp.getSlope() * perp.getSlope()));
    double ex = b.getX() + dir / magn * lineac.calcDist();
    DoublePoint e(ex, perp.evaluateAt(ex));
    Line lineeb = Line(e, b);
    cout << "distance: " << lineeb.calcDist() << " " << lineac.calcDist();
    // 4. Connect last point with E.  This will be one line that forms the square
    Line linede(d, e);

    // 5. Draw perpendiculars from the points in step 1 on the line from step 4
    // a and c are opposite each other, b and d are opposite each other
    Line sideade = linede.findPerpFromPoint(a);
    Line sidecde = linede.findPerpFromPoint(c);
    // 6. From the point in step 2 - B - draw a perpendicular on any of the two lines from step 5
    Line sidebac = sideade.findPerpFromPoint(b);
    // Find the intersections
    // First, b's intersections with a and d
    DoublePoint p1 = sidebac.findIntersect(sideade);
    DoublePoint p2 = sidebac.findIntersect(sidecde);
    DoublePoint p3 = linede.findIntersect(sidecde);
    DoublePoint p4 = linede.findIntersect(sideade);
    //double sidelength = Line::calcDist(p1, p2);
    vector<DoublePoint> dparr = {p1, p2, p3, p4};
    vector<Line> sidearray = {linede,sideade,sidecde,sidebac};
    return make_pair(dparr,sidearray);
}

/*
    Part 2 draws the triangle in the ppm file
*/
void part2()
{
    // Read the file
    fstream myFile("points.txt", ios::in);

    string line;
    getline(myFile, line); // get the line
    myFile.close();
    // split contents of line by " , "

    vector<DoublePoint> myPoints = split(line);
    ofstream outputFile;
    outputFile.open("output.txt");
    // for(DoublePoint dp:myPoints){

    //     outputFile << ",";
    // }
    // print to file and draw points
    Grid gr(800);

    outputFile.precision(dbl::max_digits10);
    for (int i = 0; i < 4; i++)
    {
        // outputFile << "(" << myPoints[i].getX() <<","<<myPoints[i].getY()<<")";
        outputFile << myPoints[i].toString(true);
        highLightPoint(gr, myPoints[i], 3, 2);
        if (i < 3)
        {
            outputFile << " , ";
        }
    }
    outputFile << "\n";
    // outputFile.close();
    // for i=1:3,
    // find 0, i, remaining two numbers,
    // calc with dir = 1
    // calc with dir = -1
    vector<DoublePoint> minPointArr;
    vector<Line> minSideArr;
    /// vector<DoublePoint> *currArr = NULL;
    double minArea = DBL_MAX;
    for (int i = 1; i < 4; i++)
    {
        // int a = 0;
        vector<int> vals = {0, i};
        for (int j = 1; j < 4; j++)
        {
            if (j != i)
            {
                vals.push_back(j);
            }
        }

        int a = vals[0];
        int c = vals[1];
        int b = vals[2];
        int d = vals[3];
        cout <<"\n"<< a <<" "<< b <<" "<<c<<" "<< d <<" "<< endl;
        vector<int> possDir = {-1,1};
        for (int dir:possDir)
        {
            pair<vector<DoublePoint>, vector<Line>> retVal= calculateLines(myPoints[a], myPoints[b], myPoints[c], myPoints[d], dir);
            vector<DoublePoint> currArr = get<0>(retVal);
            vector<Line> sideArr = get<1>(retVal);
            double side = (Line::calcDist(currArr[0], currArr[1]));
            double area = side * side;
            for (int i = 0; i < 4; i++)
            {
                outputFile << currArr[i].toString(true);
                if (i < 3)
                {
                    outputFile << " , ";
                }
            }
            outputFile << " Area = " << area << endl;
            if (area < minArea)
            {
                //delete minPointArr;
                minPointArr.clear();
                for (DoublePoint dp : currArr)
                {
                    minPointArr.push_back(dp);
                }
                minSideArr.clear();
                for(Line l: sideArr){
                    minSideArr.push_back(l);
                }
                minArea = area;
            }
        }
        //
    }

    // calculate the points and write them to the file.
    // if min area is less than current min, save the current array
    //
    // ofstream
    cout <<"\nminArea="<<minArea<<endl;
    //cout << ""
    for (int i = 0; i < 4; i++)
    {
        // DoublePoint &p1 = minPointArr[i];
        // DoublePoint &p2 = minPointArr[(i + 1) % 4];
        minSideArr[i].draw(gr);
    }
    for (int i = 0; i < 4; i++)
    {
        highLightPoint(gr, minPointArr[i], 6, 0);
    }
    //delete minPointArr;
    // draw the points by drawing the lines
    outputFile.close();
    gr.toPPM("output.ppm");
    //cout << "Finished drawing";
}
int main()
{

    //part1();
    part2();
}