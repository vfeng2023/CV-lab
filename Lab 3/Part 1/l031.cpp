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

using namespace std;
const int PRECISION = 23;
// typedef std::numeric_limits< double > dbl;

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

list<DoublePoint> split(string &filename)
{
    fstream myFile;
    myFile.open(filename);
    list<DoublePoint> toRet;
    string line;
    while(getline(myFile,line)){
        stringstream ss(line);
        vector<string> coords;
        string buffer;
        while(ss>>buffer){
            if(!buffer.empty())
                coords.push_back(buffer);
        }
        double x = stod(coords[0]);
        double y = stod(coords[1]);
        DoublePoint dp(x,y);
        toRet.push_back(dp);
    }
    return toRet;

}
/*
    calculates the things
*/
void highLightPoint(Grid &gr, DoublePoint &p, int rad, int color = 2)
{
    Circle c(p.getX(), p.getY(), rad*1.0/gr.getSize());
    c.drawCircle(gr, color,rad);
    c.drawCircle(gr,color,rad-1);
}

void part0(){
    cout << "Do you want to generate random points(Yes/No, y for Yes, n for No): ";
    char ans;
    cin >> ans;
    srand(time(NULL));
    if(ans=='y' || ans=='Y'){
        ofstream pointsFile;
        pointsFile.open("points.txt");
        pointsFile.precision(PRECISION);
        for(int i=0; i < 60;i++){
            double x = rand()*1.0/RAND_MAX;
            double y = rand()*1.0/RAND_MAX;
            pointsFile <<fixed<< x << "  "<< y << endl;
        }
        pointsFile.close();
    }
}
void part1(){
    //read file to a list
    string filename = "points.txt";
    list<DoublePoint> myPoints = split(filename);
    double minDist = DBL_MAX;
    DoublePoint mindp1;
    DoublePoint mindp2;
    //apply brute force
    //for it = l.begin, it!=l.end, it++
        //for it2 = it; ...
        //calc distance
        //if less than mindistance, save
    for(list<DoublePoint>::iterator first=myPoints.begin();first!=myPoints.end();first++){
        for(list<DoublePoint>::iterator second=next(first); second!=myPoints.end();second++){
            DoublePoint &dp1 = *first;
            DoublePoint &dp2 = *second;
            double dist = Line::calcDist(dp1,dp2);
            if(dist < minDist){
                mindp1 = dp1;
                mindp2 = dp2;
                minDist = dist;
            }
        }
    }
    //draw points in ppm
    Grid gr(800);
    for(list<DoublePoint>::iterator it=myPoints.begin();it!=myPoints.end();it++){
        highLightPoint(gr,*it,3,0);
    }
    highLightPoint(gr,mindp1,3,2);
    highLightPoint(gr,mindp2,3,2);
    cout.precision(PRECISION);
    cout << "Point 1: "<<mindp1.toString(true)<<endl;
    cout << "Point 2: "<<mindp2.toString(true)<<endl;
    cout << "Minimum Distance: "<<minDist<<endl;
    gr.toPPM("points.ppm");

    //die

}
int main(){
    part0();
    part1();
}