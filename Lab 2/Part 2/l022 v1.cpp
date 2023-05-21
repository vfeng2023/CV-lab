#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <vector>
#include <cmath>
//#include <iomanip>

using namespace std;
class IntegerPoint
{

public:
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
    string toString(){
        return "("+to_string(x)+","+to_string(y)+")";
    }
};
class DoublePoint
{
public:
    double x;
    double y;

public:
    DoublePoint(){
        x = 0;
        y = 0;
    }
    DoublePoint(double xp, double yp)
    {
        x = xp;
        y = yp;
    }
    string toString()
    {
        return ("(" + to_string(x)) + "," + to_string(y) + ")";
    }
    IntegerPoint scaleToInt(int size)
    {
        IntegerPoint intpoint((int)(round(size*x)), (int)(round(size*y)));
        return intpoint;
    }
};
class Grid
{
private:
    int size;
    int** grid;
    

public:
    Grid(int mySize)
    {
        //grid = vector<vector<int>>(mySize, vector<int>(mySize, 1));
        grid=new int*[mySize];
        size = mySize;
        for(int i=0; i < mySize; i++){
            grid[i] = new int[mySize];
        }
        for(int r=0; r < mySize; r++){
            for(int c=0; c < mySize; c++){
                grid[r][c] = 1;
            }
        }
        
    }
    ~Grid(){
        for(int i=0; i < size; i++){
            delete[] grid[i];
        }
        delete[] grid;
        cout << "Grid deleted"<<endl;
    }
    void drawPixel(IntegerPoint &p1)
    {
        if (p1.x < size && p1.y < size && p1.x >=0 && p1.y >=0)
            grid[p1.x][p1.y] = 0;
    }
    void drawPixel(int x, int y)
    {
        if (x < size && y < size && x >=0 && y >= 0)
            grid[x][y] = 0;
    }
    // void drawPixel(int x,int y,int c){
    //     if (x < size && y < size && x >=0 && y >= 0)
    //         grid[x][y] = 2;
    // }
    int getSize()
    {
        return size;
    }
    void drawLine(int x1,int y1,int x2,int y2){
        IntegerPoint p1(x1,y1);
        IntegerPoint p2(x2,y2);
        drawLine(p1,p2);
    }
    /*
        Given two unscaled points, draws line on grid that goes through these two points
    */
    void drawLine(DoublePoint& d1, DoublePoint& d2){
        double dx = d2.x - d1.x;
        double dy = d2.y-d1.y;
        if(dx==0){
            dx = 0.00000001;
        }
        double slope = dy/dx;
        //line equation: slope * (x-d1.x) + d1.y
        // double yintercept = slope *(-d1.x)+d1.y;
        // double xintercept = -d1.y/slope + d1.x;
        // int scaledY = (int)(round(yintercept*size));
        // int scaledX = (int)(round(xintercept*size));

        // cout <<"y intercept: "<<yintercept<<endl;
        // cout << "xintercept: "<<xintercept<<endl;
        cout << "\nslope " << slope<<endl;
        double y0 = slope * (-1.0*d1.x) + d1.y;
        double y1 = slope * (1.0-d1.x )+d1.y;
        DoublePoint start(0.0,y0);
        DoublePoint end(1,y1);
        IntegerPoint startInt = start.scaleToInt(size);
        IntegerPoint endInt = end.scaleToInt(size);
        drawLine(startInt,endInt);
        //drawLine(scaledX,0,0,scaledY);
        IntegerPoint i1 = d1.scaleToInt(size);
        IntegerPoint i2 = d2.scaleToInt(size);
        drawLine(i1,i2);
    }
    void drawLine(IntegerPoint& p1, IntegerPoint& p2)
    {
        // determine case
        int dx = p2.x - p1.x;
        int dy = p2.y - p1.y;
        if (abs(dx) >= abs(dy))
        {
            if (dx > 0 && dy > 0)
            {
                drawDAxPlus(p1.x, p1.y, p2.x, p2.y);
            }
            else if (dx < 0 && dy < 0)
            {
                drawDAxPlus(p2.x, p2.y, p1.x, p1.y);
            }
            else if (dy == 0)
            {
                drawHoriz(p1.x, p1.y, p2.x, p2.y);
            }
            else if (dx > 0 && dy < 0)
            {
                drawDAxMinus(p1.x, p1.y, p2.x, p2.y);
            }
            else if (dx < 0 && dy > 0)
            {
                drawDAxMinus(p2.x, p2.y, p1.x, p1.y);
            }
        }
        else
        {
            if (dx > 0 && dy > 0)
            {
                drawDAyPlus(p1.x, p1.y, p2.x, p2.y);
            }
            else if (dx < 0 && dy < 0)
            {
                drawDAyPlus(p2.x, p2.y, p1.x, p1.y);
            }
            else if (dx == 0)
            {
                drawVert(p1.x, p1.y, p2.x, p2.y);
            }
            else if (dx < 0 && dy > 0)
            {
                drawDAyMinus(p1.x, p1.y, p2.x, p2.y);
            }
            else if (dx > 0 && dy < 0)
            {
                drawDAyMinus(p2.x, p2.y, p1.x, p1.y);
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
   void toPPM(){
        ofstream myfile;
        myfile.open("triangle.ppm");
        myfile <<"P3"<<" "<<size<<" "<<size<<" "<<1 <<endl;
        for(int i=0; i < size; i++){
            for(int j=0; j < size;j++){
                if(grid[i][j]==1){
                    myfile <<"1 1 1 ";
                }else if(grid[i][j]==0){
                    myfile << "0 0 0 ";
                }else{
                    myfile << "1 0 0";
                }
            }
            myfile<<endl;
        }
        myfile.close();
   }
};
class Circle
{
private:
    DoublePoint *center;
    double radius;

public:
    Circle(double x, double y, double rad)
    {
        center = new DoublePoint(x, y);
        radius = rad;
    }
    ~Circle(){
        delete center;
        cout << "Circle deleted"<<endl;
    }
    DoublePoint& getCenter(){
        return *center;
    }
    double getRadius(){
        return radius;
    }
    string toString()
    {
        return "Center: " + center->toString() + "; radius=" + to_string(radius);
    }
    /*
        Draws the circle on a grid;
    */
    void drawCircle(Grid &gr)
    {
        int size = gr.getSize();
        IntegerPoint intCenter = center->scaleToInt(size);
        //cout << center->toString()<<" size: " << size<<endl;
        //cout << intCenter.toString() <<endl;
        /*
            Set initial values of (xc, yc) and (x, y)
            Set decision parameter d to d = 3 – (2 * r). 
            
            call drawCircle(int xc, int yc, int x, int y) function.
            Repeat steps 5 to 8 until x < = y
            Increment value of x.
            If d < 0, set d = d + (4*x) + 6
            Else, set d = d + 4 * (x – y) + 10 and decrement y by 1.
            call drawCircle(int xc, int yc, int x, int y) function
        */
       
        // int r = (int)(radius*size+0.5);
        // int x = 0;
        // int y = r;
        // int d = 3 - 2*r;
        // while(y >= x){
        //     gr.drawPixel(x+intCenter.x, y+intCenter.y);
        //     //cout << x+intCenter.x<<" "<<y+intCenter.y<<endl;
        //     gr.drawPixel(x+intCenter.x, -y+intCenter.y);
        //     gr.drawPixel(-x+intCenter.x, y+intCenter.y);
        //     gr.drawPixel(-x+intCenter.x, -y+intCenter.y);
        //     gr.drawPixel(y+intCenter.x, x+intCenter.y);
        //     gr.drawPixel(y+intCenter.x, -x+intCenter.y);
        //     gr.drawPixel(-y+intCenter.x, x+intCenter.y);
        //     gr.drawPixel(-y+intCenter.x, -x+intCenter.y);
        //     x ++;
        //     if(d < 0){
        //         d += (4*x+6);
        //     }else{
        //         d+= (4*(x-y)+10);
        //         y-= 1;
        //     }
        // }


        
        //midpoint circle algorithm
        int r = (int)(round(radius * size));
        int x, y, y2, y2_new, ty;
        vector<IntegerPoint> lastPts = vector<IntegerPoint>();
        //xmax = (int)(round(radius*size * 0.70710678)); // maximum x at radius/sqrt(2) (go to 46 degrees )
        y = r;
        y2 = y * y;
        ty = (2 * y) - 1;
        //ty = 3-2*y;
        y2_new = y2;
        //int yfinal = y;
        for (x = 0; x <= y; x++)
        {
            if ((y2 - y2_new) >= ty)
            {
                y2 -= ty;
                y -= 1;
                ty -= 2;
            }
            //shift center from given algorithm
            gr.drawPixel(x+intCenter.x, y+intCenter.y);
            //cout << x+intCenter.x<<" "<<y+intCenter.y<<endl;
            gr.drawPixel(x+intCenter.x, -y+intCenter.y);
            gr.drawPixel(-x+intCenter.x, y+intCenter.y);
            gr.drawPixel(-x+intCenter.x, -y+intCenter.y);
            gr.drawPixel(y+intCenter.x, x+intCenter.y);
            gr.drawPixel(y+intCenter.x, -x+intCenter.y);
            gr.drawPixel(-y+intCenter.x, x+intCenter.y);
            gr.drawPixel(-y+intCenter.x, -x+intCenter.y);
            y2_new -= (2 * x) - 3;
            //yfinal = y;
            //cout << gr.toString()<<endl;
        }
        //Solution to gaps: fill the start and end points with lines
        //gr.drawLine(xmax+intCenter.x,yfinal+intCenter.y,yfinal+intCenter.x,xmax+intCenter.y);
        //gr.drawLine(xmax+intCenter.x,-yfinal+intCenter.y,yfinal+intCenter.x,-xmax+intCenter.y);
        //gr.drawLine(-xmax+intCenter.x,yfinal+intCenter.y,-yfinal+intCenter.x,xmax+intCenter.y);
        //gr.drawLine(-xmax+intCenter.x,-yfinal+intCenter.y,-yfinal+intCenter.x,-xmax+intCenter.y);
    }
};
class Triangle
{
private:
    DoublePoint *p1=NULL;
    DoublePoint *p2=NULL;
    DoublePoint *p3=NULL;
    double l1;
    double l2;
    double l3;
    double semiperimeter;
    double ccR;
    double icR;

public:
    /*
        Randomly instantiates three points and finds length of edges
    */
    Triangle()
    {   
        double x1 = double(rand()) / RAND_MAX;
        double y1 = double(rand()) / RAND_MAX;
        double x2 = double(rand()) / RAND_MAX;
        double y2 = double(rand()) / RAND_MAX;
        double x3 = double(rand()) / RAND_MAX;
        double y3 = double(rand()) / RAND_MAX;
        p1 = new DoublePoint(x1,y1);
        p2 = new DoublePoint(x2,y2);
        p3 = new DoublePoint(x3,y3);
        l1 = sqrt((p3->x - p2->x) * (p3->x - p2->x) + (p3->y - p2->y) * (p3->y - p2->y));
        l2 = sqrt((p1->x - p3->x) * (p1->x - p3->x) + (p1->y - p3->y) * (p1->y - p3->y));
        l3 = sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));
        // the check below ensures that triangle is never equilateral or have same points
        semiperimeter = 0.5 * (l1 + l2 + l3);
        icR = sqrt((semiperimeter - l1) * (semiperimeter - l2) * (semiperimeter - l3) / semiperimeter);
        ccR = l1 * l2 * l3 / (4 * semiperimeter * icR);
        //generate points until not equilateral
        while (l1 < 0.000000001 || l2 < 0.00000001 || l3 < 0.000000001 || (abs(l1 - l2) < 0.000000001 && abs(l2 - l3) < 0.000000001) || icR < 0.00000000001)
        {
            delete p1;
            delete p2;
            delete p3;
            double x1 = double(rand()) / RAND_MAX;
            double y1 = double(rand()) / RAND_MAX;
            double x2 = double(rand()) / RAND_MAX;
            double y2 = double(rand()) / RAND_MAX;
            double x3 = double(rand()) / RAND_MAX;
            double y3 = double(rand()) / RAND_MAX;
            p1 = new DoublePoint(x1,y1);
            p2 = new DoublePoint(x2,y2);
            p3 = new DoublePoint(x3,y3);
            l1 = sqrt((p3->x - p2->x) * (p3->x - p2->x) + (p3->y - p2->y) * (p3->y - p2->y));
            l2 = sqrt((p1->x - p3->x) * (p1->x - p3->x) + (p1->y - p3->y) * (p1->y - p3->y));
            l3 = sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));
            semiperimeter = 0.5 * (l1 + l2 + l3);
            icR = sqrt((semiperimeter - l1) * (semiperimeter - l2) * (semiperimeter - l3) / semiperimeter);
            ccR = l1 * l2 * l3 / (4 * semiperimeter * icR);
        }
    }
    ~Triangle(){
        delete p1;
        delete p2;
        delete p3;
        cout << "Triangle Deleted"<<endl;
    }
    /*
        Finds the circumcenter of the triangle, returns a Circle Object with double center and radii.
        The circum center a circle which passes through all the vertices of a triangle
    */
    Circle findCircumCenter() 
    {
        // center (x1,y1) and (x2,y2) and slope1 slope2 ref perp bisectors  p1, p2 and p2, p3
        double x1 = (p1->x + p2->x) / 2;
        double y1 = (p1->y + p2->y) / 2;
        double x2 = (p3->x + p2->x) / 2;
        double y2 = (p3->y + p2->y) / 2;
        double dx1 = (p1->x - p2->x);
        double dy1 = (p1->y - p2->y);
        double dx2 = (p3->x - p2->x);
        double dy2 = (p3->y - p2->y);

        if (abs(dy1) == 0)
        {
            dy1 = 0.000000001;
        }
        if (abs(dy2) == 0)
        {
            dy2 = 0.000000001;
        }
        double m1 = -1 * dx1 / dy1;
        double m2 = -1 * dx2 / dy2;
        double xCent = (m1 * x1 - m2 * x2 - y1 + y2) / (m1 - m2);
        // cout<<xCent;
        double yCent = m1 * (xCent - x1) + y1;
        // cout << yCent;
        // radius
        Circle newCircle(xCent, yCent, ccR);
        return newCircle;
    }
    /*
        Returns the incenter of triangle, center of inscribed circle
    */
    Circle findInCenter()
    {
        double x = (l1 * p1->x + l2 * p2->x + l3 * p3->x) / (l1 + l2 + l3);
        double y = (l1 * p1->y + l2 * p2->y + l3 * p3->y) / (l1 + l2 + l3);
        return Circle(x, y, icR);
    }
    /*
        Calculate the orthocenter of the triangle(intersection of altitudes)
    */
    DoublePoint findOrthoCenter(){
                // center (x1,y1) and (x2,y2) and slope1 slope2 ref altitudes to segments p1, p2 and p2, p3
        double x1 = p3->x; // oppsite point is p3
        double y1 = p3->y; 
        double x2 = p1->x; //opposite point is p1
        double y2 = p1->y;
        double dx1 = (p1->x - p2->x);
        double dy1 = (p1->y - p2->y);
        double dx2 = (p3->x - p2->x);
        double dy2 = (p3->y - p2->y);

        if (abs(dy1) == 0)
        {
            dy1 = 0.000000001;
        }
        if (abs(dy2) == 0)
        {
            dy2 = 0.000000001;
        }
        double m1 = -1 * dx1 / dy1;
        double m2 = -1 * dx2 / dy2;
        double xCent = (m1 * x1 - m2 * x2 - y1 + y2) / (m1 - m2);
        // cout<<xCent;
        double yCent = m1 * (xCent - x1) + y1;
        return DoublePoint(xCent,yCent);
        // cout << yCent;
        // radius
        

    }
    /*
        The Nine-point circle passes through each of the mid points and altitude foots
    */
    Circle findNinePoint(DoublePoint& ortho,Circle& circumcirle){
        double nineptx = 0.5*(ortho.x+circumcirle.getCenter().x);
        double ninpty = 0.5*(ortho.y+circumcirle.getCenter().y);
        double radius = 0.5*circumcirle.getRadius();
        return Circle(nineptx,ninpty,radius);
    }
    string toString()
    {
        return "Triangle:" + p1->toString() + ", " + p2->toString() + ", " + p3->toString();
    }

    void drawTriangle(Grid &gr){
        int size = gr.getSize();
        IntegerPoint intp1 = p1->scaleToInt(size);
        IntegerPoint intp2 = p2->scaleToInt(size);
        IntegerPoint intp3 = p3->scaleToInt(size);
        gr.drawLine(intp1,intp2);
        gr.drawLine(intp2,intp3);
        gr.drawLine(intp3,intp1);
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
//     double crossProd = (p1.x-v1.x)*(p2.y - p1.y)* p2.x;
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
double dotProduct(DoublePoint &v1, DoublePoint &v2){
    return v1.x * v2.x + v1.y*v2.y;
}
bool sameSideChecker(DoublePoint p, DoublePoint a, DoublePoint b, DoublePoint c){
    DoublePoint v0(c.x-a.x,c.y-a.y);
    DoublePoint v1(b.x-a.x,b.y-a.y);
    DoublePoint v2(p.x-a.x,p.y-a.y);
    double dot00, dot01, dot02, dot11, dot12;
    dot00 = dotProduct(v0,v0);
    dot01 = dotProduct(v0,v1);
    dot02 = dotProduct(v0,v2);
    dot11 = dotProduct(v1,v1);
    dot12 = dotProduct(v1,v2);
    double denom = 1.0/(dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * denom;
    double v = (dot00 * dot12 - dot01 * dot02) * denom;
    return u >=0 && v >=0 && (u+v) < 1;
}
/*
    - Uses the determinant algorithm to check containment
    for each point, determine the sign of c1, c2, c3 are the same
    - return true if the triangle contains the point at i
*/
bool checkIntersect(const int &i, DoublePoint (&dplist)[4]){
    vector<DoublePoint> tri = vector<DoublePoint>();
    for(int j=0; j < 4; j++){
        if(j!=i)
            tri.push_back(dplist[j]);
    }
    return sameSideChecker(dplist[i],tri[0],tri[1],tri[2]);
    //math from https://blackpawn.com/texts/pointinpoly/default.html

}
void part1(){
    srand(time(NULL));
    ofstream myLog;
    myLog.open("log.txt");
    DoublePoint dplist[4];
    for(int i=0; i < 4; i++){
        double x = rand()*1.0/RAND_MAX;
        double y = (double)(rand())/RAND_MAX;
        dplist[i] = DoublePoint(x,y);
    }
    for(int i=0; i < 3; i++){
        std::cout << dplist[i].toString();
        myLog << dplist[i].toString();
        if(i < 2){
            std::cout << " , ";
            myLog << " , ";
        }
    }
    // dplist[0].x = 0;
    // dplist[0].y = .5;
    // dplist[1].x=.5;
    // dplist[1].y=1;
    // dplist[2].x = 1;
    // dplist[2].y = .4;
    //DoublePoint p4(rand()/RAND_MAX,rand()/RAND_MAX);
    bool inTriangle = true;
    /*
        PIP algorithm - count the number of times horizontal ray of p4 crosses the sides.
        If Even number of times, then point is outside triangle
    */
   //while point overlaps triangle
   //1. find intersection of x and the segment
   //2. if POI is greater than x, within bounds of segment, then increment count
   //if count is even, break
   cout << "\n";
   myLog << "\n";
    while(inTriangle){
        //int crosscount = 0;
        dplist[3].x=(double)(rand())/RAND_MAX;
        dplist[3].y = (double)(rand())/RAND_MAX;
        cout << "testing point "<<dplist[3].toString()<<endl;
        myLog << "testing point "<<dplist[3].toString()<<endl;
        bool intersect = false;
        for(int i=0; i < 4;i++){
            if(checkIntersect(i,dplist)){
                intersect = true;
                break;
            }
        }
        if(!intersect){
            break;
        }
        // for(int i=0; i < 4; i++){
        //     std::cout << dplist[i].toString()<<",";
        // }
        // std::cout << "\n";
        // string s;

        // cin >> s;
        
    }
    cout << "\n"<<"Done. Points: "<<endl;
    myLog << "\n"<<"Done. Points: "<<endl;
    //cout.precision(15);
    //myLog.precision(15);
    ofstream myPoints;
    myPoints.open("points.txt");
    myPoints.precision(17);
    for(int i=0; i < 4; i++){
        std::cout << "("<<dplist[i].x<<","<<dplist[i].y<<")";
        myLog << "("<<dplist[i].x<<","<<dplist[i].y<<")";
        myPoints << fixed << "("<<dplist[i].x<<","<<dplist[i].y<<")";
        if(i < 3){
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
vector<DoublePoint> split(string& myLine){
    // int start = 0;
    // vector<string> toRet;
    // size_t delimLen = delim.size();
    
    // while(myLine.find(delim)!=)
    string toNum = "";
    //get the list of doubles
    vector<double> ptvals;
    for(char ch: myLine){
        if(ch=='.'|| isdigit(ch)){
            toNum += ch;
        }
        else{
            ptvals.push_back(stod(toNum));
            toNum = "";
        }
    }
    for(int i=0; i < ptvals.size(); i+=2){
    }
    

}
void part2(){
    //Read the file
    fstream myFile("points.txt",ios::in);
    
    string line;
    getline(myFile,line); //get the line
    //split contents of line by " , "
    //cout << line << "end of line";

    //calculate the points
    //draw the points
}
int main(){
   
    //part1();
    part2();
}