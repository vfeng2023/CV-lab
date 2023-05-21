#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <time.h>
#include <fstream>


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
int main()
{
    Grid gr(800);
    // IntegerPoint p1(3, 5);
    // IntegerPoint p2(3, 0);
    // gr.drawLine(p1, p2);
    srand(time(NULL));
    Triangle t;
    cout << t.toString() << endl;
    Circle circumcircle = t.findCircumCenter(); 
    cout << "CircumCenter: "<<circumcircle.toString() << endl;
    Circle incircle = t.findInCenter();
    cout << "Incenter: " << incircle.toString() << endl;
    DoublePoint ortho = t.findOrthoCenter();
    cout << "OrthoCenter: "<<ortho.toString()<<endl;
    Circle ninepoint = t.findNinePoint(ortho,circumcircle);
    cout << "Center of Ninepoint: " << ninepoint.toString();
    t.drawTriangle(gr);
    circumcircle.drawCircle(gr);
    incircle.drawCircle(gr);
    ninepoint.drawCircle(gr);
    //Euler line. Goes through circumcenter and nine-point center 
    gr.drawLine(circumcircle.getCenter(),ninepoint.getCenter());
    // Circle c(0.5,0.5,0.4);
    // c.drawCircle(gr);
    //issues to ask about: ninepoint circle float is correct up to 5 decimal places, but drawing fails to intersect triangle
    //circle drawing algorithm is inprecise
    gr.toPPM();




    // cout <<gr.toString();
    return 0;
}