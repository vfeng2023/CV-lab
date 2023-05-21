#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <vector>
//#include <iomanip>

using namespace std;


class DoublePoint
{
public:
    double x;
    double y;

public:
    DoublePoint(){}
    DoublePoint(double xp, double yp)
    {
        x = xp;
        y = yp;
    }
    // void set(double xp, double yp){
    //     x = xp;
    //     y = yp;
    // }
    ~DoublePoint(){}
    string toString()
    {
        return ("(" + to_string(x)) + "," + to_string(y) + ")";
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
int main(){
   
    part1();
}