#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <iomanip>
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
double crossProduct(DoublePoint &p1, DoublePoint &p2){
    double crossProd = p1.x*p2.y - p1.y* p2.x;
    cout << "p1: "<<p1.toString()<< " ";
    cout << " p2 "<<p2.toString() <<" " << crossProd << endl;
    return crossProd;
}
/*
    - Uses the determinant algorithm to check containment
    for each point, determine the sign of c1, c2, c3 are the same
    - return true if the triangle contains the point at i
*/
bool checkIntersect(const int &i, DoublePoint (&dplist)[4]){
    bool pos = true;
    bool neg = true;
    for(int j=0; j < 4; j++){
        if(i!=j){
            double myProd = crossProduct(dplist[i],dplist[j]);
            pos = pos && (myProd > 0);
            neg = neg && (myProd < 0);
        }
    }
    return pos||neg;
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
        cout << "testing the point "<<dplist[3].toString()<<endl;
        myLog << "testing the point "<<dplist[3].toString()<<endl;
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
    myPoints.precision(15);
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
   
    //part1();
     DoublePoint dplist[4];
    dplist[0].x = 0;
    dplist[0].y = .5;
    dplist[1].x=.5;
    dplist[1].y=1;
    dplist[2].x = 1;
    dplist[2].y = .4;
    dplist[3] = DoublePoint(.5,0);
    // for(int i=0; i < 4; i++){
    //     double x = rand()*1.0/RAND_MAX;
    //     double y = (double)(rand())/RAND_MAX;
    //     dplist[i] = DoublePoint(x,y);
    // }
    // //pts: (0.55391,0.379061) , (0.884164,0.279246) , (0.102818,0.299981) , (0.798166,0.667823)
    // dplist[0].x = 0.55391;
    // dplist[0].y = 0.379061;
    // dplist[1].x=0.884164;
    // dplist[1].y=0.279246;
    // dplist[2].x = 0.102818;
    // dplist[2].y = 0.299981;
    // dplist[3].x = 0.798166;
    // dplist[3].y = 0.667823;
    bool intersect=false;
    for(int i=0; i < 1; i++){
        cout << i<<endl;
        if(checkIntersect(i,dplist)){
            intersect=true;
        }
    }
    cout << "Intersect? " << intersect<<endl;
    
    // //cout << "Hello World";
    // return 0;
}