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
bool checkIntersect(int i, DoublePoint (&dplist)[4]){
    int crosscount = 0;
    for(int j=0; j < 4; j++){
            //maximum
            // double miny = min(dplist[j].y,dplist[(j+1)%3].y);
            // double maxy = max(dplist[j].y, dplist[(j+1)%3].y);
            // cout << miny << " "<<maxy<<endl;
            // if(dplist[3].y <= maxy && dplist[3].y >=miny){
            //     crosscount ++;
            // }
            if(j==i) continue;
            double x1 = dplist[j].x;
            double y1 = dplist[j].y;
            int next = j + 1;
            //this checks for the next value being equal to i
            if(next == i){
                next++;
            }
            double x2 = dplist[(next)%4].x;
            double y2 = dplist[(next)%4].y;
            double dx = x2 - x1;
            double dy = y2 - y1;
            double xinter;
            if(abs(dy) < 1e-15){
                xinter=x1;
                continue;
            }
            double slope = dy/dx;
            if(abs(slope) < 1.0e-15){
                continue;
            }

            xinter = (dplist[i].y - y1)/slope + x1;
            //std::cout <<"Integersection: "<<xinter << " "<<dplist[3].y<<endl;
            if(xinter <= max(x1,x2) && xinter >= min(x1,x2)){
                if(xinter > dplist[i].x){
                    std::cout <<xinter << " "<<dplist[i].y<<endl;
                    crosscount ++;
                }
            }        
    }
    return crosscount%2!=0;
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
    for(int i=0; i < 4; i++){
        double x = rand()*1.0/RAND_MAX;
        double y = (double)(rand())/RAND_MAX;
        dplist[i] = DoublePoint(x,y);
    }
    //pts: (0.55391,0.379061) , (0.884164,0.279246) , (0.102818,0.299981) , (0.798166,0.667823)
    dplist[0].x = 0.55391;
    dplist[0].y = 0.379061;
    dplist[1].x=0.884164;
    dplist[1].y=0.279246;
    dplist[2].x = 0.102818;
    dplist[2].y = 0.299981;
    dplist[3].x = 0.798166;
    dplist[3].y = 0.667823;
    bool intersect=false;
    for(int i=0; i < 4; i++){
        cout << i<<endl;
        if(checkIntersect(i,dplist)){
            intersect=true;
        }
    }
    cout << "Intersect?" << intersect<<endl;
    
    //cout << "Hello World";
    return 0;
}