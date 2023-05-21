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

#define _USE_MATH_DEFINES


using namespace std::chrono;
using namespace std;

const int PRECISION = 23;
// typedef std::numeric_limits< double > dbl;
// ofstream results;
class Grid
{
private:
    int width;
    int height;
    int maxIntensity;
    vector<vector<int>> grid;

public:
    Grid(int width, int height,int maxIntensity)
    {
        grid = vector<vector<int>>(height, vector<int>(width, 0));
        //grid = new int *[mySize];
        //size = mySize;
        this->height = height;
        this->width=width;
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
    void setPixel(int r, int c, int val){
        grid[r][c] = val;
    }
    int getPixel(int r, int c){
        return grid[r][c];
    }
    int getHeight(){
        return height;
    }
    int getWidth(){
        return width;
    }
    void toPPM(string filename)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "P3"<< " " << width << " " << height << " " << maxIntensity << endl;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                myfile << grid[i][j] << " "<< grid[i][j] << " "<<grid[i][j]<<" ";
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
Grid* readImage(string filename){
    //open image file
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
    Grid* gr = new Grid(width,height,intensity);
    for(int r=0; r < height; r++){
        for(int c=0; c < width; c++){
            int red,blue,green;
            imgFile >> red >> blue >> green;
            int gray = (red + blue + green)/3;
            gr->setPixel(r,c,gray);
        }
    }
    imgFile.close();
    return gr;


}
/**
 * @brief applies the sobel operator to the image.
 * 
 * @param gr Pointer to grid on the heap
 * @param threshold threshold to evaluate |G| = sqrt(Gx^2 + Gy^2)
 */
void sobelOperator(Grid *gr,int threshold){
    Grid& myGrid = *gr;
    vector<vector<int>> gx = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> gy = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> xdirsobel = 
    {{-1,0,1},
    {-2,0,2},
    {-1,0,1}};
    vector<vector<int>> ydirsobel = 
    {{-1,-2,-1},
    {0,0,0},
    {1,2,1}};
    //fill gx
    for(size_t r=1; r < gx.size()-1; r++){
        for(size_t c=1; c < gx[0].size()-1; c++){
            int val = 0;
            int gxstartrow = r-1;
            int gxstartcol = c-1;
            for(size_t i=0; i < xdirsobel.size();i++){
                for(size_t j=0; j < xdirsobel.size();j++){
                    int temp =myGrid.getPixel(gxstartrow+i,gxstartcol+j)*xdirsobel[i][j];
                    val += temp;
                }
            }
            gx[r][c] = val;
        }
    }
    //fill gy
    for(size_t r=1; r < gy.size()-1; r++){
        for(size_t c=1; c < gy[0].size()-1; c++){
            int val = 0;
            int gystartrow = r-1;
            int gystartcol = c-1;
            for(size_t i=0; i < ydirsobel.size();i++){
                for(size_t j=0; j < ydirsobel.size();j++){
                    int temp =myGrid.getPixel(gystartrow+i,gystartcol+j)*ydirsobel[i][j];
                    val += temp;
                }
            }
            gy[r][c] = val;
        }
    }
    //compute magG
    vector<vector<int>> magG = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    for(size_t r=0; r < magG.size();r++){
        for(size_t c=0; c < magG[0].size();c++){
            magG[r][c] = gx[r][c]*gx[r][c] + gy[r][c] * gy[r][c];
        }
    }
    //apply threshold
    int thresholdsquared = threshold*threshold;
    Grid edges(myGrid.getWidth(),myGrid.getHeight(),1);
    for(size_t r=0; r < magG.size();r++){
        for(size_t c=0; c < magG[0].size();c++){
            if(magG[r][c] < thresholdsquared){
                edges.setPixel(r,c,0);
            }else{
                edges.setPixel(r,c,1);
            }
        }
    }
    edges.toPPM("imagem.ppm");


}
void part1(){
    // readImage();
    Grid* gr = readImage("image.ppm");
    gr->toPPM("imageg.ppm");
    sobelOperator(gr,275);
    delete gr;
}


//PART 2 starts here
void sobelOperatorPt2(Grid *gr,int low,int high, string of, string f2, string ff);
void areaFill(vector<vector<int>>& gr,int row,int col);
void computeTangents(vector<vector<int>> &gx, vector<vector<int>> &gy, vector<vector<int>> &result);
//void convertToGrid(vector<vector<int>> *v);

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

void sobelOperatorPt2(Grid *gr,int low,int high, string of, string f2, string ff){
    Grid& myGrid = *gr;
    vector<vector<int>> gx = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> gy = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> xdirsobel = 
    {{-1,0,1},
    {-2,0,2},
    {-1,0,1}};
    vector<vector<int>> ydirsobel = 
    {{-1,-2,-1},
    {0,0,0},
    {1,2,1}};
    //fill gx
    for(size_t r=1; r < gx.size()-1; r++){
        for(size_t c=1; c < gx[0].size()-1; c++){
            int val = 0;
            int gxstartrow = r-1;
            int gxstartcol = c-1;
            for(size_t i=0; i < xdirsobel.size();i++){
                for(size_t j=0; j < xdirsobel.size();j++){
                    int temp =myGrid.getPixel(gxstartrow+i,gxstartcol+j)*xdirsobel[i][j];
                    val += temp;
                }
            }
            gx[r][c] = val;
        }
    }
    //fill gy
    for(size_t r=1; r < gy.size()-1; r++){
        for(size_t c=1; c < gy[0].size()-1; c++){
            int val = 0;
            int gystartrow = r-1;
            int gystartcol = c-1;
            for(size_t i=0; i < ydirsobel.size();i++){
                for(size_t j=0; j < ydirsobel.size();j++){
                    int temp =myGrid.getPixel(gystartrow+i,gystartcol+j)*ydirsobel[i][j];
                    val += temp;
                }
            }
            gy[r][c] = val;
        }
    }
    //compute magG
    vector<vector<int>> magG = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    for(size_t r=0; r < magG.size();r++){
        for(size_t c=0; c < magG[0].size();c++){
            magG[r][c] = gx[r][c]*gx[r][c] + gy[r][c] * gy[r][c];
        }
    }
    //grid filled
    int lowsq = low*low;
    int highsq = high*high;
    //apply double threshold, returning pointer to vector on the heap
    vector<vector<int>> edges = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> &ed = edges;
    for(size_t r=0; r < magG.size(); r++){
        for(size_t c=0; c < magG[0].size(); c++){
            if(magG[r][c] < lowsq){
                ed[r][c] = 0;
            }else if(magG[r][c] < highsq){
                ed[r][c] = 1;
            }else{
                ed[r][c] = 2;
            }
        }
    }
    cout << "Grid filled" << endl;
    //apply area fill to this
    for(size_t r=0; r < ed.size(); r++){
        for(size_t c=0; c < ed.size(); c++){
            if(ed[r][c]==2){
                areaFill(ed,r,c);
            }
        }
    }
    cout << "Finished areaFill"<< endl;

    //apply nonmax suppression
        //compute angles
    vector<vector<int>> angles = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    computeTangents(gx,gy,angles);
        //produce image2
    Grid img2 = Grid(myGrid.getWidth(),myGrid.getHeight(),1);
    for(size_t r=1; r < angles.size()-1; r++){
        for(size_t c=1; c < angles[0].size()-1; c++){
            int p1 = magG[r][c];
            int p0=0;
            int p2 = 0;
            if(abs(angles[r][c])==0 or abs(angles[r][c])==180){
                p0 = magG[r][c-1];
                p2 = magG[r][c+1];
            }else if(abs(angles[r][c])==90){
                p0 = magG[r-1][c];
                p2 = magG[r+1][c];
            }else if(abs(angles[r][c])==45){
                p0 = magG[r-1][c+1];
                p2 = magG[r+1][c-1];
            }else if(abs(angles[r][c])==135){
                p0 = magG[r-1][c-1];
                p2 = magG[r+1][c+1];
            }
            if(p1 >= p0 && p1 >= p2){
                img2.setPixel(r,c,1);
            }
        }
    }
    img2.toPPM(f2);
        //place in final grid
    //grid will be filled with 0, 1, 2, 3. This is where part 2 formerly was used to create image1.ppm
    Grid newGrid = Grid(myGrid.getWidth(),myGrid.getHeight(),1);
    for(int i=0; i < newGrid.getHeight(); i++){
        for(int j=0; j < newGrid.getWidth(); j++){
            if(ed[i][j] == 0|| ed[i][j] == 1){
                newGrid.setPixel(i,j,0);
            }else{
                newGrid.setPixel(i,j,1);
            }
        }
    }
    newGrid.toPPM(of);

    //create final image by doing and of the values in img2 and newGrid
    Grid finalImage = Grid(myGrid.getWidth(),myGrid.getHeight(),1);
    for(int r=0; r < myGrid.getHeight();r++){
        for(int c=0; c < myGrid.getWidth();c++){
            if(img2.getPixel(r,c)==1 && newGrid.getPixel(r,c)==1){
                finalImage.setPixel(r,c,1);
            }
        }
    }
    //newGrid.toPPM(of);
    finalImage.toPPM(ff);
    
}



/**
 * @brief accordingly changes weak edges in the vector gr to 0 and 2. Marks visited squares with 3
 * 
 * @param gr a reference to the vector
 * @param row starting row
 * @param col starting column
 */
void areaFill(vector<vector<int>>& gr,int row,int col){
    //cout << row << " "<< endl;
    if(row < 0 || col < 0 || row >= (int)gr.size()|| col >= (int)gr[0].size()||gr[row][col] > 2|| gr[row][col] == 0){
        return;
    }
    // for(int i=0; i < gr.size(); i++){
    //     for(int j = 0; j < gr[0].size(); j++){
    //         cout << gr[i][j] << " ";
    //     }
    //     cout << "\n";
    // }
    // cout << "end iteration";
    if(gr[row][col]==1){
        gr[row][col] = 2;
    }else if(gr[row][col]==2){
        gr[row][col] = 3;
    }
    areaFill(gr,row+1,col);
    areaFill(gr,row+1,col+1);
    areaFill(gr,row, col+1);
    areaFill(gr,row-1,col+1);
    areaFill(gr,row-1,col);
    areaFill(gr,row-1,col-1);
    areaFill(gr,row,col-1);
    areaFill(gr,row+1,col-1);
}

void part3(const int &argc, char** &argv){
    string imageName="image.ppm";
    int low = 275;
    int high = 350;
    string fg = "imageg.ppm";
    string of = "image1.ppm";
    string f2 = "image2.ppm";
    string ff = "imagef.ppm";
    //read the command line arguments
    for(int i=1; i < argc; i++){
        string flag=string(argv[i]);
        if(flag=="-f"){
            //cout << "flag read" << endl;
            imageName = argv[i+1];
            
        }else if(flag=="-lt"){
            stringstream ss(argv[i+1]);
            ss >> low;
        }else if(flag=="-ht"){
            stringstream ss(argv[i+1]);
            ss >> high;
        }else if(flag=="-of"){
            of = argv[i+1];
        }else if(flag=="-f2"){
            f2 = argv[i+1];
        }else if(flag=="-ff"){
            ff = argv[i+1];
        }else if(flag=="-fg"){
            fg = argv[i+1];
        }
        i++;
    }
    //print the values of the flags
    cout << imageName << " " << low << " " << high << " " << fg << " "<<of<<" "<< f2 << " " << ff<< endl;
    Grid* gr = readImage(imageName);
    gr->toPPM(fg);
    sobelOperatorPt2(gr,low, high,of,f2,ff);
    delete gr;

}
/**
 * @brief COmpute the angles of the image pixels and puts them in result
 * 
 * @param gx 
 * @param gy 
 * @param result 
 */
void computeTangents(vector<vector<int>> &gx, vector<vector<int>> &gy, vector<vector<int>> &result){
    for(size_t r=0; r < gx.size(); r++){
        for(size_t c=0; c < gx[0].size();c++){
            double angle = atan2(gy[r][c],gx[r][c])*180/M_PI;
            // int scaled = (int)((angle + 180)*8/360+0.4);
            // result[r][c] = scaled*45-180;
            //if angle <=22.5 and 
            double absangle = abs(angle);
            if(absangle <= 22.5){
                result[r][c] = 0;
            }else if(absangle <=45){
                result[r][c] =45;
            }else if(absangle <= 67.5){
                result[r][c] = 45;
            }else if(absangle < 90){
                result[r][c] = 90;
            }else if(absangle < 112.5){
                result[r][c] = 90;
            }else if(absangle < 135){
                result[r][c] = 135;
            }else if(absangle <= 157.5){
                result[r][c] = 135;
            }else{
                result[r][c] = 180;
            }
        }
    }
}

int main(int argc, char** argv)
{
    // part1();
    part3(argc, argv);
}