/**
 * @file l051.cpp
 * @author Vivian Feng
 * @brief Lab 5, Part 1 Canny Edge detection with Sobel and single threshold
 * @version 0.1
 * @date 2022-12-14
 * 
 * @copyright Copyright (c) 2022
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
 * Given two points a and b which lie on a line, and a point c, tells if point c is to right of oriented line AB
 */

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
int main()
{
    part1();
}