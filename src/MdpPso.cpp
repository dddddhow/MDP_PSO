#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <random>
#include <ctime>
#include <armadillo>

using namespace arma;
using namespace std;

/*============================================================================*/
void itoa(int i,char* aa)
{
    int power,j;
    j=i;
    for(power=1;j>=10;j/=10)
        power*=10;
    for(;power>0;power/=10)
    {
        *aa++='0'+i/power;
        i%=power;
    }
    *aa='\0';
}
/*============================================================================*/



//----------------------------------------------------------------------//
//                     Common Parameters                                //
//----------------------------------------------------------------------//
const  int dim = 2;
const  int Nx  = 200;
const  int Nz  = 8000;
const  int w   = 100;
//extern int iter_now;
//---------------------------------------------------------------------//
//              Particle Swarm Optimizatin Parameters                  //
//---------------------------------------------------------------------//
const int p_num = 70;//粒子数量
const int iters = 200;//迭代次数
//Definition: The limit of position and velocity of every partilce
const  int spdx_max = 2;
const  int spdx_min = 1;
const  int spdz_max = 10;
const  int spdz_min = -10;

const  int posx_max = Nx-1;
const  int posx_min = 0;
const  int posz_max = Nz-1;
const  int posz_min = 0;
//Definition: The position and velocity of every particle
Cube<int> pos(iters, p_num, dim);
Cube<int> spd(iters, p_num, dim);
//Definition: The best location of each particle(p_best) and all particles(g_best)
Cube<int> p_best_location(iters, p_num, dim);
Mat<int> g_best_location(iters, dim);
Mat<float>  p_best_rew(iters, p_num);
Mat<float>  g_best_rew(iters, 1);
//----------------------------------------------------------------------//
//       Markov Decision Process Parameter Definition Parameters        //
//----------------------------------------------------------------------//

//Definition: The Action(Act) and it's probability(Prob_A)
Mat<float> Act(1, 4);
Mat<float> Prob_A(1, 4);

//Definition: The Probability(Prob)
Mat<float> Prob(1, 4);

//Definition: The Reward(Rew)
Mat<float> Rew(iters, p_num);

//Definition: The seismic data
Mat<float> seis(Nx, Nz);

//----------------------------------------------------------------------//
//                      Function Declaration                            //
//----------------------------------------------------------------------//
//Initialization
void init_func();

//Look for g/p_best_rew && g/p_best_location
void rew_func(int iter_now);

//Calculate the g/p_best_rew && g/p_best_location
void fitness_func(int iter_now, int p_num_now);

//Iterative computation
void Pso_func();

//----------------------------------------------------------------------//
//                              Calculation                             //
//----------------------------------------------------------------------//

int main()
{
    //--------------------------------//
    //       Seismic Data Input       //
    //--------------------------------//
    Mat<float> V_p_temp;
    V_p_temp.load("../file/record",raw_binary);
    V_p_temp.reshape(Nz,Nx);


    for (int i=0; i<Nx; i++)
    {
        for(int j=0; j<Nz; j++)
        {
            seis(i,j) = V_p_temp(j,i);
        }
    }

    //--------------------------------//
    //       Seismic Data Process1    //
    //--------------------------------//

    Mat<float> seis_ori(Nz,Nx);
    for (int i=0; i<Nz; i++)
    {
        for(int j=0; j<Nx; j++)
        {
            seis_ori(i,j) = V_p_temp(i,j);
        }
    }

    vec seis_mul(w);
    seis_mul.zeros();
    float max_val;
    int index_max_val;
    vec pol(Nz);
    int ik = 0;

    for(int j=5750; j<5900; j=j+w)
    {

        for(int i=0;i<w;i++)
        {
            seis_mul(i) = fabs(seis_ori(i+j,0));
        }
        max_val       = max(seis_mul);
        index_max_val = index_max(seis_mul)+j;

        if( index_max_val < j+w && index_max_val > j)
        {
            if(seis_ori(index_max_val,0) > 0.0)
            {
                pol(ik) = 100;
                ik++;
                cout << index_max_val << " : ";
                cout << seis_ori(index_max_val,0)<<endl;
            }
            else if(seis_ori(index_max_val,0) < -0)
            {
                pol(ik) = -100;
                ik++;
                cout << index_max_val << " : ";
                cout << seis_ori(index_max_val,0)<<endl;
            }
        }

    }

    float temp_max = 0;
    for(int iz=index_max_val-50; iz<index_max_val+50;iz++)
    {
        temp_max = temp_max + abs(seis_ori(iz,0));
    }
    cout<<temp_max<<endl;

    //--------------------------------//
    //       Seismic Data Process2    //
    //--------------------------------//
    init_func();

    Pso_func();

    //--------------------------------//
    //        Result Data save        //
    //--------------------------------//

    //remove old file
    remove("../file/seis");
    remove("../file/g_best_location");
    remove("../file/col");


    //save file
    seis.save("../file/seis",raw_binary);

    Mat<float> seis_col(Nz,1);
    for(int i=0; i<Nz; i++){seis_col(i,0)=seis(100,i);}
    //seis_col.save("../file/col",raw_binary);

    g_best_location.save("../file/g_best_location",raw_binary);


    /*
       for(int i=0; i<10; i++)
       {

       cout<<"迭代次数"<<"      ";
       cout<<"粒子"<<"      ";
       cout<<"速度"<<"      ";
       cout<<"位置"<<"          ";
       cout<<"奖励"<<"      ";
       cout<<"自最佳位置"<<"    ";
       cout<<"自最佳奖励"<<"       ";
       cout<<"群最佳位置"<<"       ";
       cout<<"群最佳奖励";
       cout<<endl;

       for(int j=0; j<p_num; j++)
       {


       cout<<"    "<<i<<"          "<<j<<"      ";

       cout<<" [";
       cout<<spd(i,j,0);
       cout<<",";
       cout<<spd(i,j,1);
       cout<<"]    ";

       cout<<" [";
       cout<<pos(i,j,0);
       cout<<",";
       cout<<pos(i,j,1);
       cout<<"]    ";

       cout<<" ";
       cout<<Rew(i,j);
       cout<<"    ";

       cout<<"   [";
       cout<<p_best_location(i,j,0);
       cout<<",";
       cout<<p_best_location(i,j,1);
       cout<<"]    ";

       cout<<"    ";
       cout<<p_best_rew(i,j);
       cout<<"          ";


       cout<<"[";
       cout<<g_best_location(i,0);
       cout<<",";
       cout<<g_best_location(i,1);
       cout<<"]    ";

       cout<<"  ";
       cout<<g_best_rew(i,0);
       cout<<"    ";

       cout<<endl;
       }
       cout<<endl;
       }

*/

    return 0;
}


//----------------------------------------------------------------------//
//                   Initialization Function                            //
//----------------------------------------------------------------------//
void init_func()
{
    int iter_now = 0;
    for(int j=0; j<p_num; j++)
    {

        //pos(iter_now, j,0) = int(rand()%10);
        //pos(iter_now, j,1) = int(rand()%50+2000);

        //spd(iter_now, j,0) = float(rand()%2);
        //spd(iter_now, j,1) = float(rand()%10-5);


        pos(iter_now, j,0) = 0;
        pos(iter_now, j,1) = 5807+j;
        spd(iter_now, j,0) = 1;
        spd(iter_now, j,1) = 5;


        if(spd(iter_now,j,0) > spdx_max) {spd(iter_now,j,0) = spdx_max;}
        if(spd(iter_now,j,0) < spdx_min) {spd(iter_now,j,0) = spdx_min;}
        if(spd(iter_now,j,1) > spdz_max) {spd(iter_now,j,1) = spdz_max;}
        if(spd(iter_now,j,1) < spdz_min) {spd(iter_now,j,1) = spdz_min;}

        if(pos(iter_now,j,0) > posx_max) {pos(iter_now,j,0) = posx_max;}
        if(pos(iter_now,j,0) < posx_min) {pos(iter_now,j,0) = posx_min;}
        if(pos(iter_now,j,1) > posz_max) {pos(iter_now,j,1) = posz_max;}
        if(pos(iter_now,j,1) < posz_min) {pos(iter_now,j,1) = posz_min;}
    }
    rew_func(iter_now);
}

//----------------------------------------------------------------------//
//                      Calculate reward Function                       //
//----------------------------------------------------------------------//
void fitness_func(int iter_now, int p_num_now)
{
    float temp      = 0;
    float fitness   = 0;

    int x_now        = pos(iter_now,p_num_now,0);
    int z_now        = pos(iter_now,p_num_now,1);
    int x_pro,z_pro;

    /*
       if(iter_now-1 >= 0)
       {
       x_pro        = pos(iter_now-1,p_num_now,0);
       z_pro        = pos(iter_now-1,p_num_now,1);
       }
       else
       {
       x_pro        = x_now;
       z_pro        = z_now;
       }
       int mod=1;
       */
    /*if     (x_now-x_pro == 0 && z_now-z_pro == 0){ mod=0;}*/
    //else if(x_now-x_pro == 0 && z_now-z_pro == 1){ mod=1;}
    //else if(x_now-x_pro == 1 && z_now-z_pro == 0){ mod=2;}
    //else if(x_now-x_pro == 1 && z_now-z_pro == 1){ mod=3;}
    /*else {cout << "The Way of particle is wrong!!" <<endl; }*/

    int ix,iz;
    for(int j = 0; j < w; j++)
    {
        ix=x_now;
        iz=z_now+j;

        if(iz >=0 && ix>=0 && iz<Nz && ix<Nx)
        {
            temp = temp + abs(seis(ix,iz));
        }
        else
        {
            temp = temp;
        }
    }

    //fitness = Prob_A(0,mod)*Prob(0,mod)*temp;
    if(temp < 8.52*2.0 && temp >8.52*0.5)
    {
        fitness = temp*10;
    }
    else
    {
        fitness = 0.1;
    }

    //Rew(iter_now,p_num_now) = fitness;

    if(iter_now >= 1)
    {
        Rew(iter_now,p_num_now) = Rew(iter_now-1,p_num_now) + fitness;
    }
    else if(iter_now == 0)
    {
        Rew(iter_now,p_num_now) = fitness;
    }
}


//----------------------------------------------------------------------//
//                  Looking fow max reward Function                     //
//----------------------------------------------------------------------//
void rew_func(int iter_now)
{
    //calculation of fitness function
    for(int i=0; i<p_num; i++)
    {
        fitness_func(iter_now,i);
    }

    //p_best & p_best_location

    float max_rew;
    int max_rew_x, max_rew_z;

    for(int i=0; i<p_num; i++)
    {
        max_rew   = Rew(0,i);
        max_rew_x = pos(0,i,0);
        max_rew_z = pos(0,i,1);

        for(int j=0; j<=iter_now; j++)
        {
            if(max_rew < Rew(j,i))
            {
                max_rew   = Rew(j,i);
                max_rew_x = pos(j,i,0);
                max_rew_z = pos(j,i,1);
            }
        }

        p_best_rew(iter_now,i)        = max_rew;
        p_best_location(iter_now,i,0) = max_rew_x;
        p_best_location(iter_now,i,1) = max_rew_z;
    }

    //g_best & g_best_location
    max_rew   = Rew(iter_now,0);
    max_rew_x = pos(iter_now,0,0);
    max_rew_z = pos(iter_now,0,1);



    for(int i=0; i<p_num; i++)
    {
        if(max_rew < Rew(iter_now,i))
        {
            max_rew   = Rew(iter_now,i);
            max_rew_x = pos(iter_now,i,0);
            max_rew_z = pos(iter_now,i,1);
        }
    }
    g_best_rew(iter_now,0)      = max_rew;
    g_best_location(iter_now,0) = max_rew_x;
    g_best_location(iter_now,1) = max_rew_z;
}


//----------------------------------------------------------------------//
//                       Interation Function                            //
//----------------------------------------------------------------------//
void Pso_func()
{
    int iter_now;

    //Iteration Calculation
    for(int i=1; i<iters; i++)
    {
        iter_now = i;

        for(int j=0; j<p_num; j++)
        {
            //Update of velocity & location of each particles

            //Velocity update

            //spd(iter_now,j,0) = spd(iter_now-1,j,0)+ int(
            //(p_best_location(i-1,j,0) - pos(i-1,j,0))
            //+ (g_best_location(i-1,0) - pos(i-1,j,0)));

            //spd(iter_now,j,1) = spd(iter_now-1,j,1)+ int(
            //(p_best_location(i-1,j,1) - pos(i-1,j,1))
            //+ (g_best_location(i-1,1) - pos(i-1,j,1)));



            float zc = (5950 - pos(i-1,j,1))*1.0/(199 - pos(i-1,j,0));

            float spd_z = float(rand()%20-10);

            spd(iter_now,j,0) = spd(iter_now-1,j,0)+ int(
                    (p_best_location(i-1,j,0) - pos(i-1,j,0))
                    + (g_best_location(i-1,0) - pos(i-1,j,0)));

            spd(iter_now,j,1) = spd(iter_now-1,j,1)
                + int(
                        0.2*spd_z
                        + 0.5*(p_best_location(i-1,j,1) - pos(i-1,j,1))
                        + 0.5*(g_best_location(i-1,1) - pos(i-1,j,1))
                        + 0.4*zc);


            if(spd(iter_now,j,0) > spdx_max) {spd(iter_now,j,0) = spdx_max;}
            if(spd(iter_now,j,0) < spdx_min) {spd(iter_now,j,0) = spdx_min;}
            if(spd(iter_now,j,1) > spdz_max) {spd(iter_now,j,1) = spdz_max;}
            if(spd(iter_now,j,1) < spdz_min) {spd(iter_now,j,1) = spdz_min;}

            //Location update
            pos(i,j,0) = pos(i-1,j,0) + spd(i-1,j,0);
            pos(i,j,1) = pos(i-1,j,1) + spd(i-1,j,1);

            if(pos(i,j,0) > posx_max) {pos(i,j,0) = posx_max;}
            if(pos(i,j,0) < posx_min) {pos(i,j,0) = posx_min;}
            if(pos(i,j,1) > posz_max) {pos(i,j,1) = posz_max;}
            if(pos(i,j,1) < posz_min) {pos(i,j,1) = posz_min;}
        }
        rew_func(iter_now);
    }
    /*
       int ittemp=0;
       for(int i=0; i<p_num; i++)
       {
       string str_i;
       str_i.resize(20);
       itoa(i+1000,&(str_i.at(0) ) );

       cout<<ittemp<<":"<<str_i.substr(2)<<": [";
       cout<<p_best_location(ittemp,i,0)<<","<<p_best_location(ittemp,i,1)<<"]"<<endl;
       }
       */
}










