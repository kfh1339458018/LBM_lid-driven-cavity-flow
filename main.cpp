#include <stdio.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <iostream>


#define Nx 256
#define Ny 256
#define q 9
#define EndTime 500000
#define WriteInterval 10000

int Re=1000; // 雷诺数


//分布函数的概率分布
static double prob[q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
//D2Q9网格规则(格子的长度为1)
int cx[9]={0,1,0,-1,0,1,-1,-1,1};
int cy[9]={0,0,1,0,-1,1,1,-1,-1};


double func_eq(int direc,double vx,double vy,double rho);

//将数字转换成数组
std::string numberToString(double num)
{
    std::ostringstream ss;
    ss << std::setprecision(6) << num;
    return ss.str();
}

int main() {

        //初始化
        //速度的初始条件定义
        //定义空的速度数组
        static double vx[Nx][Ny];//x方向
        static double vy[Nx][Ny];//y方向
        //定义任意坐标点；
        int x,y;
        int direc;//direc取0~8

        //判断收敛
        bool stop = false;
        static double vx_history[Nx][Ny][1000];
        static double vy_history[Nx][Ny][1000];


        //宏观流体密度初值
        static double rho[Nx][Ny];
        for(x=0;x<Nx;x++){
            for(y=0;y<Ny;y++)
            {
                rho[x][y]=1.0;

            }
        }

        //上边界vx=0.1
        double speed_upwards=0.1;
        for(x=0; x<Ny; x++)
        {
            y=Ny-1;
            vx[x][y]=speed_upwards;
        }


        static double vis;//黏度
        static double tau;//松弛时间
        vis=speed_upwards*Nx/Re;
        tau=0.5+3*vis;

        printf("tau=%lf",tau);

        //概率函数定义
        static double func[q][Nx][Ny];
        //初始条件赋值
        for(x=0;x<Nx;x++)
            for(y=0;y<Ny;y++)
                for(direc=0;direc<q;direc++)
                {
                    func[direc][x][y]=func_eq(direc,vx[x][y],vy[x][y],rho[x][y]);
                }

        //定义迭代次数
        int time;


        //主循环
        for(time = 0; time <= EndTime; time ++) {

            //计算宏观密度
            for (x = 0; x < Nx; x++)
                for (y = 0; y < Ny; y++) {
                    rho[x][y] = 0.0;
                    for (direc = 0; direc < q; direc++) {
                        rho[x][y] += func[direc][x][y];
                    }
                }


            //计算宏观速度
            for (x = 0; x < Nx; x++)
                for (y = 0; y < Ny; y++) {

                    vx_history[x][y][time%1000]=vx[x][y];
                    vy_history[x][y][time%1000]=vy[x][y];

                    vx[x][y] = 0.0;
                    vy[x][y] = 0.0;
                    for (direc = 0; direc < q; direc++) {
                        vx[x][y] += func[direc][x][y] * cx[direc];
                        vy[x][y] += func[direc][x][y] * cy[direc];
                    }

                    vx[x][y] /= rho[x][y];
                    vy[x][y] /= rho[x][y];
                }

            //边界速度重置
            for (x = 0; x < Nx; x++) {
                y = Ny-1;
                vx[x][y] = speed_upwards;

            }

            //计算平衡态
            //平衡态函数
            static double Func_eq[q][Nx][Ny];
            for (x = 0; x < Nx; x++)
                for (y = 0; y < Ny; y++)
                    for (direc = 0; direc < q; direc++) {
                        Func_eq[direc][x][y] = func_eq(direc, vx[x][y], vy[x][y], rho[x][y]);
                    }

            //碰撞
            //碰撞后函数
            static double func_out[q][Nx][Ny];
            for (x = 0; x < Nx; x++)
                for (y = 0; y < Ny; y++)
                    for (direc = 0; direc < q; direc++) {
                        func_out[direc][x][y] = func[direc][x][y] - (func[direc][x][y] - Func_eq[direc][x][y]) / tau;
                    }
            

            //边界条件-非平衡态外推
            //上下
            for(x=0;x<Nx;x++)
                for(direc=0;direc<q;direc++)
                {
                    //上边界
                    y=Ny-1;
                    rho[x][y]=2*rho[x][y-1]-rho[x][y-2];//密度
                    vx[x][y]=speed_upwards;//顶盖速度重置
                    func[direc][x][y]= func_eq(direc,vx[x][y],vy[x][y],rho[x][y])+func[direc][x][y-1]- func_eq(direc,vx[x][y-1],vy[x][y-1],rho[x][y-1]);

                    //下边界
                    y=0;
                    rho[x][y]=2*rho[x][y+1]-rho[x][y+2];//密度
                    func[direc][x][y]= func_eq(direc,vx[x][y],vy[x][y],rho[x][y])+func[direc][x][y+1]- func_eq(direc,vx[x][y+1],vy[x][y+1],rho[x][y+1]);

                }


            //左右
            for (y = 1; y < Ny - 1; y++) {
                for (direc = 0; direc < q; direc++) {
                    //左边界
                    x = 0;
                    rho[x][y] = 2*rho[x+1][y]-rho[x+2][y];
                    func[direc][x][y] = func_eq(direc,vx[x][y],vy[x][y],rho[x][y])+func[direc][x+1][y]-func_eq(direc,vx[x+1][y],vy[x+1][y],rho[x+1][y]);

                    //右边界
                    x = Nx-1;
                    rho[x][y] = 2*rho[x-1][y]-rho[x-2][y];
                    func[direc][x][y] = func_eq(direc,vx[x][y],vy[x][y],rho[x][y])+func[direc][x-1][y]-func_eq(direc,vx[x-1][y],vy[x-1][y],rho[x-1][y]);
                }
            }

            //扩散过程
            for(x=1;x<Nx-1;x++)
                for(y=1;y<Ny-1;y++)
                    for(direc=0;direc<q;direc++)
                    {
                        int x_,y_;
                        x_=x-cx[direc];
                        y_=y-cy[direc];

                        func[direc][x][y]=func_out[direc][x_][y_];

                    }


            //判断是否达到收敛条件
            // blow up
            for(x=0;x<Nx;x++)
                for(y=0;y<Ny;y++)
                {
                    if (isnan(rho[x][y]) || isnan(vx[x][y]) || isnan(vy[x][y]))
                    {
                        std::cout << "blow up point at: x=" << x << "y=" << y << std::endl; 
                        stop=true;
                    }
                    
                }
            //steady judge
            if(time>1000)
            {
                double s1=0.0,s2=0.0,L2_error=0.0;
                for(x=0;x<Nx;x++)
                    for(y=0;y<Ny;y++)
                    {
                        s1+=sqrt(pow(vx[x][y]-vx_history[x][y][time%1000],2)+pow(vy[x][y]-vy_history[x][y][time%1000],2));
                        s2+=sqrt(pow(vx[x][y],2)+pow(vy[x][y],2));
                    }
                L2_error=s1/s2;
                if(L2_error<=1e-6)
                {
                    stop=true;
                }
            }


            //输出
            if((time%WriteInterval == 0 & time!=0) || stop==true)
            {
                std::fstream fopen("result_Re="+ numberToString(Re)+"_step=" + numberToString(time) + ".dat",std::ios::out);
                fopen << "TITLE = \"Temperature Field\"" << std::endl
                      << "VARIABLES = \"X\", \"Y\", \"Rho\", \"U\", \"V\"" << std::endl
                      << "ZONE T=\"My Zone\", I=" << 256 <<", J=" << 256 << ", DATAPACKING=POINT, SOLUTIONTIME=" << time << std::endl;
                double deltaX = 1.0/256.0;
                double deltaY = 1.0/256.0;
                for (y = 0; y < 256; y++)
                {
                    for(x=0; x < 256; x++)
                    {
                        fopen << (x - 0.5)*deltaX << " " << (y - 0.5)*deltaY << " " << rho[x][y] << " " << vx[x][y] << " " << vy[x][y]<< std::endl;
                    }
                }
                std::cout << "current step=" << time << std::endl;
            }
            if(stop==true)
            {
                break;
            }
    }

    return 0;
}

double func_eq(int direc,double vx,double vy,double rho)
{
    double func_eq,a,b;

    //中间量的计算
    a=cx[direc]*vx+cy[direc]*vy;
    b=vx*vx+vy*vy;
    func_eq=prob[direc]*rho*(1.0+3.0*a+4.5*a*a-1.5*b);
    return func_eq;
}



