#include <iostream>
#include <iomanip>
#include <math.h>
#include <windows.h>
#include <string>
#include "Refrigeration_struct.h"

using namespace std;

const double PI = 3.1415926;
const double E_value = 2.71828183;

/*外部引入系统参数*/
extern double peva, teva, pcon, tcon;
extern double P_com;

/*外部引入结构参数*/
extern double V_com; //气缸容积，转速


/*外部引入模块函数声明*/

extern void compressor(struct refri point[], int num);
extern void condenser(struct refri point[], struct wetair outside, int num);
extern void capillary(struct refri point[], int num);
extern void evaporator(struct refri point[], struct wetair inside, int num);

/*压缩机模块*/
void compressor(struct refri point[], int num)
{
	double nv, nel; //容积效率，电效率
	double nv_v, nv_p, nv_t, nv_l; //容积，压力，温度，泄露系数
	double nel_i, nel_t, nel_l, nel_m, nel_mo; //指示，加热，泄露，机械，电机效率
	double epsilon_p, dp_dm;  //压比，排气阀压降
	double rpm = 3000;

	epsilon_p = point[num + 1].pressure / point[num].pressure;
	if (point[num].quality >= 1) point[num].tpProperty();
	else point[num].pSatProperty(2);
	
	point[num + 1].tem = point[num].tem * pow(epsilon_p, (point[num].k_com - 1) / point[num].k_com);
	point[num + 1].tpProperty();

	/*容积效率与流路计算*/
	nv_v = 1 - 0.015 * (pow(epsilon_p, 1 / point[num].k_com) - 1);
	nv_p = 1; nv_t = 2.57e-3 * point[num + 1].tem + 1.06e-3 * (point[num].tem - point[num].tem);
	nv_l = 0.98;
	nv = nv_v * nv_p * nv_t * nv_l;
	point[num].qm = nv * rpm * V_com * point[num].density / 60;

	/*电效率与功耗计算*/
	dp_dm = 25 * pow(point[num + 1].tem-273.15, -1.01) * pow(10, -0.15 * epsilon_p);
	nel_i = nv_t * nv_l / (1 + 1.5 * dp_dm * pow(epsilon_p, 1.0 / point[num].k_com) / point[num].density / (point[num + 1].enthalpy - point[num].enthalpy));
	nel_t = nv_t; nel_l = nv_l; nel_m = 0.8; nel_mo = 0.8;
	nel = nel_i * nel_t * nel_l * nel_m * nel_mo;

	P_com = (1 / nel) * point[num].qm * point[num].pressure * (1 / point[num].density) * (point[num].k_com / (point[num].k_com - 1)) * (pow(epsilon_p, (point[num].k_com - 1) / point[num].k_com) - 1);
} 
