#include <iostream>
#include <iomanip>
#include <math.h>
#include <windows.h>
#include <string>
#include "Refrigeration_struct.h"

using namespace std;

const double PI = 3.1415926;
const double E_value = 2.71828183;

/*�ⲿ����ϵͳ����*/
extern double peva, teva, pcon, tcon;
extern double P_com;

/*�ⲿ����ṹ����*/
extern double V_com; //�����ݻ���ת��


/*�ⲿ����ģ�麯������*/

extern void compressor(struct refri point[], int num);
extern void condenser(struct refri point[], struct wetair outside, int num);
extern void capillary(struct refri point[], int num);
extern void evaporator(struct refri point[], struct wetair inside, int num);

/*ѹ����ģ��*/
void compressor(struct refri point[], int num)
{
	double nv, nel; //�ݻ�Ч�ʣ���Ч��
	double nv_v, nv_p, nv_t, nv_l; //�ݻ���ѹ�����¶ȣ�й¶ϵ��
	double nel_i, nel_t, nel_l, nel_m, nel_mo; //ָʾ�����ȣ�й¶����е�����Ч��
	double epsilon_p, dp_dm;  //ѹ�ȣ�������ѹ��
	double rpm = 3000;

	epsilon_p = point[num + 1].pressure / point[num].pressure;
	if (point[num].quality >= 1) point[num].tpProperty();
	else point[num].pSatProperty(2);
	
	point[num + 1].tem = point[num].tem * pow(epsilon_p, (point[num].k_com - 1) / point[num].k_com);
	point[num + 1].tpProperty();

	/*�ݻ�Ч������·����*/
	nv_v = 1 - 0.015 * (pow(epsilon_p, 1 / point[num].k_com) - 1);
	nv_p = 1; nv_t = 2.57e-3 * point[num + 1].tem + 1.06e-3 * (point[num].tem - point[num].tem);
	nv_l = 0.98;
	nv = nv_v * nv_p * nv_t * nv_l;
	point[num].qm = nv * rpm * V_com * point[num].density / 60;

	/*��Ч���빦�ļ���*/
	dp_dm = 25 * pow(point[num + 1].tem-273.15, -1.01) * pow(10, -0.15 * epsilon_p);
	nel_i = nv_t * nv_l / (1 + 1.5 * dp_dm * pow(epsilon_p, 1.0 / point[num].k_com) / point[num].density / (point[num + 1].enthalpy - point[num].enthalpy));
	nel_t = nv_t; nel_l = nv_l; nel_m = 0.8; nel_mo = 0.8;
	nel = nel_i * nel_t * nel_l * nel_m * nel_mo;

	P_com = (1 / nel) * point[num].qm * point[num].pressure * (1 / point[num].density) * (point[num].k_com / (point[num].k_com - 1)) * (pow(epsilon_p, (point[num].k_com - 1) / point[num].k_com) - 1);
} 
