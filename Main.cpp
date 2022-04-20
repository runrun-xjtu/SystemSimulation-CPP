#include <iostream>
#include <iomanip>
#include <math.h>
#include <windows.h>
#include <string>
#include "Refrigeration_struct.h"

using namespace std;

struct refri point[6];
struct refri lastpoint[6];
struct wetair inside, outside;

/* 程序计算单位制 t-K, p-Pa, d-kg/m3, h-kJ/kg, s-kJ/(kg.K) vis- Pa.s cp-J/(mol.K) th_con-W/(m.K) */
/* REFPROP函数单位制 t-K, p-kPa, d-mol/L, h-J/mol, s-J/(mol.K) vis- uPa.s cp-J/(g.K) th_con-W/(m.K) */

/**********************************************全局变量定义******************************************************/
/*常量参数*/
const double PI = 3.1415926;
const double E_value = 2.71828183;
const double R_value = 8.314;

/*可变输入参数：时间步参数 + 设计参数（温度单位℃）*/
double timestepNum = 10, timestepValue = 0.1;
double tairInside_in = 25, pairInside = 101325, vairInside_in = 2, humidity_in = 0.5;
double tairOutside_in = 35, pairOutside = 101325, vairOutside_in = 3, humidity_out = 0;
double tcon_des = 45, teva_des = 10, overheat_des = 0.1, subcool_des = 5, Q_des = 2500;
double kV_eva, kV_con;
/*结构参数（初值以测试）*/
double V_com = 0.0001;                   //气缸容积，转速
double L_cap = 3, d_cap = 0.0015, epsilon = 0.00005; //毛细管长，内径，粗糙度
double tubenum_con = 5, rownum_con = 1, L_con = 2, D_con = 0.01, d_con = 0.008;//分路数，排数，管长、外径、内径
double finPitch_con = 0.002, finThickness_con = 0.00015, tubeZongPitch_con = 0.022, tubeHengPitch_con = 0.025, Hfin_con = 0.006; //翅片间距，厚度，管纵向间距，横向间距，翅片高度

double tubenum_eva = 5, rownum_eva = 1, L_eva = 1, D_eva = 0.01, d_eva = 0.008;
double finPitch_eva = 0.002, finThickness_eva = 0.00015, tubeZongPitch_eva = 0.022, tubeHengPitch_eva = 0.025, Hfin_eva = 0.006;

double thermalConduct_Wall = 384, thermalConduct_fin = 217;     //管壁热导率，翅片热导率

/*系统输出参数*/
double peva, teva, pcon, tcon;
double P_com, condensationWater, Q_eva;

/*系统计算过程参数*/
double V_con, Vliq_con, V_eva, Vliq_eva;
double pcon_lasttime, tcon_lasttime, pcon_nexttime, tcon_nexttime;
double qm_std, p_std, t_std;
double mSum_air, tSum_air, tAve_air, mFrac_H2O_out;
/****************************************************************************************************************/

/*系统主要函数声明*/

void design(struct refri point[]);
void operation(struct refri point[], struct refri lastpoint[]);
void compressor(struct refri point[], int num);
void condenser(struct refri point[], struct wetair outside, int num);
void capillary(struct refri point[], int num);
void evaporator(struct refri point[], struct wetair inside, int num);

/*主函数*/
int main()
{
	design(point);
	operation(point, lastpoint);
	cout << "仿真结束" << endl;
}

/*部件设计*/
void design(struct refri point[])
{
	//quality_des 优先级大于 overheat_des
	double qm_temp, t_temp, min_esti, max_esti;

	V_com = 0.00001; L_cap = 3; L_con = 2; L_eva = 1; //估计值

	point[0].tem = tcon_des + 273.15;
	point[0].tSatProperty(2);
	point[3].pressure = point[2].pressure = point[0].pressure;
	point[3].tem = tcon_des - subcool_des + 273.15;
	point[3].tpProperty();

	point[0].tem = teva_des + 273.15;
	point[0].tSatProperty(1);
	point[1].pressure = point[5].pressure = point[4].pressure = point[0].pressure;
	point[1].tem = teva_des + overheat_des + 273.15;
	point[1].tpProperty();

	point[4].enthalpy = point[3].enthalpy;
	point[4].phProperty();

	point[2].pressure = point[3].pressure;
	point[2].tem = point[1].tem * pow(point[2].pressure / point[1].pressure, (point[1].k_com - 1) / point[1].k_com);
	point[2].tpProperty();

	qm_std = qm_temp = point[1].qm = point[2].qm = point[3].qm = point[4].qm = Q_des / (1000 * (point[1].enthalpy - point[4].enthalpy));
	p_std = point[1].pressure; t_std = point[1].tem;

	/*压缩机设计： V_com */
	/*
	do
	{
		V_com = 2 * V_com;
		compressor(point, 1);
	} while (qm_temp > point[1].qm);
	min_esti = 0; max_esti = V_com;
	do
	{
		V_com = 0.5 * (min_esti + max_esti);
		compressor(point, 1);
		if (qm_temp < point[1].qm)
		{
			min_esti = min_esti;
			max_esti = V_com;
		}
		else if (qm_temp > point[1].qm)
		{
			min_esti = V_com;
			max_esti = max_esti;
		}
		else break;
	} while (abs(qm_temp - point[1].qm) > 0.00001);
	cout << "压缩机 V_com(m3) = " << V_com << endl;
    */
	V_com = 1.24023e-05;

	/*毛细管设计： L_cap (d_cap = 0.001 - 0.002 e = 0.00005) */
	
	/*
	do
	{
		L_cap = 2 * L_cap;
		capillary(point, 3);
	} while (qm_temp < point[3].qm);
	min_esti = 0; max_esti = L_cap;
	do
	{
		L_cap = 0.5 * (min_esti + max_esti);
		capillary(point, 3);
		if (qm_temp > point[3].qm)
		{
			min_esti = min_esti;
			max_esti = L_cap;
		}
		else if (qm_temp < point[3].qm)
		{
			min_esti = L_cap;
			max_esti = max_esti;
		}
		else break;
	} while (abs(qm_temp - point[3].qm) > 0.0001);
	cout << "毛细管 d_cap(m) = " << setw(8) << left << d_cap << "L_cap(m) = " << L_cap << endl;
    */
    
	d_cap = 0.0015; L_cap = 0.421875;
	/*冷凝器设计： L_con  */
	/*
	outside.humidity = humidity_out;
	outside.tem = tairOutside_in + 273.15;
	outside.pressure = pairOutside;
	t_temp = point[3].tem;
	do
	{
		L_con = 2 * L_con;
		condenser(point, outside, 2);
	} while (t_temp < point[3].tem);
	min_esti = 0; max_esti = L_con;
	do
	{
		L_con = 0.5 * (min_esti + max_esti);
		condenser(point, outside, 2);
		if (t_temp > point[3].tem)
		{
			min_esti = min_esti;
			max_esti = L_con;
		}
		else if (t_temp < point[3].tem)
		{
			min_esti = L_con;
			max_esti = max_esti;
		}
		else break;
	} while (abs(t_temp - point[3].tem) > 0.01);
	cout << "冷凝器 L_con(m) = " << L_con << endl;
	*/

	L_con = 5.23047;

	/*蒸发器设计： L_con  */
	/*
	inside.humidity = humidity_in;
	inside.tem = tairInside_in + 273.15;
	inside.pressure = pairInside;
	t_temp = point[1].tem;
	do
	{
		L_eva = 2 * L_eva;
		evaporator(point, inside, 4);
	} while (t_temp > point[5].tem);
	min_esti = 0; max_esti = L_eva;
	do
	{
		L_eva = 0.5 * (min_esti + max_esti);
		evaporator(point, inside, 4);
		if (t_temp < point[5].tem)
		{
			min_esti = min_esti;
			max_esti = L_eva;
		}
		else if (t_temp > point[5].tem)
		{
			min_esti = L_eva;
			max_esti = max_esti;
		}
		else break;
	} while (abs(t_temp - point[5].tem) > 0.01);
	cout << "蒸发器 L_eva(m) = " << L_eva << endl;
    */

	L_eva = 0.7146;

	point[0].tem = tcon_des + 273.15;
	point[0].tSatProperty(2);
	point[3].pressure = point[2].pressure = point[0].pressure;
	point[3].tem = tcon_des - subcool_des + 273.15;
	point[3].tpProperty();

	point[0].tem = teva_des + 273.15;
	point[0].tSatProperty(1);
	point[1].pressure = point[5].pressure = point[4].pressure = point[0].pressure;
	point[1].tem = teva_des + overheat_des + 273.15;
	point[1].tpProperty();

	point[4].enthalpy = point[3].enthalpy;
	point[4].phProperty();

	point[2].pressure = point[3].pressure;
	point[2].tem = point[1].tem * pow(point[2].pressure / point[1].pressure, (point[1].k_com - 1) / point[1].k_com);
	point[2].tpProperty();

	qm_temp = point[1].qm = point[2].qm = point[3].qm = point[4].qm = Q_des / (1000 * (point[1].enthalpy - point[4].enthalpy));

}

/*系统运行*/
void operation(struct refri point[], struct refri lastpoint[])
{
	struct refri liq, gas;
	struct wetair Outair;
	double lastV_eva, lastV_con;
	double qmH2OAir, qmH2O, dd_out, pressure_h2oout, humidity_in_out;
	kV_eva = 0.0002; kV_con = 0.001; Vliq_con = 0;  
	lastpoint[1] = point[1]; lastpoint[2] = point[2]; lastpoint[3] = point[3]; lastpoint[4] = point[4];
	outside.humidity = humidity_out;
	outside.tem = tairOutside_in + 273.15;
	outside.pressure = pairOutside;
	inside.humidity = humidity_in;
	inside.tem = tairInside_in + 273.15;
	inside.pressure = pairInside;

	inside.tpProperty();
	outside.tpProperty();

	for (int i = 0; i < 30; i++)
	{
		cout << "------------------------------------------" << endl;
		cout << "*第" << i + 1 << "个时间步：" << i * timestepValue << " - " << (i + 1) * timestepValue << endl;

		//if (i > 0) vairInside_in = 0.5 * i + 2;
		//if (i > 15)vairInside_in = 10;

		evaporator(lastpoint, inside, 4); // 输入：上一时间步(4点状态参数 + 流量) + 当前时间步（inside风参数），输出：当前时间步1点流量
		if (lastpoint[5].quality >= 1)     //压缩机流量：蒸发器蒸发量 + 气液分离器蒸发量
		{
			point[5].quality = lastpoint[5].quality = 1;
			point[5].tem = lastpoint[5].tem;
			point[5].pressure = liq.pressure = gas.pressure = lastpoint[5].pressure;
			point[5].tpProperty(); liq.pSatProperty(1); gas.pSatProperty(2);
			point[5].qm = lastpoint[4].qm * lastpoint[5].quality * (1 + (point[5].enthalpy - gas.enthalpy) / (gas.enthalpy - liq.enthalpy));
		}
		else
		{
			point[5].quality = lastpoint[5].quality;
			gas.tem = point[5].tem = lastpoint[5].tem;
			point[5].qm = lastpoint[4].qm * lastpoint[5].quality;
		}

		/*蒸发压力：蒸发量引起蒸发压力变化*/
		V_eva = tubenum_eva * L_eva * 0.25 * PI * d_eva * d_eva - Vliq_eva;
		if (i == 0) lastV_eva = V_eva;

		//point[1].pressure = lastpoint[1].pressure * (1 + (lastV_eva - V_eva) / V_eva) + (point[1].qm * timestepValue - lastpoint[1].qm * timestepValue) / point[1].molmass * R_value * point[1].tem / (V_eva * kV_eva);
		point[1].pressure = 0.5 * (lastpoint[1].pressure + point[5].qm / qm_std * gas.tem / t_std * p_std);
		//point[1].pressure = point[1].qm / qm_std * point[1].tem / t_std * p_std;
		
		point[1].pSatProperty(2);
		peva = point[4].pressure = point[1].pressure + lastpoint[4].pressure - lastpoint[5].pressure;

		/*冷凝压力：1点3点流量差引起冷凝压力变化*/
		V_con = tubenum_con * L_con * 0.25 * PI * d_con * d_con - Vliq_con;
		if (i == 0) lastV_con = V_con;		 
		pcon = point[2].pressure = lastpoint[2].pressure * (1 + (lastV_con - V_con) / V_con) + (lastpoint[1].qm * timestepValue - lastpoint[3].qm * timestepValue) / point[1].molmass * R_value * lastpoint[2].tem / (lastV_con * kV_con);

		compressor(point, 1);            // 输入：1点状态参数 + 2点压力（冷凝压力），输出：2点状态参数
		point[2].tpProperty();
		point[2].qm = point[1].qm;
		
		condenser(point, outside, 2);    // 输入：2点状态参数 + outside风参数，输出：3点状态参数
		point[3].tpProperty();

	  	capillary(point, 3);             // 输入：3点状态参数 + 4点压力，输出：毛细管流量 + 4点状态参数
		point[4].enthalpy = point[3].enthalpy;
		point[4].phProperty();

		/*出风空气*/
		qmH2OAir = inside.density * vairInside_in * (tubeZongPitch_eva * tubenum_eva) * L_eva;
		qmH2O = qmH2OAir * (inside.dd / (inside.dd + 1));
		Outair.tem = tAve_air; Outair.pressure = pairInside; Outair.humidity = 0.1;
		Outair.tpProperty();
		dd_out = (qmH2O - condensationWater) / (qmH2OAir - condensationWater);
		pressure_h2oout = dd_out * pairInside / (0.622 + dd_out);
		humidity_in_out = pressure_h2oout / Outair.psat_h2o;

		if (abs(point[4].qm - point[1].qm) < 0.0001) cout << "系统达到稳态" << endl;
		cout << "Pcon(kPa) = " << setw(12) << left << pcon / 1000;
		cout << "   Peva(kPa) = " << peva / 1000 << endl;
		cout << "Point 1 tem(C) = " << setw(7) << left << point[1].tem - 273.15;
		cout << "   Point 3 tem(C) = " << point[3].tem - 273.15 << endl;
		cout << "Point 1 p(kPa) = " << setw(7) << left << point[1].pressure / 1000;
		cout << "   Point 3 p(kPa) = " << point[3].pressure / 1000 << endl;
		cout << "qm1(kg/s) = " << setw(12) << left << point[1].qm;
		cout << "   qm3(kg/s) = " << point[4].qm << endl;
		cout << "P_com(W) = " << setw(13) << left << P_com;
		cout << "   Cooling capacity(W) = " << Q_eva << endl;
		cout << "COP = " << setw(18) << left << Q_eva / P_com;
		cout << "   condensening water mass(kg/s) = " << condensationWater << endl;
		cout << "inlet air tem(C) = " << setw(5) << left << inside.tem - 273.15;
		cout << "   outlet air tem(C) = " << tAve_air - 273.15 << endl;
		cout << "inlet humidity(%) = " << setw(5) << left << inside.humidity;
		cout << "   outlet humidity(%) = " << humidity_in_out << endl;
		lastV_eva = V_eva; lastV_con = V_con;
		lastpoint[1] = point[1]; lastpoint[2] = point[2]; lastpoint[3] = point[3]; lastpoint[4] = point[4];
	}
}


    /* 冷凝器测试 */
	/*
	outside.humidity = humidity_out;
	point[2].pressure = 1419000;
	point[2].tem = 333;
	point[2].qm = 0.005;
	outside.tem = tairOutside_in + 273.15;
	outside.pressure = 101325;
	condenser(point, outside, 2);
	pconDynamic(point, lastpoint);
	cout << V_con << "  " << Vliq_con << endl;
	*/

	/* 蒸发器测试 */
    /*
	inside.humidity = humidity_in;
	point[4].pressure = 677880;
	point[4].enthalpy = 245;
	point[4].qm = 0.005;
	inside.tem = tairInside_in + 273.15;
	inside.pressure = 101325;

	evaporator(point, inside, 4);
	*/

	/* 毛细管测试 */
	/* 
	point[3].pressure = 1732100;
	point[4].pressure = 1020100;
	point[3].tem = 310;
	capillary(point, 3);
	cout << "出口压力  " << point[4].pressure << endl;
	cout << "流量  " << point[4].qm << endl;
    */

    /* 压缩机测试(输入吸排气压力，蒸发温度)，吸气无过热 */
	/*
	peva = 500000; pcon = 2000000;
	compressor(point, 1);
	*/

	/* REFPROP测试 */
	/*
	cout << "REFPROP测试" << endl;
	double t = 273.3;
	double p = 500000;
	double h = 21000;

	point[1].tem = t;
	point[1].pressure = p;
	point[1].tpProperty();

	cout << point[1].tem << "  温度" << endl;
	cout << point[1].pressure << "  压力" << endl;
	cout << point[1].enthalpy << "  焓" << endl;
	cout << point[1].entropy << "  熵" << endl;
	cout << point[1].quality << "  干度" << endl;
	cout << point[1].density << "  密度" << endl;
	cout << point[1].viscosity << "  粘度" << endl;
	cout << point[1].molmass << "  摩尔质量" << endl;
	*/

	/* 降压REFPROP计算测试 */
	/*
	double xxx;
	point[0].pressure = 1732100;
	point[0].enthalpy = 245.43;
	while (point[0].pressure > 10000)
	{
		point[0].pressure -= 1000;
		point[0].phProperty();
		cout << point[0].pressure << "shang  " << point[0].entropy << endl;
	}
	cin >> xxx;
	*/

	//compressor(point, 1);
	//condenser(point, 2);
	//capillary(point, 3);
	//evaporator(point, 4);

