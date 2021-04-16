#include "pch.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <windows.h>
#include <string>

using namespace std;

const char refrigerant[] = "r22.fld";     //制冷剂

/*制冷剂性质求解结构体*/
typedef struct refri
{
	/*制冷剂状态参数*/
	/* 程序内计算单位 t-K, p-Pa, d-kg/m3, h-kJ/kg, s-kJ/(kg.K) vis- Pa.s cp-J/(mol.K) th_con-W/(m.K) */
	/* REFPROP函数单位 t-K, p-kPa, d-mol/L, h-J/mol, s-J/(mol.K) vis- uPa.s cp-J/(g.K) th_con-W/(m.K) */
	double tem, pressure, p_cri, molmass;                  //温度 压力 临界压力 摩尔质量
	double quality, density, densityLiq, densityGas;       //干度 密度
	double viscosity, viscosityLiq, viscosityGas;		   //动力粘度
	double enthalpy, entropy, k_com;					   //焓 熵 绝热指数
	double HeatCapacity_p, HeatCapacity_v, thermalConduct;  //定压/定容比热容 热导率
	double qm, heatFlux;                                   //流量 热流量

	/*临时参数（REFPROP函数使用）*/
	typedef void(__stdcall* fp_SETUPdllTYPE)(long&, char*, char*, char*, long&, char*, long, long, long, long); fp_SETUPdllTYPE SETUPdll;
	typedef void(__stdcall* fp_SETMIXdllTYPE)(char*, char*, char*, long&, char*, double*, long&, char*, long, long, long, long, long); fp_SETMIXdllTYPE SETMIXdll;
	typedef void(__stdcall* fp_INFOdllTYPE)(long&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&); fp_INFOdllTYPE INFOdll;
	typedef void(__stdcall* fp_PHFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_PHFLSHdllTYPE PHFLSHdll;
	typedef void(__stdcall* fp_PSFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_PSFLSHdllTYPE PSFLSHdll;
	typedef void(__stdcall* fp_TPFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_TPFLSHdllTYPE TPFLSHdll;
	typedef void(__stdcall* fp_SATPdllTYPE)(double&, double*, long&, double&, double&, double&, double*, double*, long&, char*, long); fp_SATPdllTYPE SATPdll;
	typedef void(__stdcall* fp_SATTdllTYPE)(double&, double*, long&, double&, double&, double&, double*, double*, long&, char*, long); fp_SATTdllTYPE SATTdll;
	typedef void(__stdcall* fp_TRNPRPdllTYPE)(double&, double&, double*, double&, double&, long&, char*, long); fp_TRNPRPdllTYPE TRNPRPdll;
	double x[20], xliq[20], xvap[20], f[20]; long ierr;
	char hf[255 * 20], hrf[3 + 1], herr[255 + 1], hfmix[255 + 1];
	double t, p, dl, dv, d, q, e, h, s, cv, cp, w, b, c, eta, tcx;

	/*构造函数（初始化REFPROP模块）*/
	refri()
	{
		double wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas;
		long info_index = 1; long ii = 1;
		HINSTANCE RefpropdllInstance = LoadLibrary(L"REFPRP64.DLL"); //加载DLL

		/*可能用到的REFPROP函数*/
		SETUPdll = (fp_SETUPdllTYPE)GetProcAddress(RefpropdllInstance, "SETUPdll");
		SETMIXdll = (fp_SETMIXdllTYPE)GetProcAddress(RefpropdllInstance, "SETMIXdll");
		INFOdll = (fp_INFOdllTYPE)GetProcAddress(RefpropdllInstance, "INFOdll");
		PHFLSHdll = (fp_PHFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PHFLSHdll");
		PSFLSHdll = (fp_PSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PSFLSHdll");
		TPFLSHdll = (fp_TPFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TPFLSHdll");
		SATPdll = (fp_SATPdllTYPE)GetProcAddress(RefpropdllInstance, "SATPdll");
		SATTdll = (fp_SATTdllTYPE)GetProcAddress(RefpropdllInstance, "SATTdll");
		TRNPRPdll = (fp_TRNPRPdllTYPE)GetProcAddress(RefpropdllInstance, "TRNPRPdll");
		/*暂不启用REFPROP函数*/
		//HSFLSHdll = (fp_HSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "HSFLSHdll");
		//PDFLSHdll = (fp_PDFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PDFLSHdll");
		//PQFLSHdll = (fp_PQFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PQFLSHdll");
		//TDFLSHdll = (fp_TDFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TDFLSHdll");
		//THFLSHdll = (fp_THFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "THFLSHdll");            
		//TQFLSHdll = (fp_TQFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TQFLSHdll");
		//TSFLSHdll = (fp_TSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TSFLSHdll");
		//SATDdll = (fp_SATDdllTYPE)GetProcAddress(RefpropdllInstance, "SATDdll");
		//SATHdll = (fp_SATHdllTYPE)GetProcAddress(RefpropdllInstance, "SATHdll");        
		//SATSdll = (fp_SATSdllTYPE)GetProcAddress(RefpropdllInstance, "SATSdll");      
		//TPRHOdll = (fp_TPRHOdllTYPE)GetProcAddress(RefpropdllInstance, "TPRHOdll");

		strcpy_s(hf, refrigerant);
		strcpy_s(hfmix, "hmx.bnc");
		strcpy_s(hrf, "DEF");
		strcpy_s(herr, "Ok");
		x[0] = 1.0;
		SETUPdll(ii, hf, hfmix, hrf, ierr, herr, 255 * 20, 255, 3, 255);
		INFOdll(info_index, wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas);
		if (ierr != 0) printf("%s\n", herr);
		molmass = wm; p_cri = 1000 * pc;
	}

	/*初始化setup制冷剂求解REFPROP */
	void setuprefri()
	{
		long ii = 1;
		strcpy_s(hf, refrigerant);
		strcpy_s(hfmix, "hmx.bnc");
		strcpy_s(hrf, "DEF");
		strcpy_s(herr, "Ok");
		SETUPdll(ii, hf, hfmix, hrf, ierr, herr, 255 * 20, 255, 3, 255);
	}

	/*根据 p 饱和态求解（i=1 液 i=2 气）*/
	void pSatProperty(long i)
	{
		setuprefri();
		p = 0.001 * pressure;
		SATPdll(p, x, i, t, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.001; else if (i == 2) t += 0.001;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		if (i == 1) density = dl * molmass; else if (i == 2) density = dv * molmass;
		tem = t; pressure = 1000 * p;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass;
		k_com = cp / cv; thermalConduct = tcx;
	}

	/*根据 T 饱和态求解（i=1 液 i=2 气）*/
	void tSatProperty(long i)
	{
		setuprefri();
		t = tem;
		SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.001; else if (i == 2) t += 0.001;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		if (i == 1) density = dl * molmass; else if (i == 2) density = dv * molmass;
		tem = t; pressure = 1000 * p;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass;
		k_com = cp / cv; thermalConduct = tcx;
	}

	/*根据 T p 求解状态参数*/
	void tpProperty()
	{
		setuprefri();
		t = tem; p = 0.001 * pressure;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; density = d * molmass; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass;
		k_com = cp / cv; thermalConduct = tcx;
	}

	/*根据 p h 求解状态参数*/
	void phProperty()
	{
		setuprefri();
		p = 0.001 * pressure; h = molmass * enthalpy;
		PHFLSHdll(p, h, x, t, d, dl, dv, xliq, xvap, q, e, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		density = d * molmass; densityLiq = dl * molmass; densityGas = dv * molmass;
		if (quality > 0) density = densityGas * quality + densityLiq * (1 - quality);
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass; k_com = cp / cv;
	}

	/*根据 p s 求解状态参数*/
	void psProperty()
	{
		setuprefri();
		p = 0.001 * pressure; s = molmass * entropy;
		PSFLSHdll(p, s, x, t, d, dl, dv, xliq, xvap, q, e, h, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		density = d * molmass; densityLiq = dl * molmass; densityGas = dv * molmass;
		if (quality > 0) density = densityGas * quality + densityLiq * (1 - quality);
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass; k_com = cp / cv;
	}
}refri;

/*空气性质求解结构体*/
typedef struct wetair
{
	double tem, pressure, molmass, density;                 //温度 压力 摩尔质量 密度
	double viscosity, enthalpy, entropy, qm;		        //动力粘度 焓 熵 流量
	double HeatCapacity_p, thermalConduct;                  //定压比热容 热导率
	double humidity, dd, delta_dd;                   //湿度 含湿量 凝结量

	double pressure_h2o, psat_h2o, molmass_h2o, density_h2o;
	double viscosity_h2o, enthalpy_h2o, entropy_h2o;
	double HeatCapacity_p_h2o, thermalConduct_h2o;
	/* -2 表示湿空气冷凝后的状态 */
	double tem_2, pressure_2, pressure_h2o_2, psat_h2o_2, density_h2o_2;
	double enthalpy_2, enthalpy_h2o_2, enthalpy_air_2;
	double enthalpy_h2o_2Liq, deltaenthalpy_h2o_2;

	double pressure_air, molmass_air, density_air;
	double viscosity_air, enthalpy_air, entropy_air;
	double HeatCapacity_p_air, thermalConduct_air;

	/*临时参数（REFPROP函数使用）*/
	typedef void(__stdcall* fp_SETUPdllTYPE)(long&, char*, char*, char*, long&, char*, long, long, long, long); fp_SETUPdllTYPE SETUPdll;
	typedef void(__stdcall* fp_SETMIXdllTYPE)(char*, char*, char*, long&, char*, double*, long&, char*, long, long, long, long, long); fp_SETMIXdllTYPE SETMIXdll;
	typedef void(__stdcall* fp_INFOdllTYPE)(long&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&); fp_INFOdllTYPE INFOdll;
	typedef void(__stdcall* fp_TPFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_TPFLSHdllTYPE TPFLSHdll;
	typedef void(__stdcall* fp_TRNPRPdllTYPE)(double&, double&, double*, double&, double&, long&, char*, long); fp_TRNPRPdllTYPE TRNPRPdll;
	typedef void(__stdcall* fp_SATTdllTYPE)(double&, double*, long&, double&, double&, double&, double*, double*, long&, char*, long); fp_SATTdllTYPE SATTdll;
	double x[20], xliq[20], xvap[20], f[20]; long ierr;
	char hf[255 * 20], hrf[3 + 1], herr[255 + 1], hfmix[255 + 1];
	double t, p, dl, dv, d, q, e, h, s, cv, cp, w, b, c, eta, tcx;
	double wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas;
	long info_index = 1; long ii = 1;

	/*构造函数（初始化REFPROP模块）*/
	wetair()
	{
		HINSTANCE RefpropdllInstance = LoadLibrary(L"REFPRP64.DLL"); //加载DLL
		SETUPdll = (fp_SETUPdllTYPE)GetProcAddress(RefpropdllInstance, "SETUPdll");
		SATTdll = (fp_SATTdllTYPE)GetProcAddress(RefpropdllInstance, "SATTdll");
		//SETMIXdll = (fp_SETMIXdllTYPE)GetProcAddress(RefpropdllInstance, "SETMIXdll");
		INFOdll = (fp_INFOdllTYPE)GetProcAddress(RefpropdllInstance, "INFOdll");
		TPFLSHdll = (fp_TPFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TPFLSHdll");
		TRNPRPdll = (fp_TRNPRPdllTYPE)GetProcAddress(RefpropdllInstance, "TRNPRPdll");
	}
	void setupair()
	{
		strcpy_s(hf, "air.ppf");
		//strcpy_s(hf, "air.mix");
		strcpy_s(hfmix, "hmx.bnc");
		strcpy_s(hrf, "DEF");
		strcpy_s(herr, "Ok");
		SETUPdll(ii, hf, hfmix, hrf, ierr, herr, 255 * 20, 255, 3, 255);
		//SETMIXdll = (fp_SETMIXdllTYPE)GetProcAddress(RefpropdllInstance, "SETMIXdll");
		INFOdll(info_index, wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas);
		if (ierr != 0) printf("%s\n", herr);
		molmass_air = wm;
	}
	void setupwater()
	{
		strcpy_s(hf, "water.fld");
		strcpy_s(hfmix, "hmx.bnc");
		strcpy_s(hrf, "DEF");
		strcpy_s(herr, "Ok");
		SETUPdll(ii, hf, hfmix, hrf, ierr, herr, 255 * 20, 255, 3, 255);
		INFOdll(info_index, wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas);
		if (ierr != 0) printf("%s\n", herr);
		molmass_h2o = wm;
	}

	void tpProperty()
	{
		long i;
		setupwater();
		t = tem; i = 2;
		SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.001; else if (i == 2) t += 0.001;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		if (i == 1) density_h2o = dl * molmass_h2o; else if (i == 2) density_h2o = dv * molmass_h2o;
		tem = t; psat_h2o = 1000 * p;
		enthalpy_h2o = h / molmass_h2o; entropy_h2o = s / molmass_h2o; viscosity_h2o = 1e-6 * eta;
		HeatCapacity_p_h2o = cp / molmass_h2o; thermalConduct_h2o = tcx;

		pressure_h2o = humidity * psat_h2o; //空气压力等于大气压减去水蒸气分压力

		setupair();
		t = tem; p = 0.001 * (pressure - pressure_h2o);
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure_air = 1000 * p; density_air = d * molmass_air;
		enthalpy_air = h / molmass_air; entropy_air = s / molmass_air; viscosity_air = 1e-6 * eta;
		HeatCapacity_p_air = cp / molmass_air; thermalConduct_air = tcx;

		dd = 0.622 * (pressure_h2o / pressure_air); //含湿量计算
		density = density_air + density_h2o;
		enthalpy = (enthalpy_air + dd * enthalpy_h2o) / (1 + dd);
		entropy = (entropy_air + dd * entropy_h2o) / (1 + dd);
		viscosity = (viscosity_air + dd * viscosity_h2o) / (1 + dd);
		HeatCapacity_p = (HeatCapacity_p_air + dd * HeatCapacity_p_h2o) / (1 + dd);
		thermalConduct = (thermalConduct_air + dd * thermalConduct_h2o) / (1 + dd);
	}

	void condensation()
	{
		long i;
		double ddmax, psat_h2o_2;
		setupwater();
		t = tem_2; i = 2;
		SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.001; else if (i == 2) t += 0.001;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		psat_h2o_2 = 1000 * p; enthalpy_h2o_2 = h / molmass_h2o;

		ddmax = 0.622 * psat_h2o_2 / (pressure_2 - psat_h2o_2);
		if (ddmax < dd)
		{
			t = tem_2; i = 1;
			SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, 255);
			if (i == 1) t -= 0.001; else if (i == 2) t += 0.001;
			TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
			enthalpy_h2o_2Liq = h / molmass_h2o;
			deltaenthalpy_h2o_2 = enthalpy_h2o_2 - enthalpy_h2o_2Liq;

			pressure_h2o_2 = psat_h2o_2;
			delta_dd = dd - ddmax;
		}
		else
		{
			pressure_h2o_2 = pressure_h2o;
			delta_dd = 0;
		}

		setupair();
		t = tem_2; p = 0.001 * (pressure_2 - pressure_h2o_2);
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		enthalpy_air_2 = h / molmass_air;
		enthalpy_2 = (enthalpy_air_2 + enthalpy_h2o_2 * (dd - delta_dd)) / (1 + dd - delta_dd);
	}

}wetair;

/*
#include "pch.h"
#include<iostream>
using namespace std;

void CppDynamicLinkLibrary()
{
	cout << "Hello Fluent from cpp dll" << endl;
}
*/

/* 程序计算单位制 t-K, p-Pa, d-kg/m3, h-kJ/kg, s-kJ/(kg.K) vis- Pa.s cp-J/(mol.K) th_con-W/(m.K) */
/* REFPROP函数单位制 t-K, p-kPa, d-mol/L, h-J/mol, s-J/(mol.K) vis- uPa.s cp-J/(g.K) th_con-W/(m.K) */

/**********************************************全局变量定义******************************************************/

/*常量参数*/
const double PI = 3.1415926;
const double E_value = 2.71828183;
const double R_value = 8.314;

/*可变输入参数：时间步参数 + 设计参数（温度单位℃）*/

double tairInside_in = 25, pairInside = 101325, vairInside_in = 2, humidityAir_in = 0.8;
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
double t_space[6], p_space[6], x_space[6], h_space[6],qm_space[6];  //存储lastpoint空间
int times = 0;                                                      //调用dll次数
double timestepValue;                                               //时间步长
double V_con, Vliq_con, V_eva, Vliq_eva;                            //压强计算容积
double pcon_lasttime, tcon_lasttime, pcon_nexttime, tcon_nexttime;  
double qm_std, p_std, t_std;                                        //设计工况标准参数
double heatx_lastmin = 0, heatx_lastmax = 0;                        //加快换热器迭代计算的参数
double lastV_eva = 0, lastV_con = 0;                                //上次计算的容积（用于压强计算）
double mSum_air, tSum_air, tAve_air;                             //用于换热器出风温度、相对湿度计算

/****************************************************************************************************************/

/*系统主要函数声明*/

void design(struct refri point[], struct wetair inside, struct wetair outside);
void operation(struct refri point[], struct refri lastpoint[], struct wetair inside, struct wetair outside, double flowtime, int times);
void compressor(struct refri point[], int num);
void condenser(struct refri point[], struct wetair outside, int num);
void capillary(struct refri point[], int num);
void evaporator(struct refri point[], struct wetair inside, int num);
double heatxcon(struct refri unit, struct wetair outside, double k, double A, double unitLength);
double heatx(struct refri unit, struct wetair inside, double k, double A, double unitLength, double dpair, double* conWaterUnit, double* th_out_return);

/*主函数*/

void mainDll(double temAir_in, double VAir, double mFrac_H2O_in, double flowtime, double timeInterval, double* temAir_out, double* mFrac_H2O_out, double* P)
{
	struct refri point[6];
	struct refri lastpoint[6];
	struct wetair inside, outside;

	//double temAir_in, VAir, humidityAir_in, flowtime, timeInterval;
	//temAir_in = 298.15; VAir = 2; humidityAir_in = 0.8; flowtime = 0; timeInterval = 0.1;

	/*求解室内侧相对湿度 风速 风质量流量 空气状态*/
	double humidityEstimin = 0, humidityEstimax = 1;
	inside.tem = temAir_in;
	inside.pressure = pairInside;
	do
	{
		humidityAir_in = 0.5 * (humidityEstimin + humidityEstimax);
		inside.humidity = humidityAir_in;
		inside.tpProperty();

		if ((inside.dd / (inside.dd + 1) - mFrac_H2O_in) > 0)
		{
			humidityEstimin = humidityEstimin;
			humidityEstimax = humidityAir_in;
		}
		else if ((inside.dd / (inside.dd + 1) - mFrac_H2O_in) < 0)
		{
			humidityEstimin = humidityAir_in;
			humidityEstimax = humidityEstimax;
		}

	} while (abs(inside.dd / (inside.dd + 1) - mFrac_H2O_in) > 0.00001);
	
	inside.tpProperty();

	vairInside_in = VAir / (L_eva * tubeZongPitch_eva * tubenum_eva);
	mSum_air = inside.density * VAir;

	/*求解室外侧空气状态*/
	outside.tem = tairOutside_in + 273.15;
	outside.pressure = pairOutside;
	outside.humidity = humidity_out;
	outside.tpProperty();

	/*开始*/
	timestepValue = timeInterval;
	cout << "------------------------------------------" << endl;
	cout << "Flow time: " << flowtime << " (s) - " << flowtime + timeInterval << " (s)" << endl;

	if (times == 0)
	{
		/*第一次调用 DLL 初始化*/
		design(point, inside, outside);    // 模块的尺寸设计

		/*状态点参数存储到全局变量空间内 供下一次调用DLL使用*/
		for (int i = 0; i < 6; i++)
		{
			if (i == 1)
			{
				if (point[i].quality > 0 && point[i].quality < 1)
				{
					p_space[i] = point[i].pressure;
					h_space[i] = point[i].enthalpy;
					x_space[i] = point[i].quality;
					qm_space[i] = point[i].qm;
				}
				else
				{
					t_space[i] = point[i].tem;
					p_space[i] = point[i].pressure;
					qm_space[i] = point[i].qm;
				}
			}
			else if (i == 2)
			{
				t_space[i] = point[i].tem;
				p_space[i] = point[i].pressure;
				qm_space[i] = point[i].qm;
			}
			else if (i == 3)
			{
				p_space[i] = point[i].pressure;
				h_space[i] = point[i].enthalpy;
				qm_space[i] = point[i].qm;
			}
			else if (i == 4)
			{
				p_space[i] = point[i].pressure;
				h_space[i] = point[i].enthalpy;
				qm_space[i] = point[i].qm;
			}
		}
		times++;
	}
	else
	{
		/*非第一次调用 DLL 加载储存的状态点参数 */
		for (int i = 0; i < 6; i++)
		{
			if (i == 1)
			{
				if (x_space[i] > 0 && x_space[i] < 1)
				{
					lastpoint[i].pressure = p_space[i];
					lastpoint[i].enthalpy = h_space[i];
					lastpoint[i].qm = qm_space[i];
					lastpoint[i].phProperty();
				}
				else
				{
					lastpoint[i].tem = t_space[i];
					lastpoint[i].pressure = p_space[i];
					lastpoint[i].qm = qm_space[i];
					lastpoint[i].tpProperty();
				}
			}
			else if (i == 2)
			{
				lastpoint[i].tem = t_space[i];
				lastpoint[i].pressure = p_space[i];
				lastpoint[i].qm = qm_space[i];
				lastpoint[i].tpProperty();
			}
			else if (i == 3)
			{
				lastpoint[i].pressure = p_space[i];
				lastpoint[i].enthalpy = h_space[i];
				lastpoint[i].qm = qm_space[i];
				lastpoint[i].phProperty();
			}
			else if (i == 4)
			{
				lastpoint[i].pressure = p_space[i];
				lastpoint[i].enthalpy = h_space[i];
				lastpoint[i].qm = qm_space[i];
				lastpoint[i].phProperty();
			}
		}

		/*仿真流程主函数 输出指针变量*/
		operation(point, lastpoint, inside, outside, flowtime, times);
		*P = P_com; *temAir_out = tAve_air; *mFrac_H2O_out = (mSum_air * mFrac_H2O_in - condensationWater) / mSum_air;

		lastV_eva = V_eva; lastV_con = V_con;
		times++;

		/*状态点参数存储到全局变量空间内 供下一次调用DLL使用*/
		for (int i = 0; i < 6; i++)
		{
			if (i == 1)
			{
				if (point[i].quality > 0 && point[i].quality < 1)
				{
					p_space[i] = point[i].pressure;
					h_space[i] = point[i].enthalpy;
					x_space[i] = point[i].quality;
					qm_space[i] = point[i].qm;
				}
				else
				{
					t_space[i] = point[i].tem;
					p_space[i] = point[i].pressure;
					qm_space[i] = point[i].qm;
				}
			}
			else if (i == 2)
			{
				t_space[i] = point[i].tem;
				p_space[i] = point[i].pressure;
				qm_space[i] = point[i].qm;
			}
			else if (i == 3)
			{
				p_space[i] = point[i].pressure;
				h_space[i] = point[i].enthalpy;
				qm_space[i] = point[i].qm;
			}
			else if (i == 4)
			{
				p_space[i] = point[i].pressure;
				h_space[i] = point[i].enthalpy;
				qm_space[i] = point[i].qm;
			}
		}
	}
}

/*部件设计*/
void design(struct refri point[], struct wetair inside, struct wetair outside)
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
void operation(struct refri point[], struct refri lastpoint[], struct wetair inside, struct wetair outside, double flowtime, int times)
{
	struct refri liq, gas;
	kV_eva = 0.0002; kV_con = 0.001; Vliq_con = 0;

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
	if (times == 0 || times == 1) lastV_eva = V_eva;

	//point[1].pressure = lastpoint[1].pressure * (1 + (lastV_eva - V_eva) / V_eva) + (point[1].qm * timestepValue - lastpoint[1].qm * timestepValue) / point[1].molmass * R_value * point[1].tem / (V_eva * kV_eva);
	point[1].pressure = 0.5 * (lastpoint[1].pressure + point[5].qm / qm_std * gas.tem / t_std * p_std);
	//point[1].pressure = point[1].qm / qm_std * point[1].tem / t_std * p_std;

	point[1].pSatProperty(2);
	peva = point[4].pressure = point[1].pressure + lastpoint[4].pressure - lastpoint[5].pressure;

	/*冷凝压力：1点3点流量差引起冷凝压力变化*/
	V_con = tubenum_con * L_con * 0.25 * PI * d_con * d_con - Vliq_con;
	if (times == 0 || times == 1) lastV_con = V_con;
	pcon = point[2].pressure = lastpoint[2].pressure * (1 + (lastV_con - V_con) / V_con) + (lastpoint[1].qm * timestepValue - lastpoint[3].qm * timestepValue) / point[1].molmass * R_value * lastpoint[2].tem / (lastV_con * kV_con);

	compressor(point, 1);            // 输入：1点状态参数 + 2点压力（冷凝压力），输出：2点状态参数
	point[2].tpProperty();
	point[2].qm = point[1].qm;

	condenser(point, outside, 2);    // 输入：2点状态参数 + outside风参数，输出：3点状态参数
	point[3].tpProperty();

	capillary(point, 3);             // 输入：3点状态参数 + 4点压力，输出：毛细管流量 + 4点状态参数
	point[4].enthalpy = point[3].enthalpy;
	point[4].phProperty();

	//if (abs(point[4].qm - point[1].qm) < 0.0001) cout << "系统达到稳态" << endl;
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
}

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
	dp_dm = 25 * pow(point[num + 1].tem - 273.15, -1.01) * pow(10, -0.15 * epsilon_p);
	nel_i = nv_t * nv_l / (1 + 1.5 * dp_dm * pow(epsilon_p, 1.0 / point[num].k_com) / point[num].density / (point[num + 1].enthalpy - point[num].enthalpy));
	nel_t = nv_t; nel_l = nv_l; nel_m = 0.8; nel_mo = 0.8;
	nel = nel_i * nel_t * nel_l * nel_m * nel_mo;

	P_com = (1 / nel) * point[num].qm * point[num].pressure * (1 / point[num].density) * (point[num].k_com / (point[num].k_com - 1)) * (pow(epsilon_p, (point[num].k_com - 1) / point[num].k_com) - 1);
}

/*冷凝器模块*/
void condenser(struct refri point[], struct wetair outside, int num)
{
	int i, j, k, ii;
	double A0_unit, A0_total, A_tube, unitnum, unitLength; //内外侧面积，管内截面积
	double a0, a1, k_total, Qunit, Qtube, Qsum;     //内外和总传热系数

	double j_heatx, f_heatx, part1, part2, part3, part4, P1, P2, P3, P4, P5, P6, F1, F2, F3; //空气侧传热压降计算
	double Dc, Dh, ReD, Pr_fin; double dpair, conWaterUnit;                                  //空气侧传热系数计算
	double nf, n0, n0Arf, Ar, Af, mH, th_mH;     //肋效率计算

	struct refri lastUnit, Unit, newUnit, Liquid, Gas;                                      //管内沸腾传热系数计算
	double v_liq, Re_liq, Pr_liq, Nul_liq, G_liq, pr, Nu, xtt;                    //管内沸腾传热系数计算
	double v_gas, Re_gas, Pr_gas, Nu_gas, G_gas;                         //管内沸腾传热系数计算
	double dpunit, dpfri, f, fl, fg, factorL, factorG, dp_liq, dp_gas;   //管内压降计算

	unitnum = 100; Vliq_con = 0;
	Qunit = Qtube = Qsum = 0;

	/*基本参数计算*/
	Unit = point[num]; Unit.tpProperty();
	Unit.qm = point[num].qm / tubenum_con;    //平均分配
	unitLength = L_con / unitnum;
	A0_unit = PI * d_con * unitLength;
	A0_total = PI * d_con * L_con * tubenum_con;
	A_tube = 0.25 * PI * d_con * d_con;

	/*空气侧传热系数计算*/
	outside.tpProperty();
	Dc = D_con + 2 * finThickness_con;
	Dh = 4 * (2 * Hfin_con) * Dc / D_con;
	ReD = outside.density * vairOutside_in * Dc / outside.viscosity;
	if (rownum_con == 1)
	{
		P1 = 1.9 - 0.23 * log(ReD);
		P2 = -0.236 + 0.126 * log(ReD);
		part1 = pow(tubeHengPitch_con / tubeZongPitch_con, P1);
		part2 = pow(finPitch_con / Dc, -1.084);
		part3 = pow(finPitch_con / Dh, -0.786);
		part4 = pow(finPitch_con / tubeHengPitch_con, P2);
		j_heatx = 0.108 * pow(ReD, -0.29) * part1 * part2 * part3 * part4;
	}
	else
	{
		P3 = -0.361 - 0.042 * rownum_con / log(ReD) + 0.158 * log(rownum_con * pow(finPitch_con / Dc, 0.41));
		P4 = -1.224 - 0.076 * pow(tubeZongPitch_con / Dh, 1.42) / log(ReD);
		P5 = -0.083 + 0.058 * rownum_con / log(ReD);
		P6 = -5.735 + 1.21 * log(ReD / rownum_con);
		part1 = 0.086 * pow(ReD, P3) * pow(rownum_con, P4);
		part2 = pow(finPitch_con / Dc, P5);
		part3 = pow(finPitch_con / Dh, P6);
		part4 = pow(finPitch_con / tubeHengPitch_con, -0.93);
		j_heatx = part1 * part2 * part3 * part4;
	}
	Pr_fin = 1000 * outside.HeatCapacity_p * outside.viscosity / outside.thermalConduct;
	a1 = j_heatx * outside.density * vairOutside_in * (1000 * outside.HeatCapacity_p) / pow(Pr_fin, 2.0 / 3);

	/*肋效率计算*/
	mH = pow((2 * a1 / thermalConduct_fin / finThickness_con), 0.5) * Hfin_con;
	th_mH = (pow(E_value, mH) - pow(E_value, -1 * mH)) / (pow(E_value, mH) + pow(E_value, -1 * mH));
	nf = th_mH / mH;
	Ar = PI * D_con * L_con * tubenum_con;
	Af = 2 * (2 * Hfin_con) * (tubeZongPitch_con * (tubenum_con + 1)) * (int(L_con / finPitch_con) + 1) - 2 * 0.25 * PI * D_con * D_con * tubenum_con * rownum_con;
	n0 = (Ar + nf * Af) / (Ar + Af);
	n0Arf = (Ar + Af) * n0;

	for (ii = 0; ii < unitnum; ii++)
	{
		/*制冷剂侧传热系数与压降计算*/
		if (Unit.quality < 1 && Unit.quality > 0) // 干度判断， x < 1 时计算管内沸腾传热系数
		{
			Liquid = Gas = Unit;
			Liquid.pSatProperty(1);
			Gas.pSatProperty(2);
			v_liq = Liquid.qm / Liquid.density / A_tube;
			Re_liq = Liquid.density * v_liq * d_con / Liquid.viscosity;
			Pr_liq = (1000 * Liquid.HeatCapacity_p) * Liquid.viscosity / Liquid.thermalConduct;
			Nul_liq = 0.023 * pow(Re_liq, 0.8) * pow(Pr_liq, 0.4);
			pr = Unit.pressure / Unit.p_cri;
			Nu = Nul_liq * (pow(1 - Unit.quality, 0.8) + 3.8 * pow(Unit.quality, 0.76) * pow(1 - Unit.quality, 0.04) / pow(pr, 0.38));
			a0 = Nu * Liquid.thermalConduct / d_con;
			xtt = pow(((1 - Unit.quality) / Unit.quality), 0.9) * pow((Gas.density / Liquid.density), 0.5) * pow((Liquid.viscosity / Gas.viscosity), 0.1);

			/*压降计算*/
			factorL = factorG = 0;
			if (Re_liq > 4000) factorL = 1 + 20 / xtt + 1 / (xtt * xtt);
			else factorG = 1 + 20 * xtt + xtt * xtt;
			G_gas = Gas.qm * Gas.quality / A_tube;
			Re_gas = G_gas * d_con / Gas.viscosity;
			G_liq = Liquid.qm * (1 - Liquid.quality) / A_tube;
			if (Re_liq < 2000) fl = 64 / Re_liq;
			else fl = 0.3164 * pow(Re_liq, -0.25);
			if (Re_gas < 2000) fg = 64 / Re_gas;
			else fg = 0.3164 * pow(Re_gas, -0.25);
			dp_liq = fl * unitLength * G_liq * G_liq / (2 * d_con * Liquid.density);
			dp_gas = fg * unitLength * G_gas * G_gas / (2 * d_con * Gas.density);
			dpfri = factorL * dp_liq + factorG * dp_gas; //其中一个因子为0
			dpunit = dpfri + pow((G_gas + G_liq), 2) * (1 / Unit.density - 1 / lastUnit.density);
		}
		else // 计算管内湍流强制对流传热系数（干度 x = 1 时）
		{
			Unit.tpProperty();
			v_gas = Unit.qm / Unit.density / (0.25 * PI * d_con * d_con);
			Re_gas = Unit.density * v_gas * d_con / Unit.viscosity;
			Pr_gas = (1000 * Unit.HeatCapacity_p) * Unit.viscosity / Unit.thermalConduct;
			Nu_gas = 0.023 * pow(Re_gas, 0.8) * pow(Pr_gas, 0.3);
			a0 = Nu_gas * Unit.thermalConduct / d_con;
			/*压降计算*/
			if (Re_gas < 2000) f = 64 / Re_gas;
			else f = 0.3164 * pow(Re_gas, -0.25);
			dpunit = f * unitLength * v_gas * v_gas * Unit.density / (2 * d_con);
		}
		/*单元传热量计算*/
		k_total = 1.0 / ((1 / a0) + (0.5 * (D_con - d_con) / thermalConduct_Wall) + A0_total / (n0Arf * a1));
		Qunit = heatxcon(Unit, outside, k_total, A0_unit, unitLength);
		Qtube += Qunit;
		if (Unit.quality >= 1) Unit.quality = 1;
		else if (Unit.quality < 0) Unit.quality = 0;
		Vliq_con += unitLength * A_tube * (1 - Unit.quality) / Liquid.density * (Unit.quality / Gas.density + (1 - Unit.quality) / Liquid.density);

		lastUnit = Unit;
		newUnit = Unit;
		newUnit.pressure = Unit.pressure - dpunit;
		newUnit.enthalpy = Unit.enthalpy - Qunit / 1000 / Unit.qm;
		newUnit.phProperty();
		Unit = newUnit;
	}
	/*计算换热量与出口参数*/
	point[num + 1] = newUnit;
	Qsum = Qtube * tubenum_con;
	Vliq_con = Vliq_con * tubenum_con;
}

/*冷凝器换热器模块*/
double heatxcon(struct refri unit, struct wetair outside, double k, double A, double unitLength)
{
	double th_in, th_out, tc_in, tc_out, dtm, tc_out_reslut;
	double t_esti, min_esti, max_esti, Q_esti, Q_comp, Q_result;  // t_esti 是 tc_out
	int n = 0;
	unit.phProperty();
	outside.tpProperty();
	outside.qm = outside.density * vairOutside_in * tubeZongPitch_con * unitLength;
	th_in = unit.tem;
	tc_in = outside.tem;
	min_esti = tc_in; max_esti = th_in;

	while (n < 20)
	{
		t_esti = 0.5 * (min_esti + max_esti);
		Q_esti = 1000 * outside.qm * outside.HeatCapacity_p * (t_esti - tc_in);
		if (unit.quality < 1 && unit.quality >= 0)
		{
			th_out = th_in;
		}
		else
		{
			th_out = th_in - Q_esti / 1000 / unit.qm / unit.HeatCapacity_p;
			if (th_out <= tc_in)
			{
				min_esti = min_esti;
				max_esti = t_esti;
				continue;
			}
		}

		dtm = ((th_in - t_esti) - (th_out - tc_in)) / log((th_in - t_esti) / (th_out - tc_in));
		Q_comp = k * A * dtm;

		if (Q_esti >= Q_comp)
		{
			min_esti = min_esti;
			max_esti = t_esti;
		}
		else
		{
			min_esti = t_esti;
			max_esti = max_esti;
		}
		n++;
		if (abs(Q_esti - Q_comp) < 0.001) break;
	}
	Q_result = Q_esti;
	tc_out_reslut = t_esti;
	return Q_result;
}

void capillary(struct refri point[], int num)
{
	struct refri unit;            //毛细管计算结构体
	struct refri unitLast;        //
	struct refri unitLastLast;        //
	struct refri origin;        //
	struct refri soundVelocity, S_liq, S_gas;
	struct refri Liquid, Gas;
	double qmEsti, lengthUnit, v_Sound;      //估计流量，单元长度，音速
	double velocity, G, Re, f, A, B, dpVelcocity, dpFriction; //速度，质流密度
	int n = 100, i, j;                          //微元数
	int errorID = 0;

	qmEsti = 0;
	if (point[num].quality <= 0) point[num].tpProperty();       //初始入口物性计算
	else point[num].phProperty();
	origin = point[num];
	lengthUnit = L_cap / n;

	for (i = 0; i < 6; i++)    //6位精度计算流量
	{
		for (j = 0; j < 10; j++)
		{
			unit = unitLast = unitLastLast = origin;
			qmEsti += pow(0.1, 2 + i);
			for (int k = 0; k < n; k++)
			{
				if (unitLast.quality > 0)
				{
					Liquid = Gas = unitLast;
					Liquid.pSatProperty(1); Gas.pSatProperty(2);
				}
				if (unitLast.quality > 0)
				{
					unitLast.density = unitLast.quality * Gas.density + (1 - unitLast.quality) * Liquid.density;
					unitLast.viscosity = unitLast.quality * Gas.viscosity + (1 - unitLast.quality) * Liquid.viscosity;
				}
				/*速度、质流密度、雷诺数计算*/
				velocity = qmEsti / unitLast.density / (0.25 * PI * d_cap * d_cap);
				G = velocity * unitLast.density;
				Re = unitLast.density * velocity * d_cap / unitLast.viscosity;

				/*压降计算*/
				A = pow(2.457 * log(1 / (pow(7 / Re, 0.9) + 0.27 * epsilon / d_cap)), 16);
				B = pow(37560 / Re, 16);
				f = 8 * pow(pow(8 / Re, 12) + 1.0 / pow(A + B, 1.5), 1.0 / 12);
				dpVelcocity = G * G * (1 / unitLast.density - 1 / unitLastLast.density);
				dpFriction = f * velocity * velocity * lengthUnit * unitLast.density / (2 * d_cap);

				/*更新unit数据（ph函数）*/
				unit.enthalpy = unitLast.enthalpy;
				unit.pressure -= dpVelcocity + dpFriction;
				if (unit.pressure <= point[num + 1].pressure) unit.pressure = point[num + 1].pressure - 10;
				unit.phProperty();
				unitLastLast = unitLast;
				unitLast = unit;
			}
			/*干度大于0时音速计算，判断是否达到音速*/
			if (unit.quality > 0)
			{
				soundVelocity.entropy = unit.entropy;
				soundVelocity.pressure = unit.pressure + 50; //取压力dp=50,求d(density)
				soundVelocity.psProperty();
				S_liq = S_gas = soundVelocity;
				S_liq.pSatProperty(1); S_gas.pSatProperty(2);

				unit.density = unit.quality * Gas.density + (1 - unit.quality) * Liquid.density;
				soundVelocity.density = soundVelocity.quality * S_gas.density + (1 - soundVelocity.quality) * S_liq.density;
				v_Sound = pow(abs(soundVelocity.pressure - unit.pressure) / abs(soundVelocity.density - unit.density), 0.5);
				if (velocity > v_Sound)
				{
					qmEsti -= pow(0.1, 2 + i); break;
				}
			}
			/*出口压力计算，判断是否达到背压*/
			if (unit.pressure <= point[num + 1].pressure)
			{
				qmEsti -= pow(0.1, 2 + i);
				break;
			}
		}
		if (i == 0 && j == 10)
		{
			cout << "超出计算流量范围（qm > 0.1 kg/s）" << endl;
			errorID = 2; break;
		}
	}
	if (errorID == 0)
	{
		point[num + 1].qm = point[num].qm = qmEsti;  //输出流量
	}
}



/*蒸发器模块*/
void evaporator(struct refri point[], struct wetair inside, int num)
{
	int i, j, k, ii;
	double A0_unit, A0_total, A_tube, unitnum, unitLength; //内外侧面积，管内截面积
	double a0, a1, k_total, Qunit, Qtube;     //内外和总传热系数

	double j_heatx, f_heatx, part1, part2, part3, part4, P1, P2, P3, P4, P5, P6, F1, F2, F3; //空气侧传热压降计算
	double Dc, Dh, ReD, Pr_fin; double dpair, conWaterUnit;                                  //空气侧传热系数计算
	double nf, n0, n0Arf, Ar, Af, mH, th_mH;     //肋效率计算

	struct refri lastUnit, Unit, newUnit, Liquid, Gas;                   //管内沸腾传热系数计算
	double v_liq, Re_liq, Pr_liq, Nul_liq, G_liq, al;                    //管内沸腾传热系数计算
	double a0_esti, a0_last, atp, anb, E, S, xtt, pr, Bo, hflux_d, Re;   //管内沸腾传热系数计算
	double v_gas, Re_gas, Pr_gas, Nu_gas, G_gas, th_out;                 //管内沸腾传热系数计算
	double dpunit, dpfri, f, fl, fg, factorL, factorG, dp_liq, dp_gas;   //管内压降计算

	unitnum = 100; V_eva = 0; a0_last = 2000;
	Qunit = Qtube = Q_eva = 0;
	condensationWater = 0; tSum_air = 0;

	/*基本参数计算*/
	Unit = point[num]; Unit.phProperty();
	Unit.qm = point[num].qm / tubenum_eva;    //平均分配
	unitLength = L_eva / unitnum;
	A0_unit = PI * d_eva * unitLength;
	A0_total = PI * d_eva * L_eva * tubenum_eva;

	/*空气侧传热系数计算*/
	inside.tpProperty();
	Dc = D_eva + 2 * finThickness_eva;
	Dh = 4 * (2 * Hfin_eva) * Dc / D_eva;
	ReD = inside.density * vairInside_in * Dc / inside.viscosity;
	if (rownum_eva == 1)
	{
		P1 = 1.9 - 0.23 * log(ReD);
		P2 = -0.236 + 0.126 * log(ReD);
		part1 = pow(tubeHengPitch_eva / tubeZongPitch_eva, P1);
		part2 = pow(finPitch_eva / Dc, -1.084);
		part3 = pow(finPitch_eva / Dh, -0.786);
		part4 = pow(finPitch_eva / tubeHengPitch_eva, P2);
		j_heatx = 0.108 * pow(ReD, -0.29) * part1 * part2 * part3 * part4;
	}
	else
	{
		P3 = -0.361 - 0.042 * rownum_eva / log(ReD) + 0.158 * log(rownum_eva * pow(finPitch_eva / Dc, 0.41));
		P4 = -1.224 - 0.076 * pow(tubeZongPitch_eva / Dh, 1.42) / log(ReD);
		P5 = -0.083 + 0.058 * rownum_eva / log(ReD);
		P6 = -5.735 + 1.21 * log(ReD / rownum_eva);
		part1 = 0.086 * pow(ReD, P3) * pow(rownum_eva, P4);
		part2 = pow(finPitch_eva / Dc, P5);
		part3 = pow(finPitch_eva / Dh, P6);
		part4 = pow(finPitch_eva / tubeHengPitch_eva, -0.93);
		j_heatx = part1 * part2 * part3 * part4;
	}
	Pr_fin = 1000 * inside.HeatCapacity_p * inside.viscosity / inside.thermalConduct;
	a1 = j_heatx * inside.density * vairInside_in * (1000 * inside.HeatCapacity_p) / pow(Pr_fin, 2.0 / 3);

	/*肋效率计算*/
	mH = pow((2 * a1 / thermalConduct_fin / finThickness_eva), 0.5) * Hfin_eva;
	th_mH = (pow(E_value, mH) - pow(E_value, -1 * mH)) / (pow(E_value, mH) + pow(E_value, -1 * mH));
	nf = th_mH / mH;
	Ar = PI * D_eva * L_eva * tubenum_eva;
	Af = 2 * (2 * Hfin_eva) * (tubeZongPitch_eva * (tubenum_eva + 1)) * (int(L_eva / finPitch_eva) + 1) - 2 * 0.25 * PI * D_eva * D_eva * tubenum_eva * rownum_eva;
	n0 = (Ar + nf * Af) / (Ar + Af);
	n0Arf = (Ar + Af) * n0;

	/*空气侧阻力损失计算*/
	F1 = -0.764 + 0.739 * tubeHengPitch_eva / tubeZongPitch_eva + 0.177 * finPitch_eva / Dc - 0.00758 / rownum_eva;
	F2 = -15.689 + 64.021 / log(ReD);
	F3 = 1.696 - 15.695 / log(ReD);
	f_heatx = 0.0267 * pow(ReD, F1) * pow(tubeHengPitch_eva / tubeZongPitch_eva, F2) * pow(finPitch_eva / Dc, F3);
	dpair = f_heatx * Af * inside.density * vairInside_in * vairInside_in / (2 * L_eva * tubenum_eva * (tubeZongPitch_eva - 2 * D_eva));

	for (ii = 0; ii < unitnum; ii++)
	{
		/*制冷剂侧传热系数与压降计算*/
		if (Unit.quality < 1) // 干度判断， x < 1 时计算管内沸腾传热系数
		{
			lastUnit = Liquid = Gas = Unit;
			Liquid.pSatProperty(1);
			Gas.pSatProperty(2);
			A_tube = 0.25 * PI * d_eva * d_eva;
			v_liq = Liquid.qm * (1 - Liquid.quality) / Liquid.density / A_tube;
			Re_liq = Liquid.density * v_liq * d_eva / Liquid.viscosity;
			Pr_liq = (1000 * Liquid.HeatCapacity_p) * Liquid.viscosity / Liquid.thermalConduct;
			Nul_liq = 0.023 * pow(Re_liq, 0.8) * pow(Pr_liq, 0.4);
			al = Nul_liq * Liquid.thermalConduct / d_eva;

			xtt = pow(((1 - Unit.quality) / Unit.quality), 0.9) * pow((Gas.density / Liquid.density), 0.5) * pow((Liquid.viscosity / Gas.viscosity), 0.1);
			pr = Unit.pressure / Unit.p_cri;

			/*迭代计算管内沸腾传热系数 a0 */
			a0_esti = atp = a0_last;
			do
			{
				a0_esti = atp;
				k_total = 1.0 / ((1 / a0_esti) + (0.5 * (D_eva - d_eva) / thermalConduct_Wall) + A0_total / (n0Arf * a1)); // 内壁为基准的 k
				Unit.heatFlux = heatx(Unit, inside, k_total, A0_unit, unitLength, dpair, &conWaterUnit, &th_out);

				hflux_d = Unit.heatFlux / (A0_unit);
				Bo = hflux_d / ((Gas.enthalpy - Liquid.enthalpy) * 1000 * (Unit.qm / A_tube));
				E = 1 + 24000 * pow(Bo, 1.16) + 1.37 * pow(xtt, -0.86);
				S = 1 / (1 + 1.15 * 1e-6 * E * E * pow(Re_liq, 1.17));
				anb = 55 * pow(pr, 0.12) * pow((-log10(pr)), -0.55) * pow(Unit.molmass, -0.5) * pow(hflux_d, 0.67);
				atp = E * al + S * anb;
			} while (abs(a0_esti - atp) > 0.1);
			a0 = a0_last = a0_esti;

			/*压降计算*/
			factorL = factorG = 0;
			if (Re_liq > 4000) factorL = 1 + 20 / xtt + 1 / (xtt * xtt);
			else factorG = 1 + 20 * xtt + xtt * xtt;
			G_gas = Gas.qm * Gas.quality / A_tube;
			Re_gas = G_gas * d_eva / Gas.viscosity;
			G_liq = Liquid.qm * (1 - Liquid.quality) / A_tube;
			if (Re_liq < 2000) fl = 64 / Re_liq;
			else fl = 0.3164 * pow(Re_liq, -0.25);
			if (Re_gas < 2000) fg = 64 / Re_gas;
			else fg = 0.3164 * pow(Re_gas, -0.25);
			dp_liq = fl * unitLength * G_liq * G_liq / (2 * d_eva * Liquid.density);
			dp_gas = fg * unitLength * G_gas * G_gas / (2 * d_eva * Gas.density);
			dpfri = factorL * dp_liq + factorG * dp_gas; //其中一个因子为0
			dpunit = dpfri + pow((G_gas + G_liq), 2) * (1 / Unit.density - 1 / lastUnit.density);
		}
		else // 计算管内湍流强制对流传热系数（干度 x = 1 时）
		{
			Unit.tpProperty();
			v_gas = Unit.qm / Unit.density / (0.25 * PI * d_eva * d_eva);
			Re_gas = Unit.density * v_gas * d_eva / Unit.viscosity;
			Pr_gas = (1000 * Unit.HeatCapacity_p) * Unit.viscosity / Unit.thermalConduct;
			Nu_gas = 0.023 * pow(Re_gas, 0.8) * pow(Pr_gas, 0.4);
			a0 = Nu_gas * Unit.thermalConduct / d_eva;
			/*压降计算*/
			f = 0.3164 * pow(Re_gas, -0.25);
			dpunit = f * unitLength * v_gas * v_gas * Unit.density / (2 * d_eva);
		}
		/*单元传热量计算*/
		k_total = 1.0 / ((1 / a0) + (0.5 * (D_eva - d_eva) / thermalConduct_Wall) + A0_total / (n0Arf * a1));
		Qunit = heatx(Unit, inside, k_total, A0_unit, unitLength, dpair, &conWaterUnit, &th_out);
		Qtube += Qunit;
		condensationWater += conWaterUnit;
		tSum_air += th_out;
		if (Unit.quality >= 1) Unit.quality = 1;
		Vliq_eva += unitLength * A_tube * (1 - Unit.quality) / Liquid.density * (Unit.quality / Gas.density + (1 - Unit.quality) / Liquid.density);

		lastUnit = Unit;
		newUnit = Unit;
		newUnit.pressure = Unit.pressure - dpunit;
		newUnit.enthalpy = Unit.enthalpy + Qunit / 1000 / Unit.qm;
		newUnit.phProperty();
		Unit = newUnit;
	}
	/*计算换热量与出口参数*/
	point[num + 1] = newUnit;
	Q_eva = Qtube * tubenum_eva;
	tAve_air = tSum_air / unitnum;
}

/*蒸发器换热器模块*/
double heatx(struct refri unit, struct wetair inside, double k, double A, double unitLength, double dpair, double* conWaterUnit, double* th_out_return)
{
	double th_in, th_out, tc_in, tc_out, dtm, th_out_reslut;
	double t_esti, min_esti, max_esti, Q_esti, Q_comp, Q_result;  // t_esti 是 th_out
	int n = 0, judge = 0;

	inside.tpProperty();
	inside.qm = inside.density * vairInside_in * tubeZongPitch_eva * unitLength;

	th_in = inside.tem;
	tc_in = unit.tem;
	min_esti = tc_in; max_esti = th_in;

	if ((heatx_lastmin != 0) && (heatx_lastmax != 0)) // 借用上次计算结果，加快收敛速度
	{
		judge = 0;
		heatx_lastmin = heatx_lastmin - (heatx_lastmax - heatx_lastmin);
		heatx_lastmax = heatx_lastmax + (heatx_lastmax - heatx_lastmin);
		if (heatx_lastmin < tc_in) heatx_lastmin = tc_in;
		if (heatx_lastmax > th_in) heatx_lastmax = th_in;
		min_esti = heatx_lastmin; max_esti = heatx_lastmax;
		for (int i = 0; i < 2; i++)
		{
			if (i == 0) t_esti = heatx_lastmin;
			else t_esti = heatx_lastmax;
			inside.tem_2 = t_esti;
			inside.pressure_2 = inside.pressure - dpair;
			inside.condensation();

			*conWaterUnit = inside.delta_dd * (inside.qm / (inside.dd + 1));
			if (inside.delta_dd > 0) Q_esti = 1000 * (*conWaterUnit * inside.deltaenthalpy_h2o_2 + (inside.qm - *conWaterUnit) * (inside.enthalpy - inside.enthalpy_2));
			else Q_esti = 1000 * inside.qm * (inside.enthalpy - inside.enthalpy_2);

			if (unit.quality < 1)
			{
				tc_out = tc_in;
			}
			else
			{
				tc_out = tc_in + Q_esti / 1000 / unit.qm / unit.HeatCapacity_p;
				if (tc_out >= th_in)
				{
					judge += 2; break;
				}
			}
			dtm = ((th_in - tc_out) - (t_esti - tc_in)) / log((th_in - tc_out) / (t_esti - tc_in));
			Q_comp = k * A * dtm;

			if (Q_esti >= Q_comp)
			{
				judge++;
			}
			else judge--;
		}
		if (judge == 0)
		{
			min_esti = heatx_lastmin; max_esti = heatx_lastmax;
		}
		else
		{
			min_esti = tc_in; max_esti = th_in;
		}
	}

	while (n < 15)
	{
		t_esti = 0.5 * (min_esti + max_esti);
		inside.tem_2 = t_esti;
		inside.pressure_2 = inside.pressure - dpair;
		inside.condensation();

		*conWaterUnit = inside.delta_dd * (inside.qm / (inside.dd + 1));
		if (inside.delta_dd > 0) Q_esti = 1000 * (*conWaterUnit * inside.deltaenthalpy_h2o_2 + (inside.qm - *conWaterUnit) * (inside.enthalpy - inside.enthalpy_2));
		else Q_esti = 1000 * inside.qm * (inside.enthalpy - inside.enthalpy_2);

		//Q_esti = 1000 * inside.qm * inside.HeatCapacity_p * (th_in - t_esti);

		if (unit.quality < 1)
		{
			tc_out = tc_in;
		}
		else
		{
			tc_out = tc_in + Q_esti / 1000 / unit.qm / unit.HeatCapacity_p;
			if (tc_out >= th_in)
			{
				min_esti = t_esti;
				max_esti = max_esti;
				continue;
			}
		}
		dtm = ((th_in - tc_out) - (t_esti - tc_in)) / log((th_in - tc_out) / (t_esti - tc_in));
		Q_comp = k * A * dtm;

		if (Q_esti >= Q_comp)
		{
			min_esti = t_esti;
			max_esti = max_esti;
		}
		else
		{
			min_esti = min_esti;
			max_esti = t_esti;
		}
		n++;
		if (abs(Q_esti - Q_comp) < 0.01) break;
	}

	//
	if (Q_esti >= Q_comp)
	{
		heatx_lastmin = min_esti - (max_esti - min_esti);
		heatx_lastmax = max_esti;
	}
	else
	{
		heatx_lastmin = min_esti;
		heatx_lastmax = max_esti + (max_esti - min_esti);
	}
	//

	Q_result = Q_esti;
	*th_out_return = th_out_reslut = t_esti;
	return Q_result;
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

