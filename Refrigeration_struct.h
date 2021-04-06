#pragma once
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
	/* 程序内计算单位 t-K, p-Pa, d-kg/m3, h-kJ/kg, s-kJ/(kg.K) vis- Pa.s */
	/* REFPROP函数单位 t-K, p-kPa, d-mol/L, h-J/mol, s-J/(mol.K) vis- uPa.s */
	double tem, pressure, p_cri, molmass;                  //温度 压力 临界压力 摩尔质量
	double quality, density, densityLiq, densityGas;       //干度 密度
	double viscosity, viscosityLiq, viscosityGas;		   //动力粘度
	double enthalpy, entropy, k_com;					   //焓 熵 绝热指数
	double HeatCapacity_p, HeatCapacity_v,thermalConduct;  //定压/定容比热容 热导率
	double qm, heatFlux;                                   //流量 热流量

	/*临时参数（REFPROP函数使用）*/
	typedef void(__stdcall* fp_SETUPdllTYPE)(long&, char*, char*, char*, long&, char*, long, long, long, long); fp_SETUPdllTYPE SETUPdll;
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
		HINSTANCE RefpropdllInstance = LoadLibrary(L"refprop.dll"); //加载DLL

		/*可能用到的REFPROP函数*/
		SETUPdll = (fp_SETUPdllTYPE)GetProcAddress(RefpropdllInstance, "SETUPdll");
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
		molmass = wm; p_cri = pc;
	}

	/*根据 p 饱和态求解（i=1 液 i=2 气）*/
	void pSatProperty(long i)
	{
		p = 0.001 * pressure;
		SATPdll(p, x, i, t, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.01; else if (i == 2) t += 0.001;
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
		t = tem;
		SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.01; else if (i == 2) t += 0.001;
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
		t = tem; p = 0.001 * pressure;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; density = d * molmass;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta ;
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass;
		k_com = cp / cv; thermalConduct = tcx;
	}

	/*根据 p h 求解状态参数*/
	void phProperty()
	{
		p = 0.001 * pressure; h = molmass * enthalpy;
		PHFLSHdll(p, h, x, t, d, dl, dv, xliq, xvap, q, e, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta ;
		density = d * molmass; densityLiq = dl * molmass; densityGas = dv * molmass;
		if (quality > 0) density = densityGas * quality + densityLiq * (1 - quality);
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass; k_com = cp / cv;
	}

	/*根据 p s 求解状态参数*/
	void psProperty()
	{
		p = 0.001 * pressure; s = molmass * entropy;
		PSFLSHdll(p, s, x, t, d, dl, dv, xliq, xvap, q, e, h, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta ;
		density = d * molmass; densityLiq = dl * molmass; densityGas = dv * molmass;
		if (quality > 0) density = densityGas * quality + densityLiq * (1 - quality);
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass; k_com = cp / cv;
	}
}refri;


/*调用REFPROP需要参数*/
/*
const long refpropcharlength = 255;const long filepathlength = 255;const long lengthofreference = 3;
const long errormessagelength = 255;const long ncmax = 20;const long numparams = 72;const long maxcoefs = 50;
*/

/*暂不启用REFPROP函数*/
/*
typedef void (__stdcall *fp_ABFL1dllTYPE)(double &,double &,double *,long &,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_ABFL2dllTYPE)(double &,double &,double *,long &,long &,double &,double &,double &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double *,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_ACTVYdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_AGdllTYPE)(double &,double &,double *,double &,double &);
typedef void (__stdcall *fp_CCRITdllTYPE)(double &,double &,double &,double *,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_CP0dllTYPE)(double &,double *,double &);
typedef void (__stdcall *fp_CRITPdllTYPE)(double *,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_CSATKdllTYPE)(long &,double &,long &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_CV2PKdllTYPE)(long &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_CVCPKdllTYPE)(long &,double &,double &,double &,double &);
typedef void (__stdcall *fp_CVCPdllTYPE)(double &,double &,double *,double &,double &);
typedef void (__stdcall *fp_DBDTdllTYPE)(double &,double *,double &);
typedef void (__stdcall *fp_DBFL1dllTYPE)(double &,double &,double *,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_DBFL2dllTYPE)(double &,double &,double *,long &,double &,double &,double &,double &,double &,double *,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_DDDPdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_DDDTdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_DEFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_DHD1dllTYPE)(double &,double &,double *,double &,double &,double &,double &,double &,double &);
typedef void (__stdcall *fp_DHFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_DIELECdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_DOTFILLdllTYPE)(long &,double *,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_DPDD2dllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_DPDDKdllTYPE)(long &,double &,double &,double &);
typedef void (__stdcall *fp_DPDDdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_DPDTKdllTYPE)(long &,double &,double &,double &);
typedef void (__stdcall *fp_DPDTdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_DPTSATKdllTYPE)(long &,double &,long &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_DSFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_ENTHALdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_ENTROdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_ESFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_FGCTYdllTYPE)(double &,double &,double *,double *);
typedef void (__stdcall *fp_FPVdllTYPE)(double &,double &,double &,double *,double &);
typedef void (__stdcall *fp_GERG04dllTYPE)(long &,long &,long &,char*,long );
typedef void (__stdcall *fp_GETFIJdllTYPE)(char*,double *,char*,char*,long ,long ,long );
typedef void (__stdcall *fp_GETKTVdllTYPE)(long &,long &,char*,double *,char*,char*,char*,char*,long ,long ,long ,long ,long );
typedef void (__stdcall *fp_GIBBSdllTYPE)(double &,double &,double *,double &,double &);
typedef void (__stdcall *fp_HSFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,long &,char*,long );

typedef void (__stdcall *fp_LIMITKdllTYPE)(char*,long &,double &,double &,double &,double &,double &,double &,double &,long &,char*,long ,long );
typedef void (__stdcall *fp_LIMITSdllTYPE)(char*,double *,double &,double &,double &,double &,long );
typedef void (__stdcall *fp_LIMITXdllTYPE)(char*,double &,double &,double &,double *,double &,double &,double &,double &,long &,char*,long ,long );
typedef void (__stdcall *fp_MELTPdllTYPE)(double &,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_MELTTdllTYPE)(double &,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_MLTH2OdllTYPE)(double &,double &,double &);
typedef void (__stdcall *fp_NAMEdllTYPE)(long &,char*,char*,char*,long ,long ,long );
typedef void (__stdcall *fp_PDFL1dllTYPE)(double &,double &,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_PDFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_PEFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_PHFL1dllTYPE)(double &,double &,double *,long &,double &,double &,long &,char*,long );

typedef void (__stdcall *fp_PQFLSHdllTYPE)(double &,double &,double *,long &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_PREOSdllTYPE)(long &);
typedef void (__stdcall *fp_PRESSdllTYPE)(double &,double &,double *,double &);
typedef void (__stdcall *fp_PSFL1dllTYPE)(double &,double &,double *,long &,double &,double &,long &,char*,long );

typedef void (__stdcall *fp_PUREFLDdllTYPE)(long &);
typedef void (__stdcall *fp_QMASSdllTYPE)(double &,double *,double *,double &,double *,double *,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_QMOLEdllTYPE)(double &,double *,double *,double &,double *,double *,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_SATDdllTYPE)(double &,double *,long &,long &,double &,double &,double &,double &,double *,double *,long &,char*,long );
typedef void (__stdcall *fp_SATEdllTYPE)(double &,double *,long &,long &,long &,double &,double &,double &,long &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_SATHdllTYPE)(double &,double *,long &,long &,long &,double &,double &,double &,long &,double &,double &,double &,long &,char*,long );

typedef void (__stdcall *fp_SATSdllTYPE)(double &,double *,long &,long &,long &,double &,double &,double &,long &,double &,double &,double &,long &,double &,double &,double &,long &,char*,long );

typedef void (__stdcall *fp_SETAGAdllTYPE)(long &,char*,long );
typedef void (__stdcall *fp_SETKTVdllTYPE)(long &,long &,char*,double *,char*,long &,char*,long ,long ,long );
typedef void (__stdcall *fp_SETMIXdllTYPE)(char*,char*,char*,long &,char*,double *,long &,char*,long ,long ,long ,long ,long );
typedef void (__stdcall *fp_SETMODdllTYPE)(long &,char*,char*,char*,long &,char*,long ,long ,long ,long );
typedef void (__stdcall *fp_SETREFdllTYPE)(char*,long &,double *,double &,double &,double &,double &,long &,char*,long ,long );

typedef void (__stdcall *fp_SPECGRdllTYPE)(double &,double &,double &,double &);
typedef void (__stdcall *fp_SUBLPdllTYPE)(double &,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_SUBLTdllTYPE)(double &,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_SURFTdllTYPE)(double &,double &,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_SURTENdllTYPE)(double &,double &,double &,double *,double *,double &,long &,char*,long );
typedef void (__stdcall *fp_TDFLSHdllTYPE)(double &,double &,double *,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_TEFLSHdllTYPE)(double &,double &,double *,long &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_THERM0dllTYPE)(double &,double &,double *,double &,double &,double &,double &,double &,double &,double &,double &,double &);
typedef void (__stdcall *fp_THERM2dllTYPE)(double &,double &,double *,double &,double &,double &,double &,double &,double &,double &,double *,double &,double &,double &,double &,double &,double &,double &,double &,double &,double &,double &,double &,double &,double &);
typedef void (__stdcall *fp_THERM3dllTYPE)(double &,double &,double *,double &,double &,double &,double &,double &,double &,double &,double &,double &,double &);
typedef void (__stdcall *fp_THERMdllTYPE)(double &,double &,double *,double &,double &,double &,double &,double &,double &,double &,double &);
typedef void (__stdcall *fp_THFLSHdllTYPE)(double &,double &,double *,long &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );

typedef void (__stdcall *fp_TPRHOdllTYPE)(double &,double &,double *,long &,long &,double &,long &,char*,long );
typedef void (__stdcall *fp_TQFLSHdllTYPE)(double &,double &,double *,long &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );

typedef void (__stdcall *fp_TSFLSHdllTYPE)(double &,double &,double *,long &,double &,double &,double &,double &,double *,double *,double &,double &,double &,double &,double &,double &,long &,char*,long );
typedef void (__stdcall *fp_VIRBdllTYPE)(double &,double *,double &);
typedef void (__stdcall *fp_VIRCdllTYPE)(double &,double *,double &);
typedef void (__stdcall *fp_WMOLdllTYPE)(double *,double &);
typedef void (__stdcall *fp_XMASSdllTYPE)(double *,double *,double &);
typedef void (__stdcall *fp_XMOLEdllTYPE)(double *,double *,double &);

//Define explicit function pointers
fp_ABFL1dllTYPE ABFL1dll;
fp_ABFL2dllTYPE ABFL2dll;
fp_ACTVYdllTYPE ACTVYdll;
fp_AGdllTYPE AGdll;
fp_CCRITdllTYPE CCRITdll;
fp_CP0dllTYPE CP0dll;
fp_CRITPdllTYPE CRITPdll;
fp_CSATKdllTYPE CSATKdll;
fp_CV2PKdllTYPE CV2PKdll;
fp_CVCPKdllTYPE CVCPKdll;
fp_CVCPdllTYPE CVCPdll;
fp_DBDTdllTYPE DBDTdll;
fp_DBFL1dllTYPE DBFL1dll;
fp_DBFL2dllTYPE DBFL2dll;
fp_DDDPdllTYPE DDDPdll;
fp_DDDTdllTYPE DDDTdll;
fp_DEFLSHdllTYPE DEFLSHdll;
fp_DHD1dllTYPE DHD1dll;
fp_DHFLSHdllTYPE DHFLSHdll;
fp_DIELECdllTYPE DIELECdll;
fp_DOTFILLdllTYPE DOTFILLdll;
fp_DPDD2dllTYPE DPDD2dll;
fp_DPDDKdllTYPE DPDDKdll;
fp_DPDDdllTYPE DPDDdll;
fp_DPDTKdllTYPE DPDTKdll;
fp_DPDTdllTYPE DPDTdll;
fp_DPTSATKdllTYPE DPTSATKdll;
fp_DSFLSHdllTYPE DSFLSHdll;
fp_ENTHALdllTYPE ENTHALdll;
fp_ENTROdllTYPE ENTROdll;
fp_ESFLSHdllTYPE ESFLSHdll;
fp_FGCTYdllTYPE FGCTYdll;
fp_FPVdllTYPE FPVdll;
fp_GERG04dllTYPE GERG04dll;
fp_GETFIJdllTYPE GETFIJdll;
fp_GETKTVdllTYPE GETKTVdll;
fp_GIBBSdllTYPE GIBBSdll;
fp_HSFLSHdllTYPE HSFLSHdll;

fp_LIMITKdllTYPE LIMITKdll;
fp_LIMITSdllTYPE LIMITSdll;
fp_LIMITXdllTYPE LIMITXdll;
fp_MELTPdllTYPE MELTPdll;
fp_MELTTdllTYPE MELTTdll;
fp_MLTH2OdllTYPE MLTH2Odll;
fp_NAMEdllTYPE NAMEdll;
fp_PDFL1dllTYPE PDFL1dll;
fp_PDFLSHdllTYPE PDFLSHdll;
fp_PEFLSHdllTYPE PEFLSHdll;
fp_PHFL1dllTYPE PHFL1dll;

fp_PQFLSHdllTYPE PQFLSHdll;
fp_PREOSdllTYPE PREOSdll;
fp_PRESSdllTYPE PRESSdll;
fp_PSFL1dllTYPE PSFL1dll;

fp_PUREFLDdllTYPE PUREFLDdll;
fp_QMASSdllTYPE QMASSdll;
fp_QMOLEdllTYPE QMOLEdll;
fp_SATDdllTYPE SATDdll;
fp_SATEdllTYPE SATEdll;
fp_SATHdllTYPE SATHdll;

fp_SATSdllTYPE SATSdll;

fp_SETAGAdllTYPE SETAGAdll;
fp_SETKTVdllTYPE SETKTVdll;
fp_SETMIXdllTYPE SETMIXdll;
fp_SETMODdllTYPE SETMODdll;
fp_SETREFdllTYPE SETREFdll;

fp_SPECGRdllTYPE SPECGRdll;
fp_SUBLPdllTYPE SUBLPdll;
fp_SUBLTdllTYPE SUBLTdll;
fp_SURFTdllTYPE SURFTdll;
fp_SURTENdllTYPE SURTENdll;
fp_TDFLSHdllTYPE TDFLSHdll;
fp_TEFLSHdllTYPE TEFLSHdll;
fp_THERM0dllTYPE THERM0dll;
fp_THERM2dllTYPE THERM2dll;
fp_THERM3dllTYPE THERM3dll;
fp_THERMdllTYPE THERMdll;
fp_THFLSHdllTYPE THFLSHdll;

fp_TPRHOdllTYPE TPRHOdll;
fp_TQFLSHdllTYPE TQFLSHdll;

fp_TSFLSHdllTYPE TSFLSHdll;
fp_VIRBdllTYPE VIRBdll;
fp_VIRCdllTYPE VIRCdll;
fp_WMOLdllTYPE WMOLdll;
fp_XMASSdllTYPE XMASSdll;
fp_XMOLEdllTYPE XMOLEdll;

*/