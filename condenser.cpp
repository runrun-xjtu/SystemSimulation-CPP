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
extern double tairOutside_in, vairOutside_in;

extern double V_con, Vliq_con;

/*外部引入结构参数*/
extern double tubenum_con, rownum_con, L_con, D_con, d_con;
extern double finPitch_con, finThickness_con, tubeZongPitch_con, tubeHengPitch_con, Hfin_con;
extern double thermalConduct_Wall, thermalConduct_fin;

/*外部引入模块函数声明*/

extern void compressor(struct refri point[], int num);
extern void condenser(struct refri point[], struct wetair outside, int num);
extern void capillary(struct refri point[], int num);
extern void evaporator(struct refri point[], struct wetair inside, int num);

double heatxcon(struct refri unit, struct wetair outside, double k, double A, double unitLength);

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
        if (unit.quality < 1 && unit.quality >= 0 )
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

