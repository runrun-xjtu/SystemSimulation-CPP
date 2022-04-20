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
extern double P_com, Q_eva;
extern double tairInside_in, vairInside_in, condensationWater, V_eva, Vliq_eva;

/*外部引入结构参数*/
extern double tubenum_eva, rownum_eva, L_eva, D_eva, d_eva;
extern double finPitch_eva, finThickness_eva, tubeZongPitch_eva, tubeHengPitch_eva, Hfin_eva;
extern double thermalConduct_Wall, thermalConduct_fin;
extern double mSum_air, tSum_air, tAve_air, mFrac_H2O_out;
/*外部引入模块函数声明*/

extern void compressor(struct refri point[], int num);
extern void condenser(struct refri point[], struct wetair outside, int num);
extern void capillary(struct refri point[], int num);
extern void evaporator(struct refri point[], struct wetair inside, int num);

double heatx(struct refri unit, struct wetair inside, double k, double A, double unitLength, double dpair, double *conWaterUnit, double * th_out_return);

double heatx_lastmin = 0, heatx_lastmax = 0;

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
    double a0_esti, a0_last, atp, anb, E, S, xtt, pr, Bo,hflux_d ,Re; //管内沸腾传热系数计算
    double v_gas, Re_gas, Pr_gas, Nu_gas, G_gas, th_out;                         //管内沸腾传热系数计算
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
            Re_gas = G_gas * d_eva  / Gas.viscosity;
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

double heatx(struct refri unit, struct wetair inside, double k, double A, double unitLength, double dpair, double * conWaterUnit, double* th_out_return)
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

/* a0 求解算法 使用 *迭代法 */

/* 其他求解算法*/
/*二分搜索法*/
/*
            a0_esti = 2000;
            do
            {
                a0_esti = 2 * a0_esti;
                k_total = 1.0 / ((1 / a0_esti) + (0.5 * (D_eva - d_eva) / thermalConduct_Wall) + A0_total / (n0Arf * a1)); // 内壁为基准的 k
                Unit.heatFlux = heatx(Unit, inside, k_total, A0_unit, unitLength, dpair, &conWaterUnit);

                hflux_d = Unit.heatFlux / (A0_unit);
                Bo = hflux_d / ((Gas.enthalpy - Liquid.enthalpy) * 1000 * (Unit.qm / A_tube));
                E = 1 + 24000 * pow(Bo, 1.16) + 1.37 * pow(xtt, -0.86);
                S = 1 / (1 + 1.15 * 1e-6 * E * E * pow(Re_liq, 1.17));
                anb = 55 * pow(pr, 0.12) * pow((-log10(pr)), -0.55) * pow(Unit.molmass, -0.5) * pow(hflux_d, 0.67);
                atp = E * al + S * anb;
            } while (a0_esti < atp);
            esti_min = 0; esti_max = a0_esti;
            j = 0;
            do
            {
                j++;
                a0_esti = 0.5 * (esti_min + esti_max);
                k_total = 1.0 / ((1 / a0_esti) + (0.5 * (D_eva - d_eva) / thermalConduct_Wall) + A0_total / (n0Arf * a1)); // 内壁为基准的 k
                Unit.heatFlux = heatx(Unit, inside, k_total, A0_unit, unitLength, dpair, &conWaterUnit);

                hflux_d = Unit.heatFlux / (A0_unit);
                Bo = hflux_d / ((Gas.enthalpy - Liquid.enthalpy) * 1000 * (Unit.qm / A_tube));
                E = 1 + 24000 * pow(Bo, 1.16) + 1.37 * pow(xtt, -0.86);
                S = 1 / (1 + 1.15 * 1e-6 * E * E * pow(Re_liq, 1.17));
                anb = 55 * pow(pr, 0.12) * pow((-log10(pr)), -0.55) * pow(Unit.molmass, -0.5) * pow(hflux_d, 0.67);
                atp = E * al + S * anb;
                if (a0_esti < atp)
                {
                    esti_min = a0_esti;
                    esti_max = esti_max;
                }
                else if (a0_esti > atp)
                {
                    esti_min = esti_min;
                    esti_max = a0_esti;
                }
                else break;
            } while (abs(a0_esti - atp) > 0.1);
            cout << j << endl;
            a0 = a0_esti;
            */

/*10分搜索 + 迭代搜索法*/
            /*
            for (i = 0; i < 5; i++)  //精度 1e-1 = 0.1
            {
                for (j = 0; j < 300; j++)
                {
                    a0_compare += pow(10, 4 - i);

                    k_total = 1.0 / ((1 / a0_esti) + (0.5 * (D_eva-d_eva) / thermalConduct_Wall) + A0_total / (n0Arf * a1)); // 内壁为基准的 k

                    Unit.heatFlux = heatx(Unit, inside, k_total, A0_unit, unitLength, dpair, &conWaterUnit);
                    hflux_d = Unit.heatFlux / (A0_unit);
                    Bo = hflux_d / ((Gas.enthalpy - Liquid.enthalpy) * 1000 * (Unit.qm / A_tube));
                    E = 1 + 24000 * pow(Bo, 1.16) + 1.37 * pow(xtt, -0.86);
                    S = 1 / (1 + 1.15 * 1e-6 * E * E * pow(Re_liq, 1.17));
                    anb = 55 * pow(pr, 0.12) * pow((-log10(pr)), -0.55) * pow(Unit.molmass, -0.5) * pow(hflux_d, 0.67);
                    atp = E * al + S * anb;

                    if (atp >= a0_esti && atp <= a0_compare)
                    {
                        a0_compare -= pow(10, 4 - i);
                        a0 = atp; break;
                    }
                    a0_esti += pow(10, 4 - i);
                    a0_compare += (atp - a0_esti) * 0.8;
                    a0_esti += (atp - a0_esti) * 0.8;
                }
                //cout << j << endl;
                if (i == 0 && j == 100)
                {
                    cout << "管内沸腾传热系数超出计算域" << endl;
                    break;
                }
            }
            */