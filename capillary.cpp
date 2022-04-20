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
extern double L_cap, d_cap, epsilon; // 毛细管长，内径，粗糙度

/*外部引入模块函数声明*/

extern void compressor(struct refri point[], int num);
extern void condenser(struct refri point[], struct wetair outside, int num);
extern void capillary(struct refri point[], int num);
extern void evaporator(struct refri point[], struct wetair inside, int num);

/*毛细管模块*/
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