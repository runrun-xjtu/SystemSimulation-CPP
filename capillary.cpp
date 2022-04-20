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
extern double L_cap, d_cap, epsilon; // ëϸ�ܳ����ھ����ֲڶ�

/*�ⲿ����ģ�麯������*/

extern void compressor(struct refri point[], int num);
extern void condenser(struct refri point[], struct wetair outside, int num);
extern void capillary(struct refri point[], int num);
extern void evaporator(struct refri point[], struct wetair inside, int num);

/*ëϸ��ģ��*/
void capillary(struct refri point[], int num)
{
	struct refri unit;            //ëϸ�ܼ���ṹ��
	struct refri unitLast;        //
	struct refri unitLastLast;        //
	struct refri origin;        //
	struct refri soundVelocity, S_liq, S_gas;
	struct refri Liquid, Gas;
	double qmEsti, lengthUnit, v_Sound;      //������������Ԫ���ȣ�����
	double velocity, G, Re, f, A, B, dpVelcocity, dpFriction; //�ٶȣ������ܶ�
	int n = 100, i, j;                          //΢Ԫ��
	int errorID = 0;

	qmEsti = 0;
	if (point[num].quality <= 0) point[num].tpProperty();       //��ʼ������Լ���
	else point[num].phProperty();
	origin = point[num];
	lengthUnit = L_cap / n;

	for (i = 0; i < 6; i++)    //6λ���ȼ�������
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
				/*�ٶȡ������ܶȡ���ŵ������*/
				velocity = qmEsti / unitLast.density / (0.25 * PI * d_cap * d_cap);
				G = velocity * unitLast.density;
				Re = unitLast.density * velocity * d_cap / unitLast.viscosity;

				/*ѹ������*/
				A = pow(2.457 * log(1 / (pow(7 / Re, 0.9) + 0.27 * epsilon / d_cap)), 16);
				B = pow(37560 / Re, 16);
				f = 8 * pow(pow(8 / Re, 12) + 1.0 / pow(A + B, 1.5), 1.0 / 12);
				dpVelcocity = G * G * (1 / unitLast.density - 1 / unitLastLast.density);
				dpFriction = f * velocity * velocity * lengthUnit * unitLast.density / (2 * d_cap);

				/*����unit���ݣ�ph������*/
				unit.enthalpy = unitLast.enthalpy;
				unit.pressure -= dpVelcocity + dpFriction;
				if (unit.pressure <= point[num + 1].pressure) unit.pressure = point[num + 1].pressure - 10;
				unit.phProperty();
				unitLastLast = unitLast;
				unitLast = unit;
			}
			/*�ɶȴ���0ʱ���ټ��㣬�ж��Ƿ�ﵽ����*/
			if (unit.quality > 0)
			{
				soundVelocity.entropy = unit.entropy;
				soundVelocity.pressure = unit.pressure + 50; //ȡѹ��dp=50,��d(density)
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
			/*����ѹ�����㣬�ж��Ƿ�ﵽ��ѹ*/
			if (unit.pressure <= point[num + 1].pressure)
			{
				qmEsti -= pow(0.1, 2 + i);
				break;
			}
		}
		if (i == 0 && j == 10)
		{
			cout << "��������������Χ��qm > 0.1 kg/s��" << endl;
			errorID = 2; break;
		}
	}
	if (errorID == 0)
	{
		point[num + 1].qm = point[num].qm = qmEsti;  //�������
	}
}