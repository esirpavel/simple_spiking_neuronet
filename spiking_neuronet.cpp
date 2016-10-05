/*
 * neuronet.cpp
 *
 *  Created on: 04.11.2013
 *      Author: Pavel Esir
 */
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

const float h     = .5f; 	  // временной шаг интегрирования
const int   Tsim  = 1000/.5f; // время симуляции в дискретных отсчетах
const int   Nexc  = 100; 	  // Количество возбуждающих (excitatory) нейронов
const int   Ninh  = 25;  	  // Количество тормозных (inhibitory) нейронов
const int   Nneur = Nexc + Ninh;
const int   Ncon  = Nneur*Nneur*0.1f; // Количество сязей,
			   // 0.1 это вероятность связи между 2-мя случайными нейронами

float Vms[Nneur][Tsim]; // мембранные потенциалы
float Ums[Nneur][Tsim]; // вспомогательные переменные модели Ижикевича
float Iex[Nneur]; 		// внешний постоянный ток приложенный к нейрону
float Isyn[Nneur]; 		// синаптический ток на каждый нейрон

int pre_conns[Ncon];    // индексы пресинаптических нейронов
int post_conns[Ncon];   // индексы постсинаптических нейронов
float weights[Ncon]; 	// веса связей
float y[Ncon][Tsim];    // переменная модулирующая синаптический ток в зависимости от спайков на пресинапсе

float psc_excxpire_time = 4.0f; // характерное вермя спадания постсинаптического тока, мс
float minWeight = 50.0f; // веса, размерность пкА
float maxWeight = 100.0f;

// Параметры нейрона
float Iex_max = 40.0f; // максимальный приложенный к нейрону ток 50 пкА
float a		  = 0.02f;
float b		  = 0.5f;
float c		  = -40.0f; // значение мембранного потенциала до которого он сбрасываеться после спайка
float d 	  = 100.0f;
float k		  = 0.5f;
float Vr	  = -60.0f;
float Vt	  = -45.0f;
float Vpeak	  = 35.0f;  // максимальное значение мембранного потенциала, при котором происходит сброс до значения с
float V0	  = -60.0f; // начальное значение для мембранного потенциала
float U0	  = 0.0f;   // начальное значение для вспомогательной переменной
float Cm      = 50.0f;  // электрическая ёмкость нейрона, размерность пкФ

float spike_times[Nneur*Tsim]; // времена спайков
int spike_neurons[Nneur*Tsim]; // соответвующие номера нейронов
int spike_num = 0;

using namespace std;

void init_connections(){
	for (int con_idx = 0; con_idx < Ncon; ){
		// случайно выбираем постсипантические и пресинаптические нейроны
		int pre = rand() % Nneur;
		int post = rand() % Nneur;
		pre_conns[con_idx] = pre;
		post_conns[con_idx] = post;
		weights[con_idx] = (rand() % ((int)(maxWeight - minWeight)*10))/10.0f + minWeight;
		if (pre >= Nexc){
			// если пресинаптический нейрон тормозный то вес связи идет со знаком минус
			weights[con_idx] = -weights[con_idx];
		}
		y[con_idx][0] = 0.0f;
		con_idx++;
	}
}

void init_neurons(){
	for (int neur_idx = 0; neur_idx < Nneur; neur_idx++){
		// случайно разбрасываем приложенные токи
		Iex[neur_idx] = ((float) rand() / RAND_MAX) * Iex_max;
		Isyn[neur_idx] = 0.0f;
		Vms[neur_idx][0] = V0;
		Ums[neur_idx][0] = U0;
	}
}

float izhik_Vm(int neuron, int time){
	return (k*(Vms[neuron][time] - Vr)*(Vms[neuron][time] - Vt) - Ums[neuron][time] + Iex[neuron] + Isyn[neuron])/Cm;
}

float izhik_Um(int neuron, int time){
	return a*(b*(Vms[neuron][time] - Vr) - Ums[neuron][time]);
}

void save2file(){
	ofstream res_file;
	res_file.open("rastr.csv");
	for (int k = 0; k < spike_num; k++){
		res_file << spike_times[k] << "; " << spike_neurons[k] + 1 << "; " << endl;
	}
	res_file.close();

	// Вычисление среднего по всей сети мембранного потенциала в каждый момент времени
	// нечто наподобие электроэнцефалографии
	res_file.open("oscill.csv");
	for (int t = 0; t < Tsim; t++){
		float Vm_mean= 0.0f;
		for (int m = 0; m < Nneur; m++){
			Vm_mean += Vms[m][t];
		}
		Vm_mean /= Nneur;
		res_file << t*h << "; " << Vm_mean << "; " << endl;
	}
	res_file.close();
}

int main(){
	init_connections();
	init_neurons();
	float expire_coeff = exp(-h/psc_excxpire_time);
	for (int t = 1; t < Tsim; t++){
		// проходим по всем нейронам
		for (int neur = 0; neur < Nneur; neur++){
			Vms[neur][t] = Vms[neur][t-1] + h*izhik_Vm(neur, t-1);
			Ums[neur][t] = Ums[neur][t-1] + h*izhik_Um(neur, t-1);
			Isyn[neur] = 0.0f;

			if (Vms[neur][t-1] >Vpeak){
				Vms[neur][t] = c;
				Ums[neur][t] = Ums[neur][t-1] + d;
				spike_times[spike_num] = t*h;
				spike_neurons[spike_num] = neur;
				spike_num++;
			}
		}

		// проходим по всем связям
		for (int con = 0; con < Ncon; con++){
			y[con][t] = y[con][t-1]*expire_coeff;

			if (Vms[pre_conns[con]][t-1] > Vpeak){
				y[con][t] += 1.0f;
			}
			Isyn[post_conns[con]] += y[con][t]*weights[con];
		}
	}
	save2file();
	return 0;
}
