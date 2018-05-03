#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<assert.h>
#define PI 3.1415926
#define sizepop 50
#define maxgen 500
#define pcross 0.6
#define pmutation 0.1
#define lenchrom 16

int chrom[sizepop][lenchrom];
int fitness[sizepop];
double fitness_prob[sizepop];
int bestfitness[maxgen];
int gbest_pos[lenchrom];
double average_best[maxgen + 1];
int gbest;
int gbest_index;
long int bound_up[lenchrom] = {100, 100, 500, 2000, 2000, 1000, 500,  2000, 200,  50,   50,  50,  200,  50,  50,  50}; //FIXME: improve readable.
long int bound_down[lenchrom] = {1,  1,   1,   10,   10,   1,   100,   501,  10, -50,  -50,  -50,  10, -50, -50, -50};

int fit_func(int *arr)
{
	char run_buf[100];
	char return_buf[100];
	char *shell_process = "./RUN-genetic-algorithms.sh 1m";
	FILE *fp;

	sprintf(run_buf, "%s %d %d %d %d %d  %d %d %d %d %d  %d %d %d %d %d  %d",
		shell_process, arr[0], arr[1], arr[2],
		arr[3], arr[4], arr[5], arr[6],
		arr[7], arr[8], arr[9], arr[10],
		arr[11], arr[12], arr[13], arr[14],
		arr[15]);
	fp = popen(run_buf, "r");
	printf("run_buf: %s\n", run_buf);

	for (int i = 0; i < lenchrom; i++) {
		fgets(return_buf, sizeof(return_buf), fp);
		printf("return_buf: %s", return_buf);
	}
	fgets(return_buf, sizeof(return_buf), fp);
	int return_number = atoi(return_buf);
	printf("return_number: %d\n", return_number);

	pclose(fp);
	return return_number;
}

int sum(int *fitness)
{
	int sum_fit = 0;
	for (int i = 0; i < sizepop; i++)
		sum_fit += *(fitness + i);
	return sum_fit;
}

int *min(int *fitness)
{
	int min_fit = *fitness;
	int min_index = 0;
	static int arr[2];
	for (int i = 1; i < sizepop; i++) {
		if (*(fitness + i) < min_fit) {
			min_fit = *(fitness + i);
			min_index = i;
		}
	}
	arr[0] = min_index;
	arr[1] = min_fit;
	return arr;
}

void init_chrom()
{
	double pick_number = 0;
	for (int i = 0; i < sizepop; i++) {
		for (int j = 0; j < lenchrom; j++) {
			pick_number = bound_down[j] +
				(rand()) % (bound_up - bound_down);
			chrom[i][j] = (int)(pick_number + 0.5);
		}
		fitness[i] = fit_func(chrom[i]);
	}
}

void Select(int chrom[sizepop][lenchrom])
{
	int index[sizepop];

	int sum_fitness = 0;
	for (int i = 0; i < sizepop; i++) {
		sum_fitness += fitness[i];
	}

	for (int i = 0; i < sizepop; i++) {
		fitness_prob[i] = ((double)fitness[i]) / sum_fitness;
	}

	for (int i = 0; i < sizepop; i++) {
		double pick = ((double)rand()) / RAND_MAX;
		while (pick < 0.0001)
			pick = ((double)rand()) / RAND_MAX;
		for (int j = 0; j < sizepop; j++) {
			pick = pick - fitness_prob[j];
			if (pick <= 0) {
				index[i] = j;
				break;
			}
		}
	}

	int tmp_chrom[sizepop][lenchrom];
	int tmp_fitness[sizepop];
	for (int i = 0; i < sizepop; i++) {
		for (int j = 0; j < lenchrom; j++) {
			tmp_chrom[i][j] = chrom[index[i]][j];
		}
		tmp_fitness[i] = fitness[index[i]];
	}
	for (int i = 0; i < sizepop; i++) {
		for (int j = 0; j < lenchrom; j++) {
			chrom[i][j] = tmp_chrom[i][j];
		}
		fitness[i] = tmp_fitness[i];
	}
}

void Cross(int chrom[sizepop][lenchrom])
{
	for (int i = 0; i < sizepop; i++) {
		double pick1 = ((double)rand()) / RAND_MAX;
		double pick2 = ((double)rand()) / RAND_MAX;
		int choice1 = (int)(pick1 * sizepop);
		int choice2 = (int)(pick2 * sizepop);
		while (choice1 > sizepop - 1) {
			pick1 = ((double)rand()) / RAND_MAX;
			choice1 = (int)(pick1 * sizepop);
		}
		while (choice2 > sizepop - 1) {
			pick2 = ((double)rand()) / RAND_MAX;
			choice2 = (int)(pick2 * sizepop);
		}

		double pick = ((double)rand()) / RAND_MAX;

		if (pick > pcross)
			continue;

		int flag = 0;
		while (flag == 0) {
			double pick = ((double)rand()) / RAND_MAX;
			int pos = (int)(pick * lenchrom);
			while (pos > lenchrom - 1) {
				double pick = ((double)rand()) / RAND_MAX;
				pos = (int)(pick * lenchrom);
			}

			double r = ((double)rand()) / RAND_MAX;
			int v1 = chrom[choice1][pos];
			int v2 = chrom[choice2][pos];
			chrom[choice1][pos] = (int)(r * v2 + (1 - r) * v1);
			chrom[choice2][pos] = (int)(r * v1 + (1 - r) * v2);
			if (chrom[choice1][pos] >= bound_down[pos]
			    && chrom[choice1][pos] <= bound_up[pos]
			    && chrom[choice2][pos] >= bound_down[pos]
			    && chrom[choice2][pos] <= bound_up[pos])
				flag = 1;
		}
	}
}

void Mutation(int chrom[sizepop][lenchrom])
{
	for (int i = 0; i < sizepop; i++) {
		double pick = ((double)rand()) / RAND_MAX;
		int choice = (int)(pick * sizepop);
		while (choice > sizepop - 1) {
			pick = ((double)rand()) / RAND_MAX;
			choice = (int)(pick * sizepop);
		}

		pick = ((double)rand()) / RAND_MAX;
		if (pick > pmutation)
			continue;

		pick = ((double)rand()) / RAND_MAX;
		int pos = (int)(pick * lenchrom);
		while (pos > lenchrom - 1) {
			pick = ((double)rand()) / RAND_MAX;
			pos = (int)(pick * lenchrom);
		}
		int v = chrom[i][pos];
		int v1 = v - bound_up[pos];
		int v2 = bound_down[pos] - v;
		double r = ((double)rand()) / RAND_MAX;
		double r1 = ((double)rand()) / RAND_MAX;
		if (r >= 0.5) {
			chrom[i][pos] = (int)(v - v1 * r1 * (1 - ((double)i) / maxgen)
					* (1 - ((double)i) / maxgen));
		} else {
			chrom[i][pos] = (int)(v + v2 * r1 * (1 - ((double)i) / maxgen)
					* (1 - ((double)i) / maxgen));
		}
	}
}

int main(void)
{
	time_t start, finish;
	start = clock();
	srand((unsigned)time(NULL));
	init_chrom();
	int *best_fit_index = min(fitness);
	int best_index = (int)(*best_fit_index);
	gbest = *(best_fit_index + 1);
	gbest_index = 0;
	average_best[0] = (double)(sum(fitness)) / sizepop;
	for (int i = 0; i < lenchrom; i++)
		gbest_pos[i] = chrom[best_index][i];

	for (int i = 0; i < maxgen; i++) {
		Select(chrom);
		Cross(chrom);
		Mutation(chrom);
		for (int j = 0; j < sizepop; j++) {
			fitness[j] = fit_func(chrom[j]);
		}
		double sum_fit = sum(fitness); // need or not?
		average_best[i + 1] = sum_fit / sizepop;
		int *arr = min(fitness);
		int new_best_index = *arr;
		int new_best = *(arr + 1);
		if (new_best < gbest) {
			gbest = new_best;
			for (int j = 0; j < lenchrom; j++) {
				gbest_pos[j] = chrom[new_best_index][j];
			}
			gbest_index = i + 1;
		}
	}

	finish = clock();	// 程序计算结束
	double duration = ((double)(finish - start)) / CLOCKS_PER_SEC;
	printf("程序计算耗时:%lf秒\n.", duration);
	printf("遗传算法进化了%d次，最优值为:%d,"
		"最优值在第%d代取得,此代的平均最优值为%lf.\n",
	     maxgen, gbest, gbest_index, average_best[gbest_index]);
	printf("取得最优值的地方为(%d,%d,%d,%d,%d, %d,%d,%d,%d,%d, %d,%d,%d,%d,%d, %d).\n",
	       gbest_pos[0], gbest_pos[1], gbest_pos[2], gbest_pos[3],
		gbest_pos[4], gbest_pos[5], gbest_pos[6], gbest_pos[7],
		gbest_pos[8], gbest_pos[9], gbest_pos[10], gbest_pos[11],
		gbest_pos[12], gbest_pos[13], gbest_pos[14], gbest_pos[15]);
	return 0;
}
