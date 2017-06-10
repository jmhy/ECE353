/*  Purpose: Calculate definite integral using the trapezoidal rule.
 *
 * Input:   a, b, n
 * Output:  Estimate of integral from a to b of f(x)
 *          using n trapezoids.
 *
 * Author: Naga Kandasamy, Michael Lui
 * Date: 6/22/2016
 *
 * Compile: gcc -o trap trap.c -lpthread -lm
 * Usage:   ./trap
 *
 *
 * Note:    The function f(x) is hardwired.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <pthread.h>

#define LEFT_ENDPOINT 5
#define RIGHT_ENDPOINT 1000
#define NUM_TRAPEZOIDS 100000000
#define NUM_THREADs 16 /* Number of threads to run. */


/*------------------------------------------------------------------
 * Function:    func
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 * Output: (x+1)/sqrt(x*x + x + 1)
 */
__attribute__((const)) float func(float x) 
{
	return (x + 1)/sqrt(x*x + x + 1);
}  


double compute_gold(float, float, int, float (*)(float));
double compute_using_pthreads(float, float, int, float (*)(float));

int main(int argc, char *argv[])
{
	int n = NUM_TRAPEZOIDS;
	float a = LEFT_ENDPOINT;
	float b = RIGHT_ENDPOINT;
	
	clock_t t;
	double time;
	
	t = clock();
	double reference = compute_gold(a, b, n, func);
	t = clock() - t;
	time = ((double)t)/CLOCKS_PER_SEC;
	printf("Reference solution computed on the CPU = %f \n", reference);
	printf("Reference took %f seconds\n", time);

	t = clock();
	double pthread_result = compute_using_pthreads(a, b, n, func); /* Write this function using pthreads. */
	t = clock() - t;
	time = ((double)t)/CLOCKS_PER_SEC;
	printf("Solution computed using pthreads = %f \n", pthread_result);
	printf("Solution took %f seconds\n", time);

} 

/*------------------------------------------------------------------
 * Function:    Trap
 * Purpose:     Estimate integral from a to b of f using trap rule and
 *              n trapezoids
 * Input args:  a, b, n, f
 * Return val:  Estimate of the integral 
 */
double compute_gold(float a, float b, int n, float(*f)(float))
{
	float h = (b-a)/(float)n; /* 'Height' of each trapezoid. */

	double integral = (f(a) + f(b))/2.0;
	
	for (int k = 1; k <= n-1; k++) 
		integral += f(a+k*h);
	
	integral = integral*h;
	
	return integral;
}

struct trap_data
{
	float a;
	float b;
	int n;
	float(*f)(float);
	double area;
};

void *pass_args(void* args)
{
	struct trap_data act_args = *((struct trap_data*)args);
	((struct trap_data*)args)->area = compute_gold(act_args.a, act_args.b, act_args.n, act_args.f);
	return NULL;
}

double compute_using_pthreads(float a, float b, int n, float(*f)(float))
{
	struct trap_data areas[NUM_THREADs];
	pthread_t threads[NUM_THREADs];
	float w = (b-a) / NUM_THREADs;
	
	for (int i = 0; i < NUM_THREADs; i++)
	{
		areas[i].a = a + i*w;
		areas[i].b = a + (i+1)*w;
		areas[i].n = n / NUM_THREADs;
		areas[i].f = f;
		areas[i].area = 0;
		pthread_create(&threads[i], NULL, pass_args, &areas[i]);
	}
	
	double total_area = 0;
	for (int i = 0; i < NUM_THREADs; i++)
	{
		pthread_join(threads[i], NULL);
		total_area = total_area + areas[i].area;
	}
	return total_area;
}
