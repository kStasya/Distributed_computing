#include <iostream>
#include <omp.h>
using namespace std;

void task1(int number_of_threads) 
{ 
    #pragma omp parallel num_threads(number_of_threads)
    {
        int id = omp_get_thread_num();  // идентификатор
        int count = omp_get_num_threads();  // количество потоков всего
        //cout << "thread number: " << id << ", threads: " << count << ". Hello world!" << endl;  
        printf("Thread number: %d, threads: %d. Hello world!\n", id, count);
    }
}

void task2(int number_of_threads)
{
    int N = 16000;
    int a[N];
    #pragma omp parallel for num_threads(number_of_threads) schedule(guided, 1000) // 1) schedule(static, 1000), 2) schedule(dynamic, 1000), 3) schedule(guided)
    for (int i=0; i<N; i++)
    {
        printf("Thread number: %d, threads: %d, current number: %d.\n", omp_get_thread_num(), omp_get_num_threads(), i);
        a[i] = i;
    }

    printf("\n\n\n================-----------=================\n\n\n\n");  //разграничитель при выводе, для личного удобства
  
    float b[N];
    #pragma omp parallel for num_threads(number_of_threads) schedule(guided, 1000)
    for (int i=0; i<N; i++)
    {
        printf("Thread number: %d, threads: %d, current number: %d.\n", omp_get_thread_num(), omp_get_num_threads(), i);
        if (i == 0)
            b[i] = (a[i+1] + a[i]) / 2.0;
        else if (i == N-1)
            b[i] = (a[i] + a[i-1]) / 2.0;
        else b[i] = (a[i] + a[i-1] + a[i+1]) / 2.0;
    }
}

int main()
{
    //task1(8);
    task2(8);
}
