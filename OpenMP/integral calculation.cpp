#include <iostream>
#include <omp.h>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono;

double func(double x)   // значение подынтегральной функции 
{
    return (pow(x, 3) - 3);  // например, x^3 - 3
}

double sequential(double sum, double n, double h, double a)
{
    for (int i=1; i<n; i++)
        sum += func(a + h*i);
    return sum;
}

double parallel(double sum, double n, double h, double a, int number_of_threads)
{
    int i;
    #pragma omp parallel private(i) shared(n, h) reduction(+:sum) num_threads(number_of_threads)
    {
        int threads = omp_get_num_threads();
        int id = omp_get_thread_num();
        for (i=id+1; i<n; i+=threads)
            sum += func(a + h*i);
    }
    return sum;
}

int main()
{
    int i;
    double a, b;
    cout << "Введите пределы интегрирования." << endl;
    cout << "Нижний: ";
    cin >> a;
    cout << "Верхний: ";
    cin >> b;

    double n = 10000, sum, h, integral;
    h = (b - a) / n;
    sum = (func(a) + func(b)) / 2; // значения на концах других частичных отрезков считаются дважды

    auto start = high_resolution_clock::now();
    double sum_seq = sequential(sum, n, h, a);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "sequential (function): " << duration.count() << " microseconds" << endl;

    start = high_resolution_clock::now();
    double sum_threads2 = parallel(sum, n, h, a, 2);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "parallel (2 threads, function): " << duration.count() << " microseconds" << endl;

    start = high_resolution_clock::now();
    double sum_threads4 = parallel(sum, n, h, a, 4);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "parallel (4 threads, function): " << duration.count() << " microseconds" << endl;

    start = high_resolution_clock::now();
    double sum_threads8 = parallel(sum, n, h, a, 8);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "parallel (8 threads, function): " << duration.count() << " microseconds" << endl;

    double integral_seq = sum_seq * h;
    double integral_2 = sum_threads2 * h;
    double integral_4 = sum_threads4 * h;
    double integral_8 = sum_threads8 * h;

    cout << "integral (sequential program): " << integral_seq << endl;
    cout << "integral (parallel, 2 threads): " << integral_2 << endl;
    cout << "integral (parallel, 4 threads): " << integral_4 << endl;
    cout << "integral (parallel, 8 threads): " << integral_8 << endl;
}
