#include <iostream>
#include <omp.h>
#include <chrono>
using namespace std;
using namespace std::chrono;

void task4(const int N, const int M, const int K)
{
    // объявление и инициализация двумерного массива
    int arr1[N][M];
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<M; j++)
        {
            arr1[i][j] = i + j + 1;  // для примера
        }
    }

    int arr2[M][K];
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<K; j++)
        {
            arr2[i][j] = i + j + 1;
        }
    }
    
    int result[N][K];
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<K; j++)
        {
            result[i][j] = 0;
        }
    }

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<M; j++)
        {
            for (int k=0; k<K; k++)
            {
                result[i][k] += arr1[i][j] * arr2[j][k];
            }
            
        }
    }

    /* //проверка результата
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<K; j++)
        {
            cout << result[i][j] << "\t";
        }
        cout << endl;
    }
    */


}

void task4_parr(const int N, const int M, const int K, int number_of_threads)
{
    int i, j, k;
    // объявление и инициализация двумерного массива
    int arr1[N][M];
    //#pragma omp parallel for num_threads(number_of_threads)
    #pragma omp parallel for private(i, j) collapse(2) num_threads(number_of_threads) 
    for (i=0; i<N; i++)
    {
        for (j=0; j<M; j++)
        {
            arr1[i][j] = i + j + 1;  // для примера
        }
    }

    int arr2[M][K];
    //#pragma omp parallel for num_threads(number_of_threads)
    #pragma omp parallel for private(i, j) collapse(2) num_threads(number_of_threads)
    for (i=0; i<M; i++)
    {
        for (j=0; j<K; j++)
        {
            arr2[i][j] = i + j + 1;
        }
    }
    
    int result[N][K];
    //#pragma omp parallel for num_threads(number_of_threads)
    #pragma omp parallel for private(i, j) collapse(2) num_threads(number_of_threads)
    for (i=0; i<N; i++)
    {
        for (j=0; j<K; j++)
        {
            result[i][j] = 0;
        }
    }
    //#pragma omp parallel for num_threads(number_of_threads)
    #pragma omp parallel for collapse(3) private(i, j, k) num_threads(number_of_threads)
    for (i=0; i<N; i++)
    {
        for (j=0; j<M; j++)
        {
            for (k=0; k<K; k++)
            {
                result[i][k] += arr1[i][j] * arr2[j][k];
            }
            
        }
    }
}

int main()
{
    const int a = 15;
    const int b = 25;
    const int c = 17;
    cout << "Size of matrices: (" << a << ", " << b << "), (" << b << ", " << c << ")" << endl;

    auto start = high_resolution_clock::now();
    task4(a, b, c); 
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Sequential: " << duration.count() << " microseconds" << endl;

    start = high_resolution_clock::now();
    task4_parr(a, b, c, 2);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Parallel (2 threads): " << duration.count() << " microseconds" << endl;

    start = high_resolution_clock::now();
    task4_parr(a, b, c, 4);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Parallel (4 threads): " << duration.count() << " microseconds" << endl;

    start = high_resolution_clock::now();
    task4_parr(a, b, c, 8);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Parallel (8 threads): " << duration.count() << " microseconds" << endl;

}
