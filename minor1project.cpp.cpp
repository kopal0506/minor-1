#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <ctime>
#include <random>
#include <map>

#define MAX_DISPLAY_SIZE 10
#define MATRIX_ENTRY_SPACE 12
#define SAMPLING_STRATEGY 2

using namespace std;

namespace matrix_Operations{
    
    inline void sum_Matrix(float **A, float **B, float **C, int size_l, int size_m){
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                C[i][j] = A[i][j] + B[i][j];
            }
        }
    }

    inline void difference_Matrix(float **A, float **B, float **C, int size_l, int size_m){
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                C[i][j] = A[i][j] - B[i][j];
            }
        }
    }


    inline float frobenius_Norm(float **A, int size_l, int size_m){
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                norm += A[i][j] * A[i][j];
            }
        }
        
        return sqrt(norm);
    }


    float euclideanNorm_Rows(float **A, int index, int size_l){
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            norm += A[i][index] * A[i][index];
        }
        return sqrt(norm);
    }


    float euclideanNorm_Cols(float **array, int index, int size_l){
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            norm += array[index][i] * array[index][i];
        }
        return sqrt(norm);
    }


    
    
    void subsample_Matrix(float **inputMatrix1, float **inputMatrix2, float **outputMatrix1, float **outputMatrix2, int l, int m, int n, int red_m, int sampling){
        
        if(sampling == SAMPLING_STRATEGY){
            
            float sum = 0.0;
            float *A_row = new float[m];
            float *B_col = new float[m];
            int *value = new int[red_m];
            
            for(int i = 0; i < m; i++){
                A_row[i] = euclideanNorm_Rows(inputMatrix1, i, l);
                B_col[i] = euclideanNorm_Cols(inputMatrix2, i, n);
                sum += (A_row[i] * B_col[i]);
            }
            for(int i = 0; i < m; i++){
                A_row[i] = (A_row[i] * B_col[i]) / sum;
            }
            
            random_device rd;
            mt19937 eng(rd());
            discrete_distribution<> distr(A_row[0], A_row[m]);
            
            for(int i = 0; i < red_m; i++){
                value[i] = distr(eng);
                for(int j = 0; j < n; j++){
                    outputMatrix2[value[i]][j] = inputMatrix2[value[i]][j] / (A_row[value[i]] * red_m);
                }
            }
            for(int i = 0; i < l; i++){
                for(int j = 0; j < red_m; j++){
                    outputMatrix1[i][value[j]] = inputMatrix1[i][value[j]];
                }
            }
            
            delete[] A_row;
            delete[] B_col;
            delete[] value;

        }
        else {
            
            int *value = new int[red_m];
            
            random_device rd;
            mt19937 eng(rd());
            uniform_int_distribution<> distr(0, red_m - 1);
            
            for(int i = 0; i < red_m; i++){
                value[i] = int(distr(eng));
                for(int j = 0; j < n; j++){
                    outputMatrix2[value[i]][j] = inputMatrix2[value[i]][j];
                }
            }
            for(int i = 0; i < l; i++){
                for(int j = 0; j < red_m; j++){
                    outputMatrix1[i][value[j]] = inputMatrix1[i][value[j]];
                }
            }
            delete[] value;
        }
    }


    void outerMatrixProduct(float **A, float **B, float **C, int size_l, int size_m, int size_n){
        
        for(int k = 0; k < size_m; k++){
            for(int i = 0; i < size_l; i++){
                for(int j = 0; j < size_n; j++){
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }


    void innerMatrixProduct(float **A, float **B, float **C, int size_l, int size_m, int size_n){
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_n; j++){
                C[i][j] = 0.0;
                for(int k = 0; k < size_m; k++){
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }


    void print_Matrix(float **matrix, int size_l, int size_m){
        
        cout << endl;
        for(int i = 0; i < size_l; i++){
            cout << setw(MATRIX_ENTRY_SPACE);
            for (int j = 0; j < size_m; j++){
                cout << matrix[i][j] << setw(MATRIX_ENTRY_SPACE);
            }
            cout << endl;
        }
    }


    void fill_Matrix(float **matrix, int size_l, int size_m){
        
        for(int i = 0; i < size_l; i++){
            for (int j = 0; j < size_m; j++){
                matrix[i][j] = double(rand()) / (RAND_MAX);
            }
        }
    }
}

