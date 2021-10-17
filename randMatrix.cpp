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

namespace matrixOperations{
    
    inline void sumMatrix(float **A, float **B, float **C, int size_l, int size_m){
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                C[i][j] = A[i][j] + B[i][j];
            }
        }
    }

    inline void differenceMatrix(float **A, float **B, float **C, int size_l, int size_m){
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                C[i][j] = A[i][j] - B[i][j];
            }
        }
    }


    inline float frobeniusNorm(float **A, int size_l, int size_m){
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                norm += A[i][j] * A[i][j];
            }
        }
        
        return sqrt(norm);
    }


    float euclideanNormRows(float **A, int index, int size_l){
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            norm += A[i][index] * A[i][index];
        }
        return sqrt(norm);
    }


    float euclideanNormCols(float **array, int index, int size_l){
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            norm += array[index][i] * array[index][i];
        }
        return sqrt(norm);
    }


    
    
    void subsampleMatrix(float **inputMatrix1, float **inputMatrix2, float **outputMatrix1, float **outputMatrix2, int l, int m, int n, int red_m, int sampling){
        
        if(sampling == SAMPLING_STRATEGY){
            
            float sum = 0.0;
            float *A_row = new float[m];
            float *B_col = new float[m];
            int *value = new int[red_m];
            
            for(int i = 0; i < m; i++){
                A_row[i] = euclideanNormRows(inputMatrix1, i, l);
                B_col[i] = euclideanNormCols(inputMatrix2, i, n);
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


    void printMatrix(float **matrix, int size_l, int size_m){
        
        cout << endl;
        for(int i = 0; i < size_l; i++){
            cout << setw(MATRIX_ENTRY_SPACE);
            for (int j = 0; j < size_m; j++){
                cout << matrix[i][j] << setw(MATRIX_ENTRY_SPACE);
            }
            cout << endl;
        }
    }


    void fillMatrix(float **matrix, int size_l, int size_m){
        
        for(int i = 0; i < size_l; i++){
            for (int j = 0; j < size_m; j++){
                matrix[i][j] = double(rand()) / (RAND_MAX);
            }
        }
    }
}

int main(){
    

    bool consoleMode;
    cout << "Row of experiments (type: 0) or singular experiment (type: 1)?" << endl;
    cin >> consoleMode;
    cout << "------------------------------" << endl;
    
    double elapsed_secs_inner_exact, elapsed_secs_outer_exact, elapsed_secs_outer_approx;
    float error_outer_exact, error_outer_approx;
    
    srand(time(NULL));
    
    if(consoleMode){
        
        int l, m, n, red_m, subsampling;
        float subsample;
        cout << "Rowcount of matrix A: ";
        cin >> l;
        cout << "Row resp. columncount of matrix A resp. B: ";
        cin >> m;
        cout << "Columncount of matrix B: ";
        cin >> n;
        cout << "Factor of subsampling: ";
        cin >> subsample;
        red_m = ceil(m*subsample);
        cout << "Using " << red_m << " rows/cols in the probabilistic product instead of " << m << "." << endl;
        cout << "Use uniform sampling (1) or custom sampling (2): ";
        cin >> subsampling;
        cout << "------------------------------" << endl;
        
        
        float **A, **B, **A_red, **B_red, **C_exact_inner, **C_exact_outer, **C_approx_outer;
        
        A = new float *[l];
        B = new float *[m];
        
        A_red = new float *[l];
        B_red = new float *[red_m];
        
        C_exact_inner = new float *[l];
        C_exact_outer = new float *[l];
        C_approx_outer = new float *[l];
        
        for(int i = 0; i < l; i++){
            A[i] = new float[m];
            A_red[i] = new float[red_m];
            C_exact_inner[i] = new float[n];
            C_exact_outer[i] = new float[n];
            C_approx_outer[i] = new float[n];
        }
        for(int i = 0; i <m; i++){
            B[i] = new float[n];
        }
        for(int i = 0; i < red_m; i++){
            B_red[i] = new float[n];
        }
        
        matrixOperations::fillMatrix(A, l, m);
        matrixOperations::fillMatrix(B, m, n);
        
        if(l <= MAX_DISPLAY_SIZE && m <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            cout << "A = ";
            matrixOperations::printMatrix(A, l, m);
            cout << endl;
            cout << "B = ";
            matrixOperations::printMatrix(B, m, n);
            cout << endl;
        }
        else{
            cout << "A = " << endl;
            cout << "To large to display." << endl;
            cout << "B = " << endl;
            cout << "To large to display." << endl;
        }
        cout << "------------------------------" << endl;
        
        clock_t begin_inner_exact = clock();
        matrixOperations::innerMatrixProduct(A, B, C_exact_inner, l, m, n);
        clock_t end_inner_exact = clock();
        elapsed_secs_inner_exact = (double(end_inner_exact - begin_inner_exact) / CLOCKS_PER_SEC) * 1000.0;
        cout << "INNER MATRIX PRODUCT" << endl << endl;
        cout << "A*B = ";
        if(l <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            matrixOperations::printMatrix(C_exact_inner, l, n);
            cout << endl;
        }
        else{
            cout << "To large to display." << endl;
        }
        cout << "Computed the exact inner matrix product in " << elapsed_secs_inner_exact << " milliseconds" << endl;
        cout << "------------------------------" << endl;
        
        clock_t begin_outer_exact = clock();
        matrixOperations::outerMatrixProduct(A, B, C_exact_outer, l, m, n);
        clock_t end_outer_exact = clock();
        elapsed_secs_outer_exact = (double(end_outer_exact - begin_outer_exact) / CLOCKS_PER_SEC) * 1000.0;
        cout << "OUTER MATRIX PRODUCT" << endl << endl;
        cout << "A*B = ";
        if(l <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            matrixOperations::printMatrix(C_exact_outer, l, n);
            cout << endl;
        }
        else{
            cout << "To large to display." << endl;
        }
        cout << "Computed the exact product in " << elapsed_secs_outer_exact << " milliseconds" << endl;
        matrixOperations::differenceMatrix(C_exact_inner, C_exact_outer, C_exact_outer, l, n);
        error_outer_exact = matrixOperations::frobeniusNorm(C_exact_outer, l, n);
        cout << "Error of the exact outer matrix product: " << error_outer_exact << ". Surprise ;-)" << endl;
        cout << "------------------------------" << endl;
        
        clock_t begin_outer_approx = clock();
        matrixOperations::subsampleMatrix(A, B, A_red, B_red, l, m, n, red_m, subsampling);
        matrixOperations::outerMatrixProduct(A_red, B_red, C_approx_outer, l, red_m, n);
        clock_t end_outer_approx = clock();
        elapsed_secs_outer_approx = (double(end_outer_approx - begin_outer_approx) / CLOCKS_PER_SEC) * 1000.0;
        cout << "APPROXIMATE OUTER MATRIX PRODUCT" << endl << endl;
        cout << "C*R = ";
        if(l <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            matrixOperations::printMatrix(C_approx_outer, l, n);
            cout << endl;
        }
        else{
            cout << "To large to display." << endl;
        }
        cout << "Computed the approximate outer matrix product in " << elapsed_secs_outer_approx << " milliseconds" << endl;
        matrixOperations::differenceMatrix(C_exact_inner, C_approx_outer, C_approx_outer, l, n);
        error_outer_approx = matrixOperations::frobeniusNorm(C_approx_outer, l, n);
        cout << "Error of the approximate outer product: " << error_outer_approx << "." << endl;
        cout << "------------------------------" << endl;
        
        for(int i = 0; i <l; i++){
            delete[] A[i];
            delete[] A_red[i];
            delete[] C_exact_inner[i];
            delete[] C_exact_outer[i];
            delete[] C_approx_outer[i];
        }
        for(int i = 0; i < m; i++){
            delete[] B[i];
        }
        for(int i = 0; i < red_m; i++){
            delete[] B_red[i];
        }
        delete[] A;
        delete[] B;
        delete[] A_red;
        delete[] B_red;
        delete[] C_exact_inner;
        delete[] C_exact_outer;
        delete[] C_approx_outer;
    }
    
    else{
    
        int a, b, c, l, m, n, red_m, subsampling;
        float subsample;
        
        cout << "We consider squared matrices." << endl;
        cout << "Set start dimension: ";
        cin >> a;
        cout << "Set end dimension: ";
        cin >> b;
        cout << "Set increment: ";
        cin >> c;
        cout << "Factor of subsampling: ";
        cin >> subsample;
        cout << endl;
        cout << "Use uniform sampling (1) or custom sampling (2): ";
        cin >> subsampling;
        cout << "------------------------------" << endl;
        
        ofstream result;
        result.open ("result.txt");
        result << "matrixDimension,timeInnerExact,errorInnerExact,timeOuterExact,errorOuterExact,timeOuterApprox,errorOuterApprox,subsample,subsampling\n";
        
        for(int k = a; k <= b; k += c){
            
            l = k;
            m = k;
            n = k;
            
            red_m = ceil(m*subsample);
            
            float **A, **B, **A_red, **B_red, **C_exact_inner, **C_exact_outer, **C_approx_outer;
            
            A = new float *[l];
            B = new float *[m];
            
            A_red = new float *[l];
            B_red = new float *[red_m];
            
            C_exact_inner = new float *[l];
            C_exact_outer = new float *[l];
            C_approx_outer = new float *[l];
            
            for(int i_1 = 0; i_1 < l; i_1++){
                A[i_1] = new float[m];
                A_red[i_1] = new float[red_m];
                C_exact_inner[i_1] = new float[n];
                C_exact_outer[i_1] = new float[n];
                C_approx_outer[i_1] = new float[n];
            }
            for(int i_2 = 0; i_2 <m; i_2++){
                B[i_2] = new float[n];
            }
            for(int i_3 = 0; i_3 < red_m; i_3++){
                B_red[i_3] = new float[n];
            }
            
            matrixOperations::fillMatrix(A, l, m);
            matrixOperations::fillMatrix(B, m, n);
            
            clock_t begin_inner_exact = clock();
            matrixOperations::innerMatrixProduct(A, B, C_exact_inner, l, m, n);
            clock_t end_inner_exact = clock();
            elapsed_secs_inner_exact = (double(end_inner_exact - begin_inner_exact) / CLOCKS_PER_SEC) * 1000.0;
            
            clock_t begin_outer_exact = clock();
            matrixOperations::outerMatrixProduct(A, B, C_exact_outer, l, m, n);
            clock_t end_outer_exact = clock();
            elapsed_secs_outer_exact = (double(end_outer_exact - begin_outer_exact) / CLOCKS_PER_SEC) * 1000.0;
            matrixOperations::differenceMatrix(C_exact_inner, C_exact_outer, C_exact_outer, l, n);
            error_outer_exact = matrixOperations::frobeniusNorm(C_exact_outer, l, n);
            
            clock_t begin_outer_approx = clock();
            matrixOperations::subsampleMatrix(A, B, A_red, B_red, l, m, n, red_m, subsampling);
            matrixOperations::outerMatrixProduct(A_red, B_red, C_approx_outer, l, red_m, n);
            clock_t end_outer_approx = clock();
            elapsed_secs_outer_approx = (double(end_outer_approx - begin_outer_approx) / CLOCKS_PER_SEC) * 1000.0;
            matrixOperations::differenceMatrix(C_exact_inner, C_approx_outer, C_approx_outer, l, n);
            error_outer_approx = matrixOperations::frobeniusNorm(C_approx_outer, l, n);
            
            
            result << fixed << l << "," << elapsed_secs_inner_exact << "," << 0.0 << "," << elapsed_secs_outer_exact << "," << error_outer_exact << "," << elapsed_secs_outer_approx << "," << double(error_outer_approx) << "," << subsample << "," << subsampling << "\n";
            
            
        }
        
        result.close();
        
    }
	return 0;

}
