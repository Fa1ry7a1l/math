#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>


#define elementSize 16

/**перевод последовательности unsigned short в строку*/
char *toBinary(unsigned short n, int len) {
    char *binary = (char *) malloc(sizeof(char) * len);
    int k = 0;
    for (unsigned i = (1 << (len - 1)); i > 0; i = i / 2) {
        binary[k++] = (n & i) ? '1' : '0';
    }
    binary[k] = '\0';
    return binary;
}

/**генерирует разреженный массив*/
void generateSparseArray(unsigned short *arr, int m, unsigned short resSize) {
    unsigned short i = 0;
    unsigned short a = 0;
    while (i < resSize && a < m) {
        unsigned int c = rand();
        a = a + (c % ((m - a) / (resSize - i)));
        arr[i++] = a++;
    }
}


/**преоброзует обычный вектор в "компактный"*/
void createDenseArray(unsigned short **arr, int m, const unsigned short *donor, int n) {

    //Создаем маски для получения последнего элемента
    unsigned short mask[elementSize];
    unsigned short s = 1;
    for (int i = 0; i < elementSize; i++) {
        mask[i] = s;
        s = s << 1;
    }

    /* printf("masks\n");
     for (int i = 0; i < elementSize; i++) {
         printf("%d \n", mask[i]);
     }*/

    //заполнение первого элемента
    unsigned short temp = 0;
    for (int i = 0, j = 0; i < n; i++) {
        if (donor[(n - i) % n] == 1)
            temp += mask[elementSize - 1 - (i % elementSize)];
        if ((i + 1) % elementSize == 0) {
            arr[0][j++] = temp;
            temp = 0;
        }
    }
    //поправляем окончание первой строки, если длина донора была не кратна elementSize
    int mod = n % elementSize;
    if (mod != 0)
        arr[0][m - 1] = temp;

    //печать первой строки
    /*printf("first switch\n");
    for (int i = 0; i < m; i++)
        printf("%s", toBinary(arr[0][i], elementSize));
    printf("\n");*/


    //создание маски для поиска последнего элемента
    unsigned short lastElementMask = -1;
    unsigned short lastBitMask = 1;
    lastBitMask = lastBitMask << (elementSize - mod);
    if (lastBitMask == 0)
        lastBitMask = 1;
    lastElementMask = lastElementMask >> mod;

    //printf("last element mask \n %s\n", toBinary(~lastElementMask, elementSize));


    for (int k = 1; k < elementSize; k++) {

        arr[k][m - 1] = (arr[k - 1][m - 1] >> 1);
        if (arr[k - 1][m - 2] % 2 == 1)
            arr[k][m - 1] += mask[elementSize - 1];
        //printf("last element  \n%s\n %s\n", toBinary(arr[k][m - 1], elementSize), toBinary(arr[k][m - 1] & ~lastElementMask,elementSize));
        arr[k][m - 1] = arr[k][m - 1] & ~lastElementMask;
        for (int i = m - 2; i > 0; i--) {
            arr[k][i] = arr[k - 1][i] >> 1;
            if (arr[k - 1][i - 1] % 2 == 1)
                arr[k][i] += mask[elementSize - 1];

        }
        arr[k][0] = arr[k - 1][0] >> 1;
        if ((arr[k - 1][m - 1] & lastBitMask) > 0)
            arr[k][0] += mask[elementSize - 1];
    }

    for (int i = 0; i < elementSize; i++) {
        for (int j = 0; j < m; j++) {
            printf("%s", toBinary(arr[i][j], elementSize));
        }
        printf("\n");
    }

}

/**генерирует обычный массив (b)*/
void generateDenseArray(unsigned short *arr, int m) {
    unsigned short temp = 0;
    for (int i = 0; i < m; i++) {
        temp = rand() % 20;
        if (temp > 10) {
            arr[i] = 0;
        } else {
            arr[i] = 1;
        }
    }
}

int main() {
    srand(time(0));

    /**M - количество элементов в векторе A, B , C*/
    unsigned short M = 16;
    //unsigned short M = 25546;
    /**arrSize - количество элементов в сжатом векторе A, где хранятся только номера ненулевых эементов разреженного вектора A*/
    int arrSize = 1;
    //int arrSize = 159;
    /**m - размер сжатого вектора B, C*/
    int m = (M + elementSize - 1) / elementSize;
    /**inda - сжатый вектор A, A - разреженный*/
    unsigned short inda[arrSize];
    inda[0] = 0;
    /**b - Вектор с которым происходит свертка*/
    unsigned short b[M];

    /**результирующий вевктор стандартного примера (размера M)*/
    unsigned short c[M];
    for (int i = 0; i < M; i++) {
        c[i] = 0;
    }

    /**результирующий сжатый вектор второго способа (размера m)*/
    unsigned short c2[m];
    for (int i = 0; i < m; i++) {
        c2[i] = 0;
    }


    /**сжатая версия b*/
    unsigned short **b2 = calloc(elementSize, sizeof(unsigned short *));

    for (int i = 0; i < elementSize; i++)
        b2[i] = calloc(m, sizeof(unsigned short));

    // generateSparseArray(inda, M, arrSize);

    generateDenseArray(b, M);
    createDenseArray(b2, m, b, M);


    printf("inda\n");
    for (int i = 0; i < arrSize; i++) {
        printf("%d ", inda[i]);
    }
    printf("\n");

    printf("b\n");
    for (int i = 0; i < M; i++) {
        printf("%d ", b[i]);
    }
    printf("\n");

    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < M; j++) {
            c[j] = (c[j] + b[(inda[i] - j + M) % M]) % 2;
        }
        /* printf("res\n");
         for (int j = 0; j < M; j++) {
             printf("%d", c[j]);
         }
         printf("\n");*/
    }

    /*printf("res\n");
    for (int i = 0; i < M; i++) {
        printf("%d ", c[i]);
    }
    printf("\n");*/

    /**метод для M не кратной elemntSize*/
    if (M % elementSize != 0) {
        /**маска для подсчета последнего элемента*/
        unsigned short lastElementMask = -1;
        lastElementMask = (M % elementSize) == 0 ? -1 : lastElementMask << (elementSize - (M % elementSize));
        unsigned short firstElementMask = !lastElementMask;
        printf("%s", toBinary(lastElementMask, elementSize));

        //todo добавить подсчет маски конца

        for (int i = 0; i < arrSize; i++) {
            /**номер строки определяемый остатком деления сдвига на 16*/
            int rowNumber = inda[i] % elementSize;

            /**номер стартового столбца*/
            int startColumn = inda[i] / elementSize;
            int el = 0;
            /**обрабатываем все элементы со сдвигом так, чтобы */
            for (int j = startColumn; j < m - 1; j++) {
                c2[el] = c2[el] ^ b2[rowNumber][j];
                el++;
            }

            /**обработка последнего элемента*/
            c2[el] = c2[el] ^ (lastElementMask & b2[rowNumber][m - 1]);


        }
    }


    return 0;
}