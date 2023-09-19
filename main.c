#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>


#define elementSize 16

unsigned short *mainMasks;
unsigned short *denseArrayFromCompactMasks;
unsigned short *calculateSparseAndUsual2Masks;
unsigned short *calculateSparseAndUsual2NegMasks;


int my_random_seed = 0;

/**стандартная свертка разреженного и не разреженного векторов F2
 * n - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (F2)*/
void
computationF2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);

/**стандартная свертка разреженного и не разреженного векторов в Z
 * n - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (Z)*/
void
computationZ(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);


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

/**генерирует разреженный массив, как массив индексов ненулевых элементов, упорядоченных в порядке возрастания.
 * arr - результирующий вектор индексов
 * n - максимальный индекс в векторе
 * resSize - количество элементов в результирующем векторе
 * */
void generateSparseArray(unsigned short *arr, int n, unsigned short resSize) {
    unsigned short i = 0;
    int flag;
    while (i < resSize) {
        unsigned int c = rand() % n;

        flag = 1;
        for (int j = 0; j < i; j++) {
            if (arr[i] == c) {
                flag = 0;
                break;
            }
        }
        if (flag)
            arr[i++] = c;
    }
}


/**преоброзует обычный вектор в "компактный"*/
void createDenseArray(unsigned short *res, int resSize, const unsigned short *donor, int donorSize) {
    int i;
    for (i = 0; i < donorSize / elementSize; i++) {
        unsigned short temp = 0;
        for (int j = 0; j < elementSize; j++) {
            temp = temp << 1;
            temp += donor[i * elementSize + j];
        }
        res[i] = temp;
    }
    if (donorSize % elementSize != 0) {

        unsigned short temp = 0;
        for (int j = i * elementSize; j < donorSize; j++) {
            temp = temp << 1;
            temp += donor[j];
        }
        temp = temp << (elementSize - (donorSize % elementSize));
        res[resSize - 1] = temp;
    }

}

void createDenseArrayFromCompact(unsigned short *res, int resSize, const unsigned short *donor, int donorSize) {
    for (int i = 0; i < donorSize; i++) {
        res[donor[i] / elementSize] += denseArrayFromCompactMasks[donor[i] % elementSize];
    }
}

/**генерирует обычный массив (b)*/
void generateDenseArray(unsigned short *arr, int n) {
    unsigned short temp = 0;
    for (int i = 0; i < n; i++) {
        temp = rand() % 20;
        if (temp > 10) {
            arr[i] = 0;
        } else {
            arr[i] = 1;
        }
    }
}


/** целевая функция
 * n - длина векторов, в битах
 * arrSize - количество элементов в inda
 * m - длина в short
 * inda - массив ненулевых индексов (вектор а)
 * b2 - массив short, полученный путем сжатия по 16 элементов в 1 short
 * res - результат - тоже сжатый
 * */
void calculateSparseAndUsual2(unsigned short n, int arrSize, int m, const unsigned short *inda, unsigned short *b2,
                              unsigned short *res) {


    /**modG - длина "хвоста"*/
    unsigned short modG = n % elementSize;
    /**modGNeg - количество бесполезных битов в последнем элементе*/
    unsigned short modGNeg = (elementSize - modG) % elementSize;


    if (inda[0] < elementSize) {
        int i = 0;
        for (i = 0; inda[i] < elementSize; i++) {
            unsigned short mod = inda[i] % elementSize;
            unsigned short modNeg = (elementSize - mod) % elementSize;


            res[0] = res[0] ^ (b2[0] & calculateSparseAndUsual2Masks[mod]) >> mod;

            unsigned short it = 1;
            for (int j = 1; j < m - 1; j++) {
                res[it] = res[it] ^ (((b2[j] & calculateSparseAndUsual2Masks[mod]) >> mod) +
                                     ((b2[j - 1] & calculateSparseAndUsual2NegMasks[mod]) << modNeg));
                it++;
            }

            unsigned short temp2;
            temp2 = ((b2[m - 1] >> mod) + ((b2[m - 2] & calculateSparseAndUsual2NegMasks[mod]) << modNeg));
            res[it] = res[it] ^ temp2;


            unsigned short temp3 =
                    (b2[m - 1] >> modGNeg) + ((b2[m - 2] & calculateSparseAndUsual2NegMasks[modGNeg]) << modG);
            unsigned short temp4 = (temp3 & calculateSparseAndUsual2NegMasks[mod]) << modNeg;
            res[0] = res[0] ^ temp4;
        }




        /**for по элементам, которые больше elementSize*/
        unsigned short lastElementShifted =
                (b2[m - 1] >> modGNeg) + ((b2[m - 2] & calculateSparseAndUsual2NegMasks[modGNeg]) << modG);

        for (; i < arrSize; i++) {
            unsigned short mod = inda[i] % elementSize;
            unsigned short modNeg = (elementSize - mod) % elementSize;
            unsigned short start = inda[i] / elementSize;

            res[start] =
                    res[start] ^ (((b2[0] & calculateSparseAndUsual2Masks[mod]) >> mod) +
                                  ((lastElementShifted & calculateSparseAndUsual2NegMasks[mod]) << modNeg));

            unsigned short it = start + 1;
            int j = 1;
            for (; it < m; it++, j++) {
                res[it] = res[it] ^ (((b2[j] & calculateSparseAndUsual2Masks[mod]) >> mod) +
                                     ((b2[j - 1] & calculateSparseAndUsual2NegMasks[mod]) << modNeg));
            }
            j--;

            short newMod = ((short) mod) - modG;
            int flag = 0;

            if (newMod < 0 || modG == 0) {
                flag++;
                j++;
                newMod = (newMod + elementSize) % elementSize;

            }

            /*printf("start %d\n", start);*/

            unsigned short newModNeg = (elementSize - newMod) % elementSize;

            it = 0;
            for (; j < m - 1; j++, it++) {
                res[it] =
                        res[it] ^
                        (((b2[j] & calculateSparseAndUsual2Masks[newMod]) >> newMod) +
                         ((b2[j - 1] & calculateSparseAndUsual2NegMasks[newMod]) << newModNeg));
            }

            if (flag) {
                res[it] =
                        res[it] ^
                        (((b2[j] & calculateSparseAndUsual2Masks[newMod]) >> newMod) +
                         ((b2[j - 1] & calculateSparseAndUsual2NegMasks[newMod]) << newModNeg));
                it++;
            }
        }


    } else {

        unsigned short lastElementShifted =
                (b2[m - 1] >> modGNeg) + ((b2[m - 2] & calculateSparseAndUsual2NegMasks[modGNeg]) << modG);
        for (int i = 0; i < arrSize; i++) {
            unsigned short mod = inda[i] % elementSize;
            unsigned short modNeg = (elementSize - mod) % elementSize;
            unsigned short start = inda[i] / elementSize;


            res[start] =
                    res[start] ^ (((b2[0] & calculateSparseAndUsual2Masks[mod]) >> mod) +
                                  ((lastElementShifted & calculateSparseAndUsual2NegMasks[mod]) << modNeg));

            unsigned short it = start + 1;
            int j = 1;
            for (; it < m; it++, j++) {
                res[it] = res[it] ^ (((b2[j] & calculateSparseAndUsual2Masks[mod]) >> mod) +
                                     ((b2[j - 1] & calculateSparseAndUsual2NegMasks[mod]) << modNeg));
            }
            j--;

            short newMod = ((short) mod) - modG;
            int flag = 0;

            if (newMod < 0 || modG == 0) {
                flag++;
                j++;
                newMod = (newMod + elementSize) % elementSize;

            }
            unsigned short newModNeg = (elementSize - newMod) % elementSize;

            it = 0;
            for (; j < m - 1; j++, it++) {
                res[it] =
                        res[it] ^
                        (((b2[j] & calculateSparseAndUsual2Masks[newMod]) >> newMod) +
                         ((b2[j - 1] & calculateSparseAndUsual2NegMasks[newMod]) << newModNeg));
            }

            if (flag) {
                res[it] =
                        res[it] ^
                        (((b2[j] & calculateSparseAndUsual2Masks[newMod]) >> newMod) +
                         ((b2[j - 1] & calculateSparseAndUsual2NegMasks[newMod]) << newModNeg));
                it++;
            }
        }
    }


    /* обрезка вектора тут не к чему, но создаст потенциально проблемы в скорости*/
}

/** целевая функция
 * n - длина векторов, в битах
 * arrSize - количество элементов в inda
 * m - длина в short
 * inda - массив ненулевых индексов (вектор а)
 * b2 - массив short, полученный путем сжатия по 16 элементов в 1 short
 * res - результат - тоже сжатый
 * */
void
calculateSparseAndUsual2NewVersion(unsigned short n, int arrSize, int m, const unsigned short *inda, unsigned short *b2,
                                   unsigned short *res) {


    /**modG - длина "хвоста"*/
    unsigned short modG = n % elementSize;
    /**modGNeg - количество бесполезных битов в последнем элементе*/
    unsigned short modGNeg = (elementSize - modG) % elementSize;


    unsigned short lastElementShifted =
            (b2[m - 1] >> modGNeg) + ((b2[m - 2] & calculateSparseAndUsual2NegMasks[modGNeg]) << modG);
    for (int i = 0; i < arrSize; i++) {
        unsigned short mod = inda[i] % elementSize;
        unsigned short modNeg = (elementSize - mod) % elementSize;

        int start = inda[i] / elementSize;
        res[start] =
                res[start] ^ ((b2[0] >> mod) +
                              ((lastElementShifted & calculateSparseAndUsual2NegMasks[mod]) << modNeg));
        int j = 1;
        for (int it = start + 1; it < m; it++, j++) {
            res[it] = res[it] ^ (((b2[j] /*& calculateSparseAndUsual2Masks[mod]*/) >> mod) +
                                 ((b2[j - 1] & calculateSparseAndUsual2NegMasks[mod]) << modNeg));
        }

        /*res[m - 1] = res[m - 1] ^ (((b2[j] & calculateSparseAndUsual2Masks[mod]) >> mod) +
                                   ((b2[j - 1] & calculateSparseAndUsual2NegMasks[mod]) << modNeg));*/

        /*if(inda[0] == 2 && n == 32)
        {
            printf("%s\n",toBinary(lastElementShifted, elementSize));
            for (int t = 0; t < m; t++) {
                printf("%s ", toBinary(b2[t], elementSize));
            }
            printf("\n");
            for (int t = 0; t < m; t++) {
                printf("%s ", toBinary(res[t], elementSize));
            }
            printf("\n");
        }*/

        int targetElement = (start - 1) * elementSize + n - inda[i] + 15;
        int modBack = targetElement % elementSize;
        int backStartTarget = targetElement / elementSize;
        int newMod = elementSize - modBack - 1;
        int newModNeg = (elementSize - newMod) % elementSize;
        j = backStartTarget;
        for (int it = start - 1; it >= 0; it--, j--) {
            res[it] = res[it] ^ (b2[j] >> newMod) ^ ((b2[j - 1] & calculateSparseAndUsual2NegMasks[newMod]) << newModNeg);
        }

    }

    res[m - 1] = res[m - 1] & calculateSparseAndUsual2Masks[modGNeg];


    /* обрезка вектора тут не к чему, но создаст потенциально проблемы в скорости*/
}

/**разжимает вектор*/
void decompressArray(unsigned short n, unsigned short m, unsigned short *arrTo, unsigned short *arrFrom) {
    int it = n - 1;
    unsigned short temp = arrFrom[m - 1];
    unsigned short mod = n % elementSize;
    temp = temp >> (elementSize - (n % elementSize)) % elementSize;

    if (mod == 0) {
        for (int i = m - 1; i >= 0; i--) {
            temp = arrFrom[i];
            for (int j = 0; j < elementSize; j++) {
                arrTo[it--] = temp & 1;
                temp = temp >> 1;
            }
        }
    } else {
        //обработка последнего элемента массива
        for (int i = 0; i < mod; i++) {
            arrTo[it--] = temp & 1;
            temp = temp >> 1;
        }

        //обработка основной части элементов
        for (int i = m - 2; i >= 0; i--) {
            temp = arrFrom[i];
            for (int j = 0; j < elementSize; j++) {
                arrTo[it--] = temp & 1;
                temp = temp >> 1;
            }
        }
    }

}

int compare(unsigned short n, unsigned short *arr1, unsigned short m, unsigned short *arr2) {
    unsigned short arr2Decompressed[n];
    decompressArray(n, m, arr2Decompressed, arr2);
    for (int i = 0; i < n; i++) {
        if (arr1[i] != arr2Decompressed[i]) {
            printf("error in %d\n", i);
            return 0;
        }
    }
    return 1;
}

void main4(long long seed) {
    srand(seed);

    unsigned short epoh = 100;


    /**|e1| + |e2| <= t
     * |e1| = |e2| - не обязательно
     *
     * Для примера
     * n = 4801
     * |h1| = |h2| = 45
     * |e1| + |e2| = 84
     * T = 45/2+4*/
    unsigned short hLength = 45;
    unsigned short eLength = 42;
    unsigned T = hLength / 2 + 5;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 4801;

    unsigned short m = (n + elementSize - 1) / elementSize;

    /**компактное хранение h1*/
    unsigned short h1Compact[hLength];

    /**компактное хранение h2*/
    unsigned short h2Compact[hLength];

    /**компактное хранение h1 в обратном порядке*/
    unsigned short h1TransCompact[hLength];

    /**компактное хранение h2 в обратном порядке*/
    unsigned short h2TransCompact[hLength];



    /**компактное хранение e1*/
    unsigned short e1Compact[eLength];

    /**компактное хранение e2*/
    unsigned short e2Compact[eLength];

    unsigned short e1Compressed[m];
    unsigned short e2Compressed[m];

    unsigned short c1[m];
    unsigned short c2[m];
    unsigned short s[m];

    unsigned short u[m], v[m];

    unsigned short sTemp[m];
    unsigned short upc1[n], upc2[n];
    unsigned short sTempDecompressed[n];


    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {
        clock_t begin = clock();

        for (int i = 0; i < hLength; i++) {
            h1Compact[i] = h2Compact[i] = h1TransCompact[i] = h2TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength; i++) {
            e1Compact[i] = e2Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength);
        generateSparseArray(h2Compact, n, hLength);
        generateSparseArray(e1Compact, n, eLength);
        generateSparseArray(e2Compact, n, eLength);

        for (int i = 0; i < hLength; i++) {
            h1TransCompact[i] = n - h1Compact[i];
            h2TransCompact[i] = n - h2Compact[i];
        }

        for (int i = 0; i < m; i++) {
            e1Compressed[i] = e2Compressed[i] = 0;
        }

        createDenseArrayFromCompact(e1Compressed, m, e1Compact, eLength);
        createDenseArrayFromCompact(e2Compressed, m, e2Compact, eLength);

        for (int i = 0; i < m; i++) {
            c1[i] = c2[i] = s[i] = 0;
        }


        calculateSparseAndUsual2(n, hLength, m, h1Compact, e1Compressed, c1);

        calculateSparseAndUsual2(n, hLength, m, h2Compact, e2Compressed, c2);

        for (int i = 0; i < m; i++) {
            s[i] = c1[i] ^ c2[i];
        }

        /**посчитали S*/



        for (int i = 0; i < m; i++)
            u[i] = v[i] = 0;

        for (int i = 0; i < m; i++) {
            sTemp[i] = s[i];
        }

        for (int z = 0; z < num_it; z++) {
            for (int i = 0; i < n; i++)
                upc1[i] = upc2[i] = 0;

            for (int i = 0; i < n; i++)
                sTempDecompressed[i] = 0;

            decompressArray(n, m, sTempDecompressed, sTemp);

            computationZ(n, hLength, h1TransCompact, sTempDecompressed, upc1);
            computationZ(n, hLength, h2TransCompact, sTempDecompressed, upc2);

            //todo ошибка в подсчете изменения
            for (int j = 0; j < n; j++) {
                if (upc1[j] >= T) {
                    u[j / elementSize] = u[j / elementSize] ^ mainMasks[j % elementSize];
                }

                if (upc2[j] >= T) {
                    v[j / elementSize] = v[j / elementSize] ^ mainMasks[j % elementSize];
                }

            }

            for (int i = 0; i < m; i++) {
                sTemp[i] = c1[i] = c2[i] = 0;
            }


            calculateSparseAndUsual2(n, hLength, m, h1Compact, u, c1);


            calculateSparseAndUsual2(n, hLength, m, h2Compact, v, c2);

            for (int i = 0; i < m; i++) {
                sTemp[i] = s[i] ^ c1[i] ^ c2[i];
            }

            /**проверка s` на ноль*/
            flag = 0;
            exitFlag = 0;
            for (int i = 0; i < m; i++) {
                if (sTemp[i] != 0) {
                    flag = 1;
                    break;
                }
            }
            if (!flag) {
                exitFlag = 1;
                break;
            }

        }

        if (exitFlag) {
            suc++;

        } else {
            fail++;
        }
        clock_t end = clock();
        time_spent += (double) (end - begin) / CLOCKS_PER_SEC;

    }

    printf("last mod\n");
    printf("The elapsed time is %f seconds\n", time_spent);
    printf("One run is %f seconds\n", (time_spent / epoh));
    printf("%f \n", 100.0 * suc / (fail + suc));
}

void test(int nFrom, int nTo) {
    int arrSize = 1;
    unsigned short *inda = calloc(arrSize, sizeof(unsigned short));
    for (int n = nFrom; n <= nTo; n++) {
        int m = (n + elementSize - 1) / elementSize;
        for (int i = 0; i < n; i++) {
            inda[0] = i;

            unsigned short *b = calloc(n, sizeof(unsigned short));
            generateDenseArray(b, n);
            unsigned short *b2 = calloc(m, sizeof(unsigned short));
            createDenseArray(b2, m, b, n); //сжимаем вектор для использование в нашей функции
            unsigned short *res = calloc(n, sizeof(unsigned short));
            for (int j = 0; j < n; j++)
                res[j] = 0;
            unsigned short *res2 = calloc(m, sizeof(unsigned short));
            for (int j = 0; j < m; j++)
                res2[j] = 0;

            computationF2(n, arrSize, inda, b, res);
            calculateSparseAndUsual2NewVersion(n, arrSize, m, inda, b2, res2);

            if (compare(n, res, m, res2) == 0) {
                printf("n %d\n", n);
                printf("inda %d\n", inda[0]);

                printf("b and b2\n");
                for (int j = 0; j < n; j++) {
                    printf("%d", b[j]);
                    if (j % 16 == 15) {
                        printf(" ");
                    }
                }
                printf("\n");

                for (int j = 0; j < m; j++) {
                    printf("%s ", toBinary(b2[j], elementSize));
                }
                printf("\n");

                printf("res and res2\n");
                for (int j = 0; j < n; j++) {
                    printf("%d", res[j]);
                    if (j % 16 == 15) {
                        printf(" ");
                    }
                }
                printf("\n");

                for (int j = 0; j < m; j++) {
                    printf("%s ", toBinary(res2[j], elementSize));
                }
                printf("\n");
            }

            free(b);
            free(b2);
            free(res);
            free(res2);
        }

    }

    free(inda);
}

int main() {
    long long seed = time(0);

    /*начало блока инициализации*/
    mainMasks = calloc(elementSize, sizeof(unsigned short));
    mainMasks[0] = 1 << (elementSize - 1);
    for (int i = 1; i < elementSize; i++)
        mainMasks[i] = mainMasks[i - 1] >> 1;
    denseArrayFromCompactMasks = calloc(elementSize, sizeof(unsigned short));
    denseArrayFromCompactMasks[elementSize - 1] = 1;
    for (int i = 0; i < elementSize - 1; i++) {
        denseArrayFromCompactMasks[elementSize - 2 - i] = denseArrayFromCompactMasks[elementSize - 1 - i] << 1;
    }

    calculateSparseAndUsual2Masks = calloc(elementSize, sizeof(unsigned short));
    calculateSparseAndUsual2NegMasks = calloc(elementSize, sizeof(unsigned short));
    unsigned short temp = -1;
    for (int i = 0; i < elementSize; i++) {
        calculateSparseAndUsual2Masks[i] = temp;
        calculateSparseAndUsual2NegMasks[i] = ~temp;
        temp -= 1 << (i);
    }

    /*printf("mainMasks\n");
    for (int i = 0; i < elementSize; i++) {
        printf("%s\n", toBinary(mainMasks[i], elementSize));
    }
    printf("\n");
    printf("\n");


    printf("calculateSparseAndUsual2Masks\n");
    for (int i = 0; i < elementSize; i++) {
        printf("%s\n", toBinary(calculateSparseAndUsual2Masks[i], elementSize));
    }

    printf("\n");
    printf("\n");


    printf("calculateSparseAndUsual2NegMasks\n");
    for (int i = 0; i < elementSize; i++) {
        printf("%s\n", toBinary(calculateSparseAndUsual2NegMasks[i], elementSize));

    }*/

    /*конец блока инициализации*/

    srand(time(0));

    //ручной тест
    /*int n = 35;
    int m = (n + elementSize - 1) / elementSize;
    unsigned short *b = calloc(35, sizeof(unsigned short));
    for (int i = 0; i < n; i++) {
        b[i] = rand() % 2;
    }
    unsigned short *b2 = calloc(m, sizeof(unsigned short));

    createDenseArray(b2, m, b, n);
    unsigned short *inda = calloc(1, sizeof(unsigned short));
    inda[0] = 18;

    unsigned short *res2 = calloc(m, sizeof(unsigned short));
    for (int i = 0; i < m; i++) {
        res2[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        printf("%d", b[i]);
        if (i % 16 == 15) {
            printf(" ");
        }
    }
    printf("\n");

    for (int i = 0; i < m; i++) {
        printf("%s ", toBinary(b2[i], elementSize));
    }
    printf("\n");

    calculateSparseAndUsual2NewVersion(n, 1, m, inda, b2, res2);
    for (int i = 0; i < m; i++) {
        printf("%s ", toBinary(res2[i], elementSize));
    }
    printf("\n");*/



    test(4801, 4801);



    /* main4(seed);
     //main_glukhikh_optimaz(seed);
     main3(seed);*/
    return 0;
}

void
computationF2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - inda[i] + n) % n]) % 2;
        }
    }
}

void baseComputationF2(unsigned short n, const unsigned short *a, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - i + n) % n]) % 2;
        }
    }
}

void baseComputationZ(unsigned short n, const unsigned short *a, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - i + n) % n]);
        }
    }
}

void
computationZ(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - inda[i] + n) % n]);
        }
    }
}




