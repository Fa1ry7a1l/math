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

void baseComputationF2(unsigned short n, const unsigned short *a, const unsigned short *b, unsigned short *c);

void baseComputationZ(unsigned short n, const unsigned short *a, const unsigned short *b, unsigned short *c);

void
my_computationF2_Fig2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                      unsigned short *c);


void
my_computationF3_Fig2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                      unsigned short *c);


void
my_computationZ_Fig2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                     unsigned short *c);

/**стандартная свертка разреженного и не разреженного векторов F2
 * n - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (F2)*/
void
computationF2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);

void
my_computationF2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                 unsigned short *c);

/**стандартная свертка разреженного и не разреженного векторов в Z
 * n - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (Z)*/
void
computationZ(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);

void
my_computationZ_last(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                     unsigned short *c);


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

int compare_f(const void *a, const void *b) {
    int int_a = *((unsigned short *) a);
    int int_b = *((unsigned short *) b);

    if (int_a == int_b) return 0;
    else if (int_a < int_b) return -1;
    else return 1;
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
            if (arr[j] == c) {
                flag = 0;
                break;
            }
        }
        if (flag)
            arr[i++] = c;
    }

    qsort(arr, resSize, sizeof(unsigned short), compare_f);
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
        res[donor[i] / elementSize] = res[donor[i] / elementSize] ^ denseArrayFromCompactMasks[donor[i] % elementSize];
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
void
calculateSparseAndUsual2NewVersion(unsigned short n, int arrSize, int m, const unsigned short *inda,
                                   const unsigned short *b2,
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
            res[it] =
                    res[it] ^ (b2[j] >> newMod) ^ ((b2[j - 1] & calculateSparseAndUsual2NegMasks[newMod]) << newModNeg);
        }

    }

    res[m - 1] = res[m - 1] & calculateSparseAndUsual2Masks[modGNeg];


    /* обрезка вектора тут не к чему, но создаст потенциально проблемы в скорости*/
}

/**разжимает вектор*/
void decompressArray(unsigned short n, unsigned short m, unsigned short *arrTo, const unsigned short *arrFrom) {
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

int compare(unsigned short n, unsigned short *arr1, unsigned short m, const unsigned short *arr2) {
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


void testDecoder(long long seed, int epoh) {
    srand(seed);


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
    //unsigned short n = 50;

    unsigned short m = (n + elementSize - 1) / elementSize;

    /**компактное хранение h1*/
    unsigned short *h1Compact = calloc(hLength, sizeof(unsigned short));

    /**компактное хранение h2*/
    unsigned short *h2Compact = calloc(hLength, sizeof(unsigned short));

    /**компактное хранение h1 в обратном порядке*/
    unsigned short *h1TransCompact = calloc(hLength, sizeof(unsigned short));

    /**компактное хранение h2 в обратном порядке*/
    unsigned short *h2TransCompact = calloc(hLength, sizeof(unsigned short));



    /**компактное хранение e1*/
    unsigned short *e1Compact = calloc(eLength, sizeof(unsigned short));

    /**компактное хранение e2*/
    unsigned short *e2Compact = calloc(eLength, sizeof(unsigned short));

    unsigned short *e1Compressed = calloc(m, sizeof(unsigned short));
    unsigned short *e2Compressed = calloc(m, sizeof(unsigned short));

    unsigned short *c1 = calloc(m, sizeof(unsigned short));
    unsigned short *c2 = calloc(m, sizeof(unsigned short));
    unsigned short *s = calloc(m, sizeof(unsigned short));

    unsigned short *u = calloc(m, sizeof(unsigned short));
    unsigned short *v = calloc(m, sizeof(unsigned short));;

    unsigned short *sTemp = calloc(m, sizeof(unsigned short));
    unsigned short *upc1 = calloc(n, sizeof(unsigned short));
    unsigned short *upc2 = calloc(n, sizeof(unsigned short));
    unsigned short *sTempDecompressed = calloc(n, sizeof(unsigned short));;


    //переменные для тестирования
    unsigned short *e1Temp = calloc(n, sizeof(unsigned short));
    unsigned short *e2Temp = calloc(n, sizeof(unsigned short));
    unsigned short *vTemp = calloc(n, sizeof(unsigned short));;
    unsigned short *uTemp = calloc(n, sizeof(unsigned short));
    unsigned short *c1Test = calloc(n, sizeof(unsigned short));
    unsigned short *c2Test = calloc(n, sizeof(unsigned short));
    unsigned short *sTest = calloc(n, sizeof(unsigned short));
    unsigned short *sTestTemp = calloc(n, sizeof(unsigned short));
    unsigned short *upc1Test = calloc(n, sizeof(unsigned short));
    unsigned short *upc2Test = calloc(n, sizeof(unsigned short));
    int flagTemp;

    int totalIt = 0;


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

        for (int i = 0; i < n; i++) {
            e1Temp[i] = e2Temp[i] = 0;
        }

        /**Проверка правильности заполнения данными*/
        for (int i = 0; i < eLength; i++) {
            e1Temp[e1Compact[i]] = 1;
            e2Temp[e2Compact[i]] = 1;
        }

        if (!compare(n, e1Temp, m, e1Compressed)) {
            printf("ERRRRROR\n");

            for (int i = 0; i < eLength; i++) {
                printf("%d ", e1Compact[i]);
            }
            printf("\n");


            for (int i = 0; i < n; i++) {
                printf("%d", e1Temp[i]);
                if (i % elementSize == 15)
                    printf(" ");
            }
            printf("\n");
            for (int i = 0; i < m; i++) {
                printf("%s ", toBinary(e1Compressed[i], elementSize));
            }
            printf("\n");

        }
        if (!compare(n, e2Temp, m, e2Compressed)) {
            printf("ERRRRROR\n");
        }


        for (int i = 0; i < m; i++) {
            c1[i] = c2[i] = s[i] = 0;
        }

        calculateSparseAndUsual2NewVersion(n, hLength, m, h1Compact, e1Compressed, c1);
        calculateSparseAndUsual2NewVersion(n, hLength, m, h2Compact, e2Compressed, c2);


        /**проверка свертки генерации ключей*/
        for (int i = 0; i < n; i++) {
            c2Test[i] = c1Test[i] = 0;
        }

        my_computationF2(n, hLength, h1Compact, e1Temp, c1Test);
        my_computationF2(n, hLength, h2Compact, e2Temp, c2Test);

        if (!compare(n, c1Test, m, c1)) {
            printf("error in key get convolution c1\n");
        }
        if (!compare(n, c2Test, m, c2)) {
            printf("error in key get convolution c2\n");
        }


        for (int i = 0; i < m; i++) {
            s[i] = c1[i] ^ c2[i];
        }

        /**проверка вычисления синдрома*/
        for (int i = 0; i < n; i++) {
            sTest[i] = c1Test[i] ^ c2Test[i];
        }

        if (!compare(n, sTest, m, s)) {
            printf("error in syndrome S\n");
        }

        /**посчитали S*/



        for (int i = 0; i < m; i++)
            u[i] = v[i] = 0;

        for (int i = 0; i < n; i++)
            uTemp[i] = vTemp[i] = 0;

        for (int i = 0; i < m; i++) {
            sTemp[i] = s[i];
        }


        /**для теста*/
        for (int i = 0; i < n; i++) {
            sTestTemp[i] = sTest[i];
        }

        for (int z = 0; z < num_it; z++) {
            totalIt++;
            for (int i = 0; i < n; i++)
                upc1[i] = upc2[i] = 0;

            for (int i = 0; i < n; i++)
                sTempDecompressed[i] = 0;

            decompressArray(n, m, sTempDecompressed, sTemp);

            /**тест распаковки*/
            for (int i = 0; i < n; i++) {
                if (sTempDecompressed[i] != sTestTemp[i]) {
                    printf("error in s decompress %d\n", z);
                    break;
                }
            }

            computationZ(n, hLength, h1TransCompact, sTempDecompressed, upc1);
            computationZ(n, hLength, h2TransCompact, sTempDecompressed, upc2);

            /**тестирование на совпадение целочисленной свертки после распаковки*/
            for (int i = 0; i < n; i++) {
                upc1Test[i] = upc2Test[i] = 0;
            }
            computationZ(n, hLength, h1TransCompact, sTestTemp, upc1Test);
            computationZ(n, hLength, h2TransCompact, sTestTemp, upc2Test);
            for (int i = 0; i < n; i++) {
                if (upc1Test[i] != upc1[i]) {
                    printf("Error in convolution z with upc1\n");
                    break;
                }
            }
            for (int i = 0; i < n; i++) {
                if (upc2Test[i] != upc2[i]) {
                    printf("Error in convolution z with upc2\n");
                    break;
                }
            }


            for (int j = 0; j < n; j++) {
                if (upc1[j] >= T) {
                    u[j / elementSize] = u[j / elementSize] ^ mainMasks[j % elementSize];

                    uTemp[j] = uTemp[j] ^ 1;
                }

                if (upc2[j] >= T) {
                    /*printf("need to replace element %d\n",j);
                    printf("from %s", toBinary(v[j / elementSize],elementSize));*/
                    v[j / elementSize] = v[j / elementSize] ^ mainMasks[j % elementSize];
                    /*printf(" got %s\n", toBinary(v[j / elementSize],elementSize));
*/
                    vTemp[j] = vTemp[j] ^ 1;
                }

            }

            /**проверка изменения upc*/
            if (compare(n, vTemp, m, v) == 0 || compare(n, uTemp, m, u) == 0) {
                printf("asd %d it %d\n", asd, z);

                printf("T %d\n", T);

                printf("upc\n");
                for (int i = 0; i < n; i++) {
                    printf("%d ", upc2[i]);
                    if (i % elementSize == 15)
                        printf(" ");
                }
                printf("\n");

                printf("v and vTemp\n");

                for (int j = 0; j < n; j++) {
                    printf("%d", vTemp[j]);
                    if (j % elementSize == 15)
                        printf(" ");
                }
                printf("\n");
                for (int i = 0; i < m; i++) {
                    printf("%s ", toBinary(v[i], elementSize));
                }
                printf("\n");
                printf("\n");

                printf("u and uTemp\n");

                for (int j = 0; j < n; j++) {
                    printf("%d", uTemp[j]);
                    if (j % elementSize == 15)
                        printf(" ");
                }
                printf("\n");
                for (int i = 0; i < m; i++) {
                    printf("%s ", toBinary(u[i], elementSize));
                }
                printf("\n");
            }


            for (int i = 0; i < m; i++) {
                sTemp[i] = c1[i] = c2[i] = 0;
            }

            for (int i = 0; i < n; i++) {
                sTestTemp[i] = c1Test[i] = c2Test[i] = 0;
            }


            calculateSparseAndUsual2NewVersion(n, hLength, m, h1Compact, u, c1);
            calculateSparseAndUsual2NewVersion(n, hLength, m, h2Compact, v, c2);

            /**проверки свертки f2 внутри цикла*/

            my_computationF2(n, hLength, h1Compact, uTemp, c1Test);
            my_computationF2(n, hLength, h2Compact, vTemp, c2Test);

            if (!compare(n, c1Test, m, c1)) {
                printf("error is convolution F2 in cycle with c1\n");
            }
            if (!compare(n, c2Test, m, c2)) {
                printf("error is convolution F2 in cycle with c2\n");
            }


            for (int i = 0; i < m; i++) {
                sTemp[i] = s[i] ^ c1[i] ^ c2[i];
            }
            /**тестирование вычисления конечного s*/
            for (int i = 0; i < n; i++) {
                sTestTemp[i] = sTest[i] ^ c1Test[i] ^ c2Test[i];
            }

            if (!compare(n, sTestTemp, m, sTemp)) {
                printf("error in s temp generating in cycle\n");
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

            /**проверка сравнения s`*/
            flagTemp = 0;
            for (int i = 0; i < n; i++) {
                if (sTestTemp[i] != 0) {
                    flagTemp = 1;
                    break;
                }
            }
            if (flagTemp != flag) {
                printf("error in comparing sTemp to zero\n");
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

    free(h1Compact);
    free(h2Compact);
    free(h1TransCompact);
    free(h2TransCompact);
    free(e1Compact);
    free(e2Compact);
    free(e1Compressed);
    free(e2Compressed);
    free(c1);
    free(c2);
    free(s);
    free(u);
    free(v);
    free(sTemp);
    free(upc1);
    free(upc2);
    free(sTempDecompressed);

    printf("test decoder\n");
    printf("The elapsed time is %f seconds\n", time_spent);
    printf("One run is %f seconds\n", (time_spent / epoh));
    printf("%f \n", 100.0 * suc / (double) (epoh));
    printf("Total itarations is %f \n", (double) totalIt / epoh);

}

void decoderOptimized(long long seed, int epoh) {
    srand(seed);


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
    //unsigned short n = 50;

    unsigned short m = (n + elementSize - 1) / elementSize;

    /**компактное хранение h1*/
    unsigned short *h1Compact = calloc(hLength, sizeof(unsigned short));

    /**компактное хранение h2*/
    unsigned short *h2Compact = calloc(hLength, sizeof(unsigned short));

    /**компактное хранение h1 в обратном порядке*/
    unsigned short *h1TransCompact = calloc(hLength, sizeof(unsigned short));

    /**компактное хранение h2 в обратном порядке*/
    unsigned short *h2TransCompact = calloc(hLength, sizeof(unsigned short));



    /**компактное хранение e1*/
    unsigned short *e1Compact = calloc(eLength, sizeof(unsigned short));

    /**компактное хранение e2*/
    unsigned short *e2Compact = calloc(eLength, sizeof(unsigned short));

    unsigned short *e1Compressed = calloc(m, sizeof(unsigned short));
    unsigned short *e2Compressed = calloc(m, sizeof(unsigned short));

    unsigned short *c1 = calloc(m, sizeof(unsigned short));
    unsigned short *c2 = calloc(m, sizeof(unsigned short));
    unsigned short *s = calloc(m, sizeof(unsigned short));

    unsigned short *u = calloc(m, sizeof(unsigned short));
    unsigned short *v = calloc(m, sizeof(unsigned short));;

    unsigned short *sTemp = calloc(m, sizeof(unsigned short));
    unsigned short *upc1 = calloc(n, sizeof(unsigned short));
    unsigned short *upc2 = calloc(n, sizeof(unsigned short));
    unsigned short *sTempDecompressed = calloc(n, sizeof(unsigned short));;

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

        calculateSparseAndUsual2NewVersion(n, hLength, m, h1Compact, e1Compressed, c1);
        calculateSparseAndUsual2NewVersion(n, hLength, m, h2Compact, e2Compressed, c2);

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

            my_computationZ_last(n, hLength, h1TransCompact, sTempDecompressed, upc1);
            my_computationZ_last(n, hLength, h2TransCompact, sTempDecompressed, upc2);

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

            calculateSparseAndUsual2NewVersion(n, hLength, m, h1Compact, u, c1);
            calculateSparseAndUsual2NewVersion(n, hLength, m, h2Compact, v, c2);

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

    free(h1Compact);
    free(h2Compact);
    free(h1TransCompact);
    free(h2TransCompact);
    free(e1Compact);
    free(e2Compact);
    free(e1Compressed);
    free(e2Compressed);
    free(c1);
    free(c2);
    free(s);
    free(u);
    free(v);
    free(sTemp);
    free(upc1);
    free(upc2);
    free(sTempDecompressed);

    printf("Optimized decoder\n");
    printf("The elapsed time is %f seconds\n", time_spent);
    printf("One run is %f seconds\n", (time_spent / epoh));
    printf("%f \n", 100.0 * suc / (double) (epoh));
}

void decoder(long long seed, int epoh) {
    srand(seed);


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

    unsigned short *h1Compact = calloc(hLength, sizeof(unsigned short));
    unsigned short *h2Compact = calloc(hLength, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength, sizeof(unsigned short));
    unsigned short *h2TransCompact = calloc(hLength, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength, sizeof(unsigned short));
    unsigned short *e2Compact = calloc(eLength, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));
    unsigned short *e2Full = calloc(n, sizeof(unsigned short));


    unsigned short *u = calloc(n, sizeof(unsigned short));
    unsigned short *v = calloc(n, sizeof(unsigned short));
    unsigned short *c1 = calloc(n, sizeof(unsigned short));
    unsigned short *c2 = calloc(n, sizeof(unsigned short));
    unsigned short *s = calloc(n, sizeof(unsigned short));
    unsigned short *sTemp = calloc(n, sizeof(unsigned short));
    unsigned short *upc1 = calloc(n, sizeof(unsigned short));
    unsigned short *upc2 = calloc(n, sizeof(unsigned short));



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

        for (int i = 0; i < n; i++) {
            e1Full[i] = e2Full[i] = 0;
        }

        for (int i = 0; i < eLength; i++) {
            e1Full[e1Compact[i]] = 1;
            e2Full[e2Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = c2[i] = 0;
        }


        my_computationF2(n, hLength, h1Compact, e1Full, c1);
        my_computationF2(n, hLength, h2Compact, e2Full, c2);


        for (int i = 0; i < n; i++) {
            s[i] = c1[i] ^ c2[i];
        }

        /**посчитали S*/



        for (int i = 0; i < n; i++)
            u[i] = v[i] = 0;

        for (int i = 0; i < n; i++) {
            sTemp[i] = s[i];
        }

        for (int z = 0; z < num_it; z++) {
            for (int i = 0; i < n; i++) {
                upc1[i] = upc2[i] = 0;
            }


            my_computationZ_last(n, hLength, h1TransCompact, sTemp, upc1);
            my_computationZ_last(n, hLength, h2TransCompact, sTemp, upc2);

            for (int j = 0; j < n; j++) {
                if (upc1[j] >= T) {
                    u[j] = u[j] ^ 1;
                }

                if (upc2[j] >= T) {
                    v[j] = v[j] ^ 1;
                }

            }

            for (int i = 0; i < n; i++) {
                sTemp[i] = c1[i] = c2[i] = 0;
            }


            my_computationF2(n, hLength, h1Compact, u, c1);
            my_computationF2(n, hLength, h2Compact, v, c2);

            for (int i = 0; i < n; i++) {
                sTemp[i] = s[i] ^ c1[i] ^ c2[i];
            }

            /**проверка s` на ноль*/
            flag = 0;
            exitFlag = 0;
            for (int i = 0; i < n; i++) {
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


    free(h1Compact);
    free(h2Compact);
    free(h1TransCompact);
    free(h2TransCompact);
    free(e1Compact);
    free(e2Compact);

    free(e1Full);
    free(e2Full);

    free(u);
    free(v);
    free(c1);
    free(c2);
    free(s);
    free(sTemp);
    free(upc1);
    free(upc2);


    printf("original decoder\n");
    printf("The elapsed time is %f seconds\n", time_spent);
    printf("One run is %f seconds\n", (time_spent / epoh));
    printf("%f \n", 100.0 * suc / (double) epoh);
}

void fig1DecoderTestF2(long long seed, int epoh) {
    srand(seed);


    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;

    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));
    unsigned short *h1Full = calloc(n, sizeof(unsigned short));

    unsigned short *c1 = calloc(n, sizeof(unsigned short));


    int resForPrint = 0;
    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(e1Compact, n, eLength1);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = h1Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Full[e1Compact[i]] = 1;
        }

        for (int i = 0; i < hLength1; i++) {
            h1Full[h1Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = 0;
        }
        clock_t begin = clock();

        baseComputationF2(n, h1Full, e1Full, c1);
        clock_t end = clock();
        time_spent += (double) (end - begin) ;
        resForPrint += c1[0] ^ 1;
    }
    free(h1Compact);
    free(h1TransCompact);
    free(e1Compact);

    free(e1Full);
    free(h1Full);

    free(c1);


    printf("fig1DecoderTestF2\n");


    printf("The elapsed time for fig1DecoderTestF2 is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh);
    printf("The elapsed time for fig1DecoderTestF2 in clock for decoder is %f seconds\n", time_spent/epoh);


}

void fig1DecoderTestF2V2(long long seed, int epoh) {
    srand(seed);


    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;

    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));
    unsigned short *h1Full = calloc(n, sizeof(unsigned short));

    unsigned short *c1 = calloc(n, sizeof(unsigned short));


    int resForPrint = 0;
    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(e1Compact, n, eLength1);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = h1Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Full[e1Compact[i]] = 1;
        }

        for (int i = 0; i < hLength1; i++) {
            h1Full[h1Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = 0;
        }
        for (int z = 0; z < 10; z++) {
            clock_t begin = clock();

            baseComputationF2(n, h1Full, e1Full, c1);
            clock_t end = clock();
            time_spent += (double) (end - begin) ;
            resForPrint += c1[0] ^ 1;
        }
    }
    free(h1Compact);
    free(h1TransCompact);
    free(e1Compact);

    free(e1Full);
    free(h1Full);

    free(c1);


    printf("fig1DecoderTestF2V2\n");


    printf("The elapsed time for fig1DecoderTestF2V2 is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh);
    printf("The elapsed time for fig1DecoderTestF2V2 in clock for decoder is %f seconds\n", time_spent/epoh);

}

void fig1DecoderTestZ(long long seed, int epoh) {
    srand(seed);


    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;

    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));
    unsigned short *h1Full = calloc(n, sizeof(unsigned short));

    unsigned short *c1 = calloc(n, sizeof(unsigned short));


    int resForPrint = 0;
    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(e1Compact, n, eLength1);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = h1Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Full[e1Compact[i]] = 1;
        }

        for (int i = 0; i < hLength1; i++) {
            h1Full[h1Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = 0;
        }
        clock_t begin = clock();

        baseComputationZ(n, h1Full, e1Full, c1);
        clock_t end = clock();
        time_spent += (double) (end - begin) ;
        resForPrint += c1[0] ^ 1;
    }
    free(h1Compact);
    free(h1TransCompact);
    free(e1Compact);

    free(e1Full);
    free(h1Full);

    free(c1);


    printf("fig1DecoderTestZ\n");


    printf("The elapsed time for fig1DecoderTestZ is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh);
    printf("The elapsed time for fig1DecoderTestZ in clock for decoder is %f seconds\n", time_spent/epoh);

}

void fig2DecoderTestF2(long long seed, int epoh) {
    srand(seed);


    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;

    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));

    unsigned short *c1 = calloc(n, sizeof(unsigned short));


    int resForPrint = 0;
    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(e1Compact, n, eLength1);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Full[e1Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = 0;
        }
        clock_t begin = clock();

        my_computationF2_Fig2(n, hLength1, h1Compact, e1Full, c1);
        clock_t end = clock();
        time_spent += (double) (end - begin);
        resForPrint += c1[0] ^ 1;
    }
    free(h1Compact);
    free(h1TransCompact);
    free(e1Compact);

    free(e1Full);

    free(c1);


    printf("fig2DecoderTestF2\n");
    printf("The elapsed time for fig2DecoderTestF2 is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh);
    printf("The elapsed time for fig2DecoderTestF2 in clock for decoder is %f seconds\n", time_spent/epoh);

}

void fig2DecoderTestZ(long long seed, int epoh) {
    srand(seed);
    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;

    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));

    unsigned short *c1 = calloc(n, sizeof(unsigned short));


    int resForPrint = 0;
    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(e1Compact, n, eLength1);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Full[e1Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = 0;
        }
        clock_t begin = clock();

        my_computationZ_last(n, hLength1, h1Compact, e1Full, c1);
        clock_t end = clock();
        time_spent += (double) (end - begin);
        resForPrint += c1[0] ^ 1;
    }
    free(h1Compact);
    free(h1TransCompact);
    free(e1Compact);

    free(e1Full);

    free(c1);


    printf("fig2DecoderTestZ\n");
    printf("The elapsed time for fig2DecoderTestZ is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh);
    printf("The elapsed time for fig2DecoderTestZ in clock for decoder is %f seconds\n", time_spent/epoh);

}

void fig3DecoderTestF2(long long seed, int epoh) {
    srand(seed);


    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;

    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));

    unsigned short *c1 = calloc(n, sizeof(unsigned short));


    int resForPrint = 0;
    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(e1Compact, n, eLength1);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Full[e1Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = 0;
        }
        clock_t begin = clock();

        my_computationF3_Fig2(n, hLength1, h1Compact, e1Full, c1);
        clock_t end = clock();
        time_spent += (double) (end - begin) ;
        resForPrint += c1[0] ^ 1;
    }
    free(h1Compact);
    free(h1TransCompact);
    free(e1Compact);

    free(e1Full);

    free(c1);


    printf("fig3DecoderTestF2\n");
    printf("The elapsed time for fig3DecoderTestF2 is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh);
    printf("The elapsed time for fig3DecoderTestF2 in clock for decoder is %f seconds\n", time_spent/epoh);
}

void fig4DecoderTestF2(long long seed, int epoh) {
    srand(seed);
    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;
    unsigned short m = (n + elementSize - 1) / elementSize;


    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));

    unsigned short *c1 = calloc(n, sizeof(unsigned short));


    int resForPrint = 0;
    /**Конец объявления*/
    float suc = 0, fail = 0;

    double time_spent = 0;


    for (int asd = 0; asd < epoh; asd++) {

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(e1Compact, n, eLength1);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Full[e1Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            c1[i] = 0;
        }
        clock_t begin = clock();

        calculateSparseAndUsual2NewVersion(n, hLength1, m, h1Compact, e1Full, c1);
        clock_t end = clock();
        time_spent += (double) (end - begin) ;
        resForPrint += c1[0] ^ 1;
    }
    free(h1Compact);
    free(h1TransCompact);
    free(e1Compact);

    free(e1Full);

    free(c1);


    printf("fig4DecoderTestF2\n");
    printf("The elapsed time for fig4DecoderTestF2 is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh);
    printf("The elapsed time for fig4DecoderTestF2 in clock for decoder is %f seconds\n", time_spent/epoh);
}

void decoderOnlyBase(long long seed, int epoh) {
    srand(seed);


    /**|e1| + |e2| <= t
     * |e1| = |e2| - не обязательно
     *
     * Для примера
     * n = 4801
     * |h1| = |h2| = 45
     * |e1| + |e2| = 84
     * T = 45/2+4*/
    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 10;

    unsigned short n = 12323;

    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h2Compact = calloc(hLength2, sizeof(unsigned short));
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));
    unsigned short *h2TransCompact = calloc(hLength2, sizeof(unsigned short));
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));
    unsigned short *e2Compact = calloc(eLength2, sizeof(unsigned short));

    unsigned short *e1Full = calloc(n, sizeof(unsigned short));
    unsigned short *e2Full = calloc(n, sizeof(unsigned short));
    unsigned short *h1Full = calloc(n, sizeof(unsigned short));
    unsigned short *h2Full = calloc(n, sizeof(unsigned short));
    unsigned short *h1TransFull = calloc(n, sizeof(unsigned short));
    unsigned short *h2TransFull = calloc(n, sizeof(unsigned short));


    unsigned short *u = calloc(n, sizeof(unsigned short));
    unsigned short *v = calloc(n, sizeof(unsigned short));
    unsigned short *c1 = calloc(n, sizeof(unsigned short));
    unsigned short *c2 = calloc(n, sizeof(unsigned short));
    unsigned short *s = calloc(n, sizeof(unsigned short));
    unsigned short *sTemp = calloc(n, sizeof(unsigned short));
    unsigned short *upc1 = calloc(n, sizeof(unsigned short));
    unsigned short *upc2 = calloc(n, sizeof(unsigned short));



    /**Конец объявления*/
    float suc = 0, fail = 0;

    long long time_spent = 0;
    long long time_spentFull = 0;

    int realIt = 0;


    for (int asd = 0; asd < epoh; asd++) {
        clock_t beginFull = clock();

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }
        for (int i = 0; i < hLength2; i++) {
            h2Compact[i] = h2TransCompact[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e1Compact[i]  = 0;
        }
        for (int i = 0; i < eLength2; i++) {
            e2Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(h2Compact, n, hLength2);
        generateSparseArray(e1Compact, n, eLength1);
        generateSparseArray(e2Compact, n, eLength2);

        for (int i = 0; i < n; i++) {
            c1[i] = c2[i] = 0;
            h1Full[i] = h2Full[i] = h1TransFull[i] = h2TransFull[i] = 0;
        }

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
            h1Full[h1Compact[i]] = 1;
            h1TransFull[n - h1Compact[i]] = 1;
        }

        for (int i = 0; i < hLength2; i++) {
            h2TransCompact[i] = n - h2Compact[i];
            h2Full[h2Compact[i]] = 1;
            h2TransFull[n - h2Compact[i]] = 1;
        }

        for (int i = 0; i < n; i++) {
            e1Full[i] = e2Full[i] = 0;
        }

        for (int i = 0; i < eLength1; i++) {
            e2Full[e2Compact[i]] = 1;
        }

        for (int i = 0; i < eLength2; i++) {
            e2Full[e2Compact[i]] = 1;
        }


        baseComputationF2(n, h1Full, e1Full, c1);
        baseComputationF2(n, h2Full, e2Full, c2);

        for (int i = 0; i < n; i++) {
            s[i] = c1[i] ^ c2[i];
        }

        /**посчитали S*/

        clock_t begin = clock();


        for (int i = 0; i < n; i++)
            u[i] = v[i] = 0;

        for (int i = 0; i < n; i++) {
            sTemp[i] = s[i];
        }

        for (int z = 0; z < num_it; z++) {
            realIt++;

            for (int i = 0; i < n; i++) {
                upc1[i] = upc2[i] = 0;
            }


            baseComputationZ(n, h1TransFull, sTemp, upc1);
            baseComputationZ(n, h2TransFull, sTemp, upc2);

            for (int j = 0; j < n; j++) {
                if (upc1[j] >= T) {
                    u[j] = u[j] ^ 1;
                }

                if (upc2[j] >= T) {
                    v[j] = v[j] ^ 1;
                }

            }

            for (int i = 0; i < n; i++) {
                sTemp[i] = c1[i] = c2[i] = 0;
            }


            baseComputationF2(n, h1Full, u, c1);
            baseComputationF2(n, h2Full, v, c2);

            for (int i = 0; i < n; i++) {
                sTemp[i] = s[i] ^ c1[i] ^ c2[i];
            }


            /**проверка s` на ноль*/
            flag = 0;
            exitFlag = 0;
            for (int i = 0; i < n; i++) {
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
        time_spent += (end - begin) ;
        time_spentFull += (end - beginFull);
    }


    free(h1Compact);
    free(h2Compact);
    free(h1TransCompact);
    free(h2TransCompact);
    free(e1Compact);
    free(e2Compact);

    free(e1Full);
    free(e2Full);

    free(h1Full);
    free(h2Full);
    free(h1TransFull);
    free(h2TransFull);

    free(u);
    free(v);
    free(c1);
    free(c2);
    free(s);
    free(sTemp);
    free(upc1);
    free(upc2);


    printf("decoderOnlyBase\n");
    printf("The elapsed time for decoder is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh/num_it);
    printf("The elapsed time for full is %f seconds\n", ((double)time_spentFull)/ CLOCKS_PER_SEC/epoh/num_it);
    printf("The elapsed time in clock for decoder is %lld seconds\n", time_spent/epoh/num_it);
    printf("The elapsed time in clock for full is %lld seconds\n", time_spentFull/epoh/num_it);
    printf("One run is %f seconds\n", (((double)time_spent) / epoh/num_it));
    printf("%f \n", 100.0 * suc / (double) epoh);
    printf("real iterations %f \n", realIt / (double) epoh);
}


void decoderOnlyOptimized(long long seed, int epoh) {
    srand(seed);


    /**|e1| + |e2| <= t
     * |e1| = |e2| - не обязательно
     *
     * Для примера
     * n = 4801
     * |h1| = |h2| = 45
     * |e1| + |e2| = 84
     * T = 45/2+4*/
    unsigned short hLength1 = 71;
    unsigned short hLength2 = 71;
    unsigned short eLength1 = 67;
    unsigned short eLength2 = 67;
    unsigned T = /*hLength1 / 2 + 5;*/61;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 12323;

    unsigned short m = (n + elementSize - 1) / elementSize;

    /**компактное хранение h1*/
    unsigned short *h1Compact = calloc(hLength1, sizeof(unsigned short));

    /**компактное хранение h2*/
    unsigned short *h2Compact = calloc(hLength2, sizeof(unsigned short));

    /**компактное хранение h1 в обратном порядке*/
    unsigned short *h1TransCompact = calloc(hLength1, sizeof(unsigned short));

    /**компактное хранение h2 в обратном порядке*/
    unsigned short *h2TransCompact = calloc(hLength2, sizeof(unsigned short));



    /**компактное хранение e1*/
    unsigned short *e1Compact = calloc(eLength1, sizeof(unsigned short));

    /**компактное хранение e2*/
    unsigned short *e2Compact = calloc(eLength2, sizeof(unsigned short));

    unsigned short *e1Compressed = calloc(m, sizeof(unsigned short));
    unsigned short *e2Compressed = calloc(m, sizeof(unsigned short));

    unsigned short *c1 = calloc(m, sizeof(unsigned short));
    unsigned short *c2 = calloc(m, sizeof(unsigned short));
    unsigned short *s = calloc(m, sizeof(unsigned short));

    unsigned short *u = calloc(m, sizeof(unsigned short));
    unsigned short *v = calloc(m, sizeof(unsigned short));;

    unsigned short *sTemp = calloc(m, sizeof(unsigned short));
    unsigned short *upc1 = calloc(n, sizeof(unsigned short));
    unsigned short *upc2 = calloc(n, sizeof(unsigned short));
    unsigned short *sTempDecompressed = calloc(n, sizeof(unsigned short));;

    /**Конец объявления*/
    float suc = 0, fail = 0;
    int totalIt = 0;

    long long time_spent = 0;
    long long time_spentFull = 0;


    for (int asd = 0; asd < epoh; asd++) {
        clock_t beginFull = clock();

        for (int i = 0; i < hLength1; i++) {
            h1Compact[i] = h1TransCompact[i] = 0;
        }
        for (int i = 0; i < hLength2; i++) {
            h2Compact[i] = h2TransCompact[i] = 0;
        }



        for (int i = 0; i < eLength1; i++) {
            e1Compact[i] = 0;
        }
        for (int i = 0; i < eLength2; i++) {
            e2Compact[i] = 0;
        }

        /**заполняем вектора данными*/
        generateSparseArray(h1Compact, n, hLength1);
        generateSparseArray(h2Compact, n, hLength2);
        generateSparseArray(e1Compact, n, eLength1);
        generateSparseArray(e2Compact, n, eLength2);

        for (int i = 0; i < hLength1; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }
        for (int i = 0; i < hLength2; i++) {
            h1TransCompact[i] = n - h1Compact[i];
        }

        for (int i = 0; i < m; i++) {
            e1Compressed[i] = e2Compressed[i] = 0;
        }

        createDenseArrayFromCompact(e1Compressed, m, e1Compact, eLength1);
        createDenseArrayFromCompact(e2Compressed, m, e2Compact, eLength2);


        for (int i = 0; i < m; i++) {
            c1[i] = c2[i] = s[i] = 0;
        }

        calculateSparseAndUsual2NewVersion(n, hLength1, m, h1Compact, e1Compressed, c1);
        calculateSparseAndUsual2NewVersion(n, hLength2, m, h2Compact, e2Compressed, c2);

        for (int i = 0; i < m; i++) {
            s[i] = c1[i] ^ c2[i];
        }

        /**посчитали S*/

        clock_t begin = clock();

        for (int i = 0; i < m; i++)
            u[i] = v[i] = 0;

        for (int i = 0; i < m; i++) {
            sTemp[i] = s[i];
        }

        for (int z = 0; z < num_it; z++) {
            totalIt++;
            for (int i = 0; i < n; i++)
                upc1[i] = upc2[i] = 0;

            for (int i = 0; i < n; i++)
                sTempDecompressed[i] = 0;

            decompressArray(n, m, sTempDecompressed, sTemp);

            my_computationZ_last(n, hLength1, h1TransCompact, sTempDecompressed, upc1);
            my_computationZ_last(n, hLength2, h2TransCompact, sTempDecompressed, upc2);

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

            calculateSparseAndUsual2NewVersion(n, hLength1, m, h1Compact, u, c1);
            calculateSparseAndUsual2NewVersion(n, hLength2, m, h2Compact, v, c2);

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
        time_spent +=  (end - begin) ;
        time_spentFull += (end - beginFull) ;


    }

    free(h1Compact);
    free(h2Compact);
    free(h1TransCompact);
    free(h2TransCompact);
    free(e1Compact);
    free(e2Compact);
    free(e1Compressed);
    free(e2Compressed);
    free(c1);
    free(c2);
    free(s);
    free(u);
    free(v);
    free(sTemp);
    free(upc1);
    free(upc2);
    free(sTempDecompressed);

    printf("decoderOnlyOptimized\n");
    printf("The elapsed time for decoder is %f seconds\n", ((double)time_spent)/ CLOCKS_PER_SEC/epoh/num_it);
    printf("The elapsed time for full is %f seconds\n", ((double)time_spentFull)/ CLOCKS_PER_SEC/epoh/num_it);
    printf("The elapsed time in clock for decoder is %lld clock\n", time_spent/epoh/num_it);
    printf("The elapsed time in clock for full is %lld clock\n", time_spentFull/epoh/num_it);
    printf("One run is %f seconds\n", (((double)time_spent/CLOCKS_PER_SEC) / epoh/num_it));
    printf("%f \n", 100.0 * suc / (double) epoh);
    printf("real iterations %f \n", totalIt / (double) epoh);
}


int test(int nFrom, int nTo) {
    int arrSize = 1;
    int counter = 0;
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

            //res2[0]=res2[0]^1;
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
            } else {
                counter++;
            }

            free(b);
            free(b2);
            free(res);
            free(res2);
        }

    }

    free(inda);

    return counter;
}

int main() {
    int epoh = 1;

    long long seed = time(0);

    /*начало блока инициализации*/
    mainMasks = calloc(elementSize, sizeof(unsigned short));
    mainMasks[0] = 1 << (elementSize - 1);
    for (int i = 1; i < elementSize; i++)
        mainMasks[i] = mainMasks[i - 1] >> 1;

    denseArrayFromCompactMasks = calloc(elementSize, sizeof(unsigned short));
    denseArrayFromCompactMasks[elementSize - 1] = 1;
    for (int i = elementSize - 2; i >= 0; i--) {
        denseArrayFromCompactMasks[i] = denseArrayFromCompactMasks[i + 1] << 1;
    }

    calculateSparseAndUsual2Masks = calloc(elementSize, sizeof(unsigned short));
    calculateSparseAndUsual2NegMasks = calloc(elementSize, sizeof(unsigned short));
    unsigned short temp = -1;
    for (int i = 0; i < elementSize; i++) {
        calculateSparseAndUsual2Masks[i] = temp;
        calculateSparseAndUsual2NegMasks[i] = ~temp;
        temp -= 1 << (i);
    }

    /*конец блока инициализации*/
    //fig1DecoderTestF2(seed, epoh);
    /*fig1DecoderTestZ(seed, epoh);
    fig2DecoderTestF2(seed, epoh);
    fig2DecoderTestZ(seed, epoh);
    fig3DecoderTestF2(seed, epoh);
    fig4DecoderTestF2(seed, epoh);*/

    decoderOnlyBase(seed, epoh);

    decoderOnlyOptimized(seed, epoh);
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

void
my_computationF2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                 unsigned short *c) {
    //printf("begin computation\n");
    int j;
    int new_j;
    for (unsigned short *i = inda; i != inda + arrSize; i++) {
        j = 0;
        new_j = n - *i + j;
        for (; j < *i; j++) {
            c[j] = (c[j] + b[new_j]) % 2;
            new_j++;
        }
        new_j = j - *i;
        for (; j < n; j++) {
            c[j] = (c[j] + b[new_j]) % 2;
            new_j++;
        }
    }
}

void
my_computationF2_Fig2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                      unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - inda[i] + n) % n]) % 2;
        }
    }
}

void
my_computationF3_Fig2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                      unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] ^ b[(j - inda[i] + n) % n]);
        }
    }
}

void baseComputationF2(unsigned short n, const unsigned short *a, const unsigned short *b, unsigned short *c) {
    int cSum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - i + n) % n] * a[i]) % 2;
        }
    }
}

void baseComputationZ(unsigned short n, const unsigned short *a, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - i + n) % n] * a[i]);
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

void
my_computationZ_last(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                     unsigned short *c) {
    //printf("begin computation\n");
    unsigned short *pointer_c;
    int j;
    int new_j;
    unsigned short *end = inda + arrSize;
    for (unsigned short *i = inda; i != end; i++) {
        pointer_c = c;
        j = 0;
        new_j = n - *i;
        for (; j < *i; j++) {
            *pointer_c = (*pointer_c + b[new_j]);
            new_j++;
            pointer_c++;
        }
        new_j = j - *i;
        for (; j < n; j++) {
            *pointer_c = (*pointer_c + b[new_j]);
            new_j++;
            pointer_c++;

        }
    }

}

void
my_computationZ_Fig2(unsigned short n, int arrSize, const unsigned short *inda, const unsigned short *b,
                     unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < n; j++) {
            c[j] = (c[j] + b[(j - inda[i] + n) % n]);
        }
    }
}




