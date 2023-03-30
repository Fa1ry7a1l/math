#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>


#define elementSize 16

/**стандартная свертка в F2*/
void
computationF2(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);

/**стандартная свертка в Z*/
void
computationZ(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);


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

    /*for (int i = 0; i < elementSize; i++) {
        printf("%d ", i);
        for (int j = 0; j < m; j++) {
            printf("%s ", toBinary(arr[i][j], elementSize));
        }
        printf("\n");
    }*/

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


void calculateSparsedAndUsual1(unsigned short M, int arrSize, int m, const unsigned short *inda, unsigned short **b2,
                               unsigned short *res) {
    /**подсчет вектора новым методом*/
    /**остаток деления длины вектора на размер ячейки*/
    unsigned short mod = M % elementSize;

    /**метод для M не кратной elemntSize*/
    if (mod != 0) {
        /**маска для подсчета последнего элемента*/
        int negMod = elementSize - mod;
        unsigned short lastElementMask = -1;
        unsigned short finishElementMask = -1;
        lastElementMask = lastElementMask << (elementSize - mod);
        finishElementMask = finishElementMask >> (elementSize - mod);
        unsigned short startElementMask = ~finishElementMask;
        if (inda[arrSize - 1] <= mod) {
            /**ветка содержит 2 почти одинаковых цикла, с отличиями в помеченном блоке для избежания if внутри цикла*/
            int i = 0;
            for (; inda[i] <= mod; i++) {

                int rowNumber;
                int startColumn = 0;
                int el = 0;
                /**начало отличия*/
                rowNumber = inda[i];
                /**конец отличия*/


                /**обрабатываем все элементы со сдвигом так, чтобы */
                for (int j = startColumn; j < m - 1; j++) {
                    res[el] = res[el] ^ b2[rowNumber][j];
                    el++;
                }


                /**обработка последнего элемента*/

                unsigned short nextVal;
                unsigned short startNumber;
                unsigned short finishNumber;

                /**подсчет переходного элемента*/
                startNumber = b2[rowNumber][m - 1] & lastElementMask;
                finishNumber = (b2[rowNumber][0] & startElementMask) >> mod;
                nextVal = startNumber + finishNumber;
                res[el] = res[el] ^ nextVal;
                el++;

                for (int j = 0; j < startColumn - 1; j++) {
                    res[el] = res[el] ^ (((b2[rowNumber][j] & finishElementMask) << negMod) +
                                         ((b2[rowNumber][j + 1] & startElementMask) >> mod));
                    el++;
                }
                /**отдельная обработка последнего элемента массива должна предотвратить выход за пределы*/
                res[el] = res[el] ^ ((b2[rowNumber][startColumn - 1] & finishElementMask) << negMod);
            }

            //начало второго цикла
            for (; i < arrSize; i++) {

                int rowNumber;
                int startColumn = 0;
                int el = 0;
                /**начало отличия*/
                int j = M - inda[i];
                rowNumber = (elementSize - (j % elementSize)) % elementSize;
                startColumn = (j + elementSize - 1) / elementSize;
                /**конец отличия*/


                /**обрабатываем все элементы со сдвигом так, чтобы */
                for (int j = startColumn; j < m - 1; j++) {
                    res[el] = res[el] ^ b2[rowNumber][j];
                    el++;
                }


                /**обработка последнего элемента*/

                unsigned short nextVal;
                unsigned short startNumber;
                unsigned short finishNumber;

                /**подсчет переходного элемента*/
                startNumber = b2[rowNumber][m - 1] & lastElementMask;
                finishNumber = (b2[rowNumber][0] & startElementMask) >> mod;
                nextVal = startNumber + finishNumber;
                res[el] = res[el] ^ nextVal;
                el++;

                for (int j = 0; j < startColumn - 1; j++) {
                    res[el] = res[el] ^ (((b2[rowNumber][j] & finishElementMask) << negMod) +
                                         ((b2[rowNumber][j + 1] & startElementMask) >> mod));
                    el++;
                }
                /**отдельная обработка последнего элемента массива должна предотвратить выход за пределы*/
                res[el] = res[el] ^ ((b2[rowNumber][startColumn - 1] & finishElementMask) << negMod);


            }
        } else {
            for (int i = 0; i < arrSize; i++) {

                int rowNumber;
                int startColumn = 0;
                int el = 0;
                rowNumber = inda[i];


                /**обрабатываем все элементы со сдвигом так, чтобы */
                for (int j = startColumn; j < m - 1; j++) {
                    res[el] = res[el] ^ b2[rowNumber][j];
                    el++;
                }


                /**обработка последнего элемента*/

                unsigned short nextVal;
                unsigned short startNumber;
                unsigned short finishNumber;

                /**подсчет переходного элемента*/
                startNumber = b2[rowNumber][m - 1] & lastElementMask;
                finishNumber = (b2[rowNumber][0] & startElementMask) >> mod;
                nextVal = startNumber + finishNumber;
                res[el] = res[el] ^ nextVal;
                el++;

                for (int j = 0; j < startColumn - 1; j++) {
                    res[el] = res[el] ^ (((b2[rowNumber][j] & finishElementMask) << negMod) +
                                         ((b2[rowNumber][j + 1] & startElementMask) >> mod));
                    el++;
                }
                /**отдельная обработка последнего элемента массива должна предотвратить выход за пределы*/
                res[el] = res[el] ^ ((b2[rowNumber][startColumn - 1] & finishElementMask) << negMod);


            }
        }
        res[m - 1] = res[m - 1] & lastElementMask;
    } else {
        for (int i = 0; i < arrSize; i++) {

            int rowNumber = inda[i] % elementSize;
            int startColumn = ((M - inda[i] + elementSize - 1) / elementSize) % m;
            int el = 0;

            /**обрабатываем все элементы со сдвигом так, до конца*/
            for (int j = startColumn; j < m; j++) {
                res[el] = res[el] ^ b2[rowNumber][j];
                el++;
            }
            for (int j = 0; j < startColumn; j++) {
                res[el] = res[el] ^ b2[rowNumber][j];
                el++;
            }


        }
    }
}

//todo не обрабатывает случай Gmod == 0;
void calculateSparsedAndUsual2(unsigned short M, int arrSize, int m, const unsigned short *inda, unsigned short *b2,
                               unsigned short *res) {
    unsigned short masks[elementSize];
    unsigned short negMasks[elementSize];
    unsigned short temp = 0;
    for (int i = 0; i < elementSize; i++) {
        masks[i] = temp;
        negMasks[i] = ~temp;
        temp += 1 << (elementSize - i - 1);
        printf("%s %s\n", toBinary(masks[i], elementSize), toBinary(negMasks[i], elementSize));
    }
    unsigned short Gmod = M % elementSize;
    unsigned short GmodNeg = elementSize - Gmod;
    printf("Gmod %d\n", Gmod);

    for (int j = 0; j < m; j++) {
        printf("%s ", toBinary(b2[j], elementSize));
    }
    printf("\n");


    for (int i = 0; i < arrSize; i++) {
        unsigned short pos = (M - inda[i]) % M;
        printf("pos %d\n", pos);
        unsigned short mod = pos % elementSize;
        printf("mod %d\n", mod);
        unsigned short negMod = elementSize - mod;
        printf("negMod %d\n", negMod);

        unsigned short startColumn = pos / elementSize;
        unsigned short it = 0;

        if (mod != 0) {
            for (int j = startColumn; j < m - 2; j++) {
                printf("%s\n", toBinary(((b2[j] & negMasks[mod]) << mod), elementSize));
                printf("%s\n", toBinary(((b2[j + 1] & masks[mod]) >> negMod), elementSize));
                res[it++] += ((b2[j] & negMasks[mod]) << mod) + ((b2[j + 1] & masks[negMod]) >> negMod);
                for (int z = 0; z < m; z++)
                    printf("%s ", toBinary(res[z], elementSize));
                printf("\n");
            }
            /*for(int z = 0;z<m;z++)
                printf("%s ", toBinary(b2[z],elementSize));
            printf("\n");*/
            unsigned short Nmod = abs(Gmod - mod);
            unsigned short NmodNeg = elementSize - Nmod;
            printf("Nmod %d\n", Nmod);
            printf("NmodNeg %d\n", NmodNeg);

            printf("-----------------\n");

            for (int z = 0; z < m; z++)
                printf("%s ", toBinary(res[z], elementSize));
            printf("\n");
            res[it] += ((b2[m - 2] & negMasks[mod]) << mod) + ((b2[m - 1] & masks[mod]) >> negMod);
            for (int z = 0; z < m; z++)
                printf("%s ", toBinary(res[z], elementSize));
            printf("\n");

            if (Gmod > mod) {
                it++;

                //первая часть в конце
                res[it] += ((b2[m - 1] & negMasks[mod]) << mod);

                for (int z = 0; z < m; z++)
                    printf("%s ", toBinary(res[z], elementSize));
                printf("\n");
                //добиваем до конца
                res[it] += ((b2[0] & masks[NmodNeg]) >> Nmod);
                it++;
                for (int z = 0; z < m; z++)
                    printf("%s ", toBinary(res[z], elementSize));
                printf("\n");

                for (int j = 0; it < m; it++, j++) {
                    printf("%s\n", toBinary(((b2[j] & negMasks[NmodNeg]) << NmodNeg), elementSize));
                    printf("%s\n", toBinary(((b2[j + 1] & masks[NmodNeg]) >> Nmod), elementSize));
                    res[it] += ((b2[j] & negMasks[NmodNeg]) << NmodNeg) + ((b2[j + 1] & masks[NmodNeg]) >> Nmod);
                }


            } else {


                res[it] += b2[0] & masks[Nmod] >> NmodNeg;
                it++;
                for (int z = 0; z < m; z++)
                    printf("%s ", toBinary(res[z], elementSize));
                printf("\n");

                for (int j = 0; it < m; it++, j++) {
                    printf("%s\n", toBinary(((b2[j] & negMasks[Nmod]) << Nmod), elementSize));
                    printf("%s\n", toBinary(((b2[j + 1] & masks[Nmod]) >> NmodNeg), elementSize));
                    res[it] += ((b2[j] & negMasks[Nmod]) << Nmod) + ((b2[j + 1] & masks[Nmod]) >> NmodNeg);
                }
            }

            printf("-----------------\n");


            for (int z = 0; z < m; z++)
                printf("%s ", toBinary(res[z], elementSize));
            printf("\n");
        } else {
            for (int j = startColumn; j < m - 1; j++) {
                res[it++] += b2[j];
            }


        }

    }
    printf("final\n");
    for (int z = 0; z < m; z++)
        printf("%s ", toBinary(res[z], elementSize));
    printf("\n");
}

/**разжимает вектор*/
void decompressArray(unsigned short M, unsigned short m, unsigned short *arrTo, unsigned short *arrFrom) {
    int it = M - 1;
    unsigned short temp = arrFrom[m - 1];
    unsigned short mod = M % elementSize;
    temp = temp >> (elementSize - (M % elementSize)) % elementSize;
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

int main() {

    srand(time(0));

    /**M - количество элементов в векторе A, B , C*//*
    unsigned short M = 36;
    //unsigned short M = 45371;
    *//**arrSize - количество элементов в сжатом векторе A, где хранятся только номера ненулевых эементов разреженного вектора A*//*
    int arrSize = 1;
    //int arrSize = 223;
    *//**m - размер сжатого вектора B, C*//*
    int m = (M + elementSize - 1) / elementSize;
    *//**inda - сжатый вектор A, A - разреженный*//*
    unsigned short inda[arrSize];
    inda[0] = 4;
    *//*for (int z = 0; z < M; z++) {
        printf("\n\n");
        inda[0] = z;*//*
    *//**b - Вектор с которым происходит свертка*//*
    unsigned short b[M];

    *//**результирующий вевктор стандартного примера (размера M)*//*
    unsigned short c[M];
    for (int i = 0; i < M; i++) {
        c[i] = 0;
    }

    *//**результирующий сжатый вектор второго способа (размера m)*//*
    unsigned short *c2 = calloc(m, sizeof(unsigned short));
    for (int i = 0; i < m; i++) {
        c2[i] = 0;
    }





    *//**сжатая версия b*//*
    unsigned short **b2 = calloc(elementSize, sizeof(unsigned short *));

    for (int i = 0; i < elementSize; i++)
        b2[i] = calloc(m, sizeof(unsigned short));

    // generateSparseArray(inda, M, arrSize);

    generateDenseArray(b, M);
    *//*for (int i = 0; i < M; i++)
        b[i] = 1;*//*
    createDenseArray(b2, m, b, M);

    calculateSparsedAndUsual2(M, arrSize, m, inda, b2[0], c2);


    unsigned short *bCopy = calloc(M, sizeof(unsigned short *));*/
    /*decompressArray(M, m, bCopy, b2[0]);
    printf("b and bCopy\n");
    for (int i = 0; i < M; i++)
        printf("%d", b[i]);
    printf("\n");
    for (int i = 0; i < M; i++)
        printf("%d", bCopy[i]);*/




    /*printf("inda\n");
    for (int i = 0; i < arrSize; i++) {
        printf("%d ", inda[i]);
    }
    printf("\n");

    printf("b\n");
    for (int i = 0; i < M; i++) {
        printf("%d", b[i]);
    }

    printf("\n");*/

    /**обычный метод*/
    /*computationF2(M, arrSize, inda, b, c);
    for (int i = 0; i < M; i++) {
        printf("%d", c[i]);
    }
    printf("\n");*/

    /**Печать обычного метода подсчета вектора*/
    /*printf("res\n");
    for (int i = 0; i < M; i++) {
        if (i != 0 && (i % elementSize) == 0)
            printf(" ");
        printf("%d", c[i]);
    }
    for (int i = 0; i < (elementSize - M % elementSize) % elementSize; i++) {
        printf("0");
    }
    printf("\n");*/

    /*double time_spent_generating = 0;
    double time_spent_calculating = 0;
    double time_spent_generating_compressed = 0;
    int answ = 0;
    for (int k = 0; k < 10000; k++) {
        clock_t begin_generating = clock();

        generateSparseArray(inda, M, arrSize);
        generateDenseArray(b, M);


        clock_t end_generating = clock();
        time_spent_generating += (double) (end_generating - begin_generating) / CLOCKS_PER_SEC;

        clock_t begin_generating_compressed = clock();
        createDenseArray(b2, m, b, M);
        clock_t finish_generating_compressed = clock();
        time_spent_generating_compressed +=
                (double) (finish_generating_compressed - begin_generating_compressed) / CLOCKS_PER_SEC;


        clock_t begin_calculating = clock();
        calculateSparsedAndUsual1(M, arrSize, m, inda, b2, c2);
        int mod = M % elementSize;


        unsigned short lastElementMask = -1 << (elementSize - mod);
        clock_t end_calculating = clock();
        time_spent_calculating += (double) (end_calculating - begin_calculating) / CLOCKS_PER_SEC;
        answ += ((c2[m - 1] & lastElementMask) >> (elementSize - mod)) % 2;
    }
    printf("My new Method\n");
    printf("The elapsed time generating is %f seconds\n", time_spent_generating);
    printf("One run is %f seconds\n", (time_spent_generating / 10000));
    printf("The elapsed time generating compressed vector is %f seconds\n", time_spent_generating_compressed);
    printf("One run is %f seconds\n", (time_spent_generating_compressed / 10000));
    printf("The elapsed time calculating is %f seconds\n", time_spent_calculating);
    printf("One run is %f seconds\n", (time_spent_calculating / 10000));
    printf("answ is %d\n", answ);*/

    //printf("res\n");
    /*for (int i = 0; i < m; i++) {
        printf("%s ", toBinary(c2[i], elementSize));
    }
    printf("\n");*/
    //}

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
    unsigned T = hLength / 2 + 4;

    unsigned short flag = 0;
    unsigned short exitFlag = 0;

    unsigned short num_it = 100;

    unsigned short n = 4801;

    unsigned short h1Compact[hLength];
    unsigned short h2Compact[hLength];
    unsigned short h1TransCompact[hLength];
    unsigned short h2TransCompact[hLength];
    unsigned short e1Compact[eLength];
    unsigned short e2Compact[eLength];
    generateSparseArray(h1Compact, n, hLength);
    generateSparseArray(h2Compact, n, hLength);
    generateSparseArray(e1Compact, n, eLength);
    generateSparseArray(e2Compact, n, eLength);

    for (int i = 0; i < hLength; i++) {
        h1TransCompact[hLength - 1 - i] = n - 1 - h1Compact[i];
        h2TransCompact[hLength - 1 - i] = n - 1 - h2Compact[i];
    }

    unsigned short e1[n];
    unsigned short e2[n];
    for (int i = 0; i < n; i++)
        e1[i] = e2[i] = 0;
    for (int i = 0; i < eLength; i++)
        e1[e1Compact[i]] = 1;
    for (int i = 0; i < eLength; i++)
        e2[e2Compact[i]] = 1;

    unsigned short c1[n];
    unsigned short c2[n];
    unsigned short s[n];

    computationF2(n, hLength, h1Compact, e1, c1);
    computationF2(n, hLength, h2Compact, e2, c2);
    for (int i = 0; i < n; i++) {
        s[i] = c1[i] ^ c2[i];
    }
    for (int i = 0; i < 150; i++)
        printf("%d", s[i]);
    printf("\n");

    unsigned short u[n], v[n];
    for (int i = 0; i < n; i++)
        u[i] = v[i] = 0;

    unsigned short sTemp[n];
    for (int i = 0; i < n; i++) {
        sTemp[i] = s[i];
    }

    for (int z = 0; z < num_it; z++) {
        unsigned short upc1[n], upc2[n];
        for (int i = 0; i < n; i++)
            upc1[i] = upc2[i] = 0;
        computationZ(n,hLength,h1TransCompact,sTemp,upc1);
        computationZ(n,hLength,h2TransCompact,sTemp,upc2);
        for(int j = 0;j<n;j++)
        {
            if(upc1[j] >=T)
                u[j] = u[j] ^ 1;
            if(upc2[j] >=T)
                v[j] = v[j] ^ 1;
        }

        for(int i = 0;i<n;i++)
        {
            sTemp[i]=c1[i] = c2[i] = 0;
        }
        computationF2(n,hLength,h1Compact,u,c1);
        computationF2(n,hLength,h2Compact,v,c2);
        for(int i = 0;i<n;i++)
        {
            sTemp[i] = s[i] ^ c1[i] ^ c2[i];
        }

        flag = 0;
        for(int i = 0;i<n;i++)
        {
            if(sTemp[i] != 0)
            {
                flag = 1;
                break;
            }
        }
        if(!flag)
        {
            exitFlag = 1;
            break;
        }

    }

    if(exitFlag)
    {
        printf("e1\n");
        for(int i = 0;i<n;i++)
        {
            printf("%d",e1[i]);
        }
        printf("\n");
        printf("e2\n");
        for(int i = 0;i<n;i++)
        {
            printf("%d",e2[i]);
        }
        printf("\n");
        printf("u\n");
        for(int i = 0;i<n;i++)
        {
            printf("%d",u[i]);
        }
        printf("\n");
        printf("v\n");
        for(int i = 0;i<n;i++)
        {
            printf("%d",v[i]);
        }
        printf("\n");
    }
    else
    {
        printf("fail\n");
    }


    return 0;
}

void
computationF2(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < M; j++) {
            c[j] = (c[j] + b[(inda[i] - j + M) % M]) % 2;
        }
    }
}

void
computationZ(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < M; j++) {
            c[j] = (c[j] + b[(inda[i] - j + M) % M]);
        }
    }
}

