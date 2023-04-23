#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>


#define elementSize 16


/**стандартная свертка разреженного и не разреженного векторов F2
 * M - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (F2)*/
void
computationF2(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c);

/**стандартная свертка разреженного и не разреженного векторов в Z
 * M - длина неразреженного вектора
 * arrSize - количество ненулевых элементов резреженного вектора
 * inda - вектор номеров позиций ненулевых координт разреженного вектора
 * b - неразреженный вектор F2
 * c - результирующий вектор (Z)*/
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

/**генерирует разреженный массив, как массив индексов ненулевых элементов, упорядоченных в порядке возрастания.
 * arr - результирующий вектор индексов
 * m - максимальный индекс в векторе
 * res - количество элементов в результирующем векторе
 * */
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

void calculateSparseAndUsual2(unsigned short M, int arrSize, int m, const unsigned short *inda, unsigned short *b2,
                              unsigned short *res) {
    unsigned short masks[elementSize];
    unsigned short negMasks[elementSize];
    unsigned short modG = M % elementSize;
    unsigned short modGNeg = (elementSize - modG) % elementSize;
    unsigned short temp = -1;


    for (int i = 0; i < elementSize; i++) {
        masks[i] = temp;
        negMasks[i] = ~temp;
        temp -= 1 << (i);
        //printf("%s %s\n", toBinary(masks[i], elementSize), toBinary(negMasks[i], elementSize));
    }

    //b2[m - 1] += (b2[0] & masks[modG]) >> modG;


    if (inda[0] < elementSize) {
        if (inda[arrSize - 1] >= elementSize) {
            int i = 0;
            for (i = 0; inda[i] < elementSize; i++) {
                unsigned short mod = inda[i] % elementSize;
                unsigned short modNeg = (elementSize - mod) % elementSize;

                //printf("mod %d modNeg %d\n", mod, modNeg);


                res[0] = res[0] ^ (b2[0] & masks[mod]) >> mod;
                /*printf("start\n");
                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }*/
                /*printf("\n");*/
                unsigned short it = 1;
                for (int j = 1; j < m - 1; j++) {
                    /*printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                           toBinary((b2[j] & masks[mod]), elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]), elementSize));*/
                    res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
                    it++;
                }

                /*for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");*/

                // todo подумать что будет с окончаниями элементов
                unsigned short temp2;
                temp2 = ((b2[m - 1] >> mod) + ((b2[m - 2] & negMasks[mod]) << modNeg));
                res[it] = res[it] ^ temp2;

                /*printf("%s %s\n", toBinary((b2[m - 1] >> modGNeg), elementSize),
                       toBinary(((b2[m - 2] & negMasks[modGNeg]) << modG), elementSize));*/

                unsigned short temp3 = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);
                //printf("%s\n", toBinary(temp3, elementSize));
                unsigned short temp4 = (temp3 & negMasks[mod]) << modNeg;
                res[0] = res[0] ^ temp4;


            }




            /**for по элементам, которые больше elementSize*/
            unsigned short lastElementShifted = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);
            //printf("lastElementShifted=%s\n", toBinary(lastElementShifted, elementSize));
            for (; i < arrSize; i++) {
                unsigned short mod = inda[i] % elementSize;
                unsigned short modNeg = (elementSize - mod) % elementSize;
                unsigned short start = inda[i] / elementSize;

                /*printf("mod %d modNeg %d\n", mod, modNeg);*/

                res[start] =
                        res[start] ^ (((b2[0] & masks[mod]) >> mod) + ((lastElementShifted & negMasks[mod]) << modNeg));
                /*printf("start\n");
                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");*/
                unsigned short it = start + 1;
                int j = 1;
                for (; it < m; it++, j++) {
                    /*printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                           toBinary((b2[j] & masks[mod]), elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]), elementSize));*/
                    res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
                }
                j--;

                /*printf("after first cycle \n");
                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");*/



                //todo выдаст ошибку при modG == 0
                short newMod = ((short) mod) - modG;
                int flag = 0;

                if (newMod < 0 || modG == 0) {
                    flag++;
                    /* printf("inside\n");*/
                    j++;
                    newMod = (newMod + elementSize) % elementSize;

                }

                /*printf("start %d\n", start);*/

                unsigned short newModNeg = (elementSize - newMod) % elementSize;
                /*printf("oldMod %d Gmod %d\n", mod, modG);
                printf("newMod %d newModNeg %d\n", newMod, newModNeg);
                printf("j = %d it = %d \n", j, it);*/
                it = 0;
                for (; j < m - 1; j++, it++) {
                    /*printf("for element № %d\ntakes element № %d with value %s and mask %s , what gets %s\n"
                           "takes element № %d with value %s and mask %s , what gets %s\n"
                           "we gets %s\n\n",
                           it, j, toBinary(b2[j], elementSize), toBinary(masks[newMod], elementSize),
                           toBinary(((b2[j] & masks[newMod]) >> newMod), elementSize),
                           j - 1, toBinary(b2[j - 1], elementSize), toBinary(negMasks[newMod], elementSize),
                           toBinary(((b2[j - 1] & negMasks[newMod]) << newModNeg), elementSize),
                           toBinary((((b2[j] & masks[newMod]) >> newMod) +
                                     ((b2[j - 1] & negMasks[newMod]) << newModNeg)), elementSize));*/
                    res[it] =
                            res[it] ^
                            (((b2[j] & masks[newMod]) >> newMod) + ((b2[j - 1] & negMasks[newMod]) << newModNeg));

                    /*for (int z = 0; z < m; z++) {
                        printf("%s ", toBinary(res[z], elementSize));
                    }
                    printf("\n");*/
                }

                if (flag) {
                    res[it] =
                            res[it] ^
                            (((b2[j] & masks[newMod]) >> newMod) + ((b2[j - 1] & negMasks[newMod]) << newModNeg));
                    it++;
                }
            }


        } else {
            for (int i = 0; i < arrSize; i++) {
                unsigned short mod = inda[i] % elementSize;
                unsigned short modNeg = (elementSize - mod) % elementSize;

                //printf("mod %d modNeg %d\n", mod, modNeg);


                res[0] = res[0] ^ (b2[0] & masks[mod]) >> mod;
                /*printf("start\n");
                for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }*/
                /*printf("\n");*/
                unsigned short it = 1;
                for (int j = 1; j < m - 1; j++) {
                    /*printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                           toBinary((b2[j] & masks[mod]), elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                           toBinary((b2[j - 1] & negMasks[mod]), elementSize));*/
                    res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
                    it++;
                }

                /*for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");*/

                // todo подумать что будет с окончаниями элементов
                unsigned short temp2;
                temp2 = ((b2[m - 1] >> mod) + ((b2[m - 2] & negMasks[mod]) << modNeg));
                res[it] = res[it] ^ temp2;

                /*printf("%s %s\n", toBinary((b2[m - 1] >> modGNeg), elementSize),
                       toBinary(((b2[m - 2] & negMasks[modGNeg]) << modG), elementSize));*/

                unsigned short temp3 = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);
                /*printf("%s\n", toBinary(temp3, elementSize));*/
                unsigned short temp4 = (temp3 & negMasks[mod]) << modNeg;
                res[0] = res[0] ^ temp4;
            }

        }
    } else {

        unsigned short lastElementShifted = (b2[m - 1] >> modGNeg) + ((b2[m - 2] & negMasks[modGNeg]) << modG);
        //printf("lastElementShifted=%s\n", toBinary(lastElementShifted, elementSize));
        for (int i = 0; i < arrSize; i++) {
            unsigned short mod = inda[i] % elementSize;
            unsigned short modNeg = (elementSize - mod) % elementSize;
            unsigned short start = inda[i] / elementSize;

            /*printf("mod %d modNeg %d\n", mod, modNeg);*/

            res[start] =
                    res[start] ^ (((b2[0] & masks[mod]) >> mod) + ((lastElementShifted & negMasks[mod]) << modNeg));
            /*printf("start\n");
            for (int z = 0; z < m; z++) {
                printf("%s ", toBinary(res[z], elementSize));
            }
            printf("\n");*/
            unsigned short it = start + 1;
            int j = 1;
            for (; it < m; it++, j++) {
                /*printf("%s %s %s %s\n", toBinary(((b2[j] & masks[mod]) >> mod), elementSize),
                       toBinary((b2[j] & masks[mod]), elementSize),
                       toBinary((b2[j - 1] & negMasks[mod]) << modNeg, elementSize),
                       toBinary((b2[j - 1] & negMasks[mod]), elementSize));*/
                res[it] = res[it] ^ (((b2[j] & masks[mod]) >> mod) + ((b2[j - 1] & negMasks[mod]) << modNeg));
            }
            j--;

            /*printf("after first cycle \n");
            for (int z = 0; z < m; z++) {
                printf("%s ", toBinary(res[z], elementSize));
            }
            printf("\n");*/



            //todo выдаст ошибку при modG == 0
            short newMod = ((short) mod) - modG;
            int flag = 0;

            if (newMod < 0 || modG == 0) {
                flag++;
                /* printf("inside\n");*/
                j++;
                newMod = (newMod + elementSize) % elementSize;

            }

            /*printf("start %d\n", start);*/

            unsigned short newModNeg = (elementSize - newMod) % elementSize;
            /*printf("oldMod %d Gmod %d\n", mod, modG);
            printf("newMod %d newModNeg %d\n", newMod, newModNeg);
            printf("j = %d it = %d \n", j, it);*/
            it = 0;
            for (; j < m - 1; j++, it++) {
                /*printf("for element № %d\ntakes element № %d with value %s and mask %s , what gets %s\n"
                       "takes element № %d with value %s and mask %s , what gets %s\n"
                       "we gets %s\n\n",
                       it, j, toBinary(b2[j], elementSize), toBinary(masks[newMod], elementSize),
                       toBinary(((b2[j] & masks[newMod]) >> newMod), elementSize),
                       j - 1, toBinary(b2[j - 1], elementSize), toBinary(negMasks[newMod], elementSize),
                       toBinary(((b2[j - 1] & negMasks[newMod]) << newModNeg), elementSize),
                       toBinary((((b2[j] & masks[newMod]) >> newMod) +
                                 ((b2[j - 1] & negMasks[newMod]) << newModNeg)), elementSize));*/
                res[it] =
                        res[it] ^
                        (((b2[j] & masks[newMod]) >> newMod) + ((b2[j - 1] & negMasks[newMod]) << newModNeg));

                /*for (int z = 0; z < m; z++) {
                    printf("%s ", toBinary(res[z], elementSize));
                }
                printf("\n");*/
            }

            if (flag) {
                res[it] =
                        res[it] ^
                        (((b2[j] & masks[newMod]) >> newMod) + ((b2[j - 1] & negMasks[newMod]) << newModNeg));
                it++;
            }
        }
    }


    /*todo обрезка вектора тут не к чему, но создаст потенциально проблемы в скорости*/
    /*unsigned short t;
    for (int i = 0; i < modG; i++)
        t += 1 << i;
    t = t << ((elementSize - modG + 1) % elementSize);
    printf("%s\n", toBinary(t,elementSize));
*/
    //res[m - 1] = res[m - 1] & t;

    /*for (int z = 0; z < m; z++) {
        printf("%s ", toBinary(res[z], elementSize));
    }
    printf("\n");*/


}

/**разжимает вектор*/
void decompressArray(unsigned short M, unsigned short m, unsigned short *arrTo, unsigned short *arrFrom) {
    int it = M - 1;
    unsigned short temp = arrFrom[m - 1];
    unsigned short mod = M % elementSize;
    temp = temp >> (elementSize - (M % elementSize)) % elementSize;

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

int main2() {
    srand(time(0));

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

    /**компактное хранение h1*/
    unsigned short h1Compact[hLength];

    /**компактное хранение h2*/
    unsigned short h2Compact[hLength];

    /**компактное хранение h1 в обратном порядке*/
    unsigned short h1TransCompact[hLength];

    /**компактное хранение h2 в обратном порядке*/
    unsigned short h2TransCompact[hLength];

    for (int i = 0; i < hLength; i++) {
        h1Compact[i] = h2Compact[i] = h1TransCompact[i] = h2TransCompact[i] = 0;
    }

    /**компактное хранение e1*/
    unsigned short e1Compact[eLength];

    /**компактное хранение e2*/
    unsigned short e2Compact[eLength];

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
    for (int i = 0; i < n; i++) {
        c1[i] = c2[i] = s[i] = 0;
    }


    computationF2(n, hLength, h1Compact, e1, c1);
    computationF2(n, hLength, h2Compact, e2, c2);
    for (int i = 0; i < n; i++) {
        s[i] = c1[i] ^ c2[i];
    }

    /**посчитали S*/
    /*for (int i = 0; i < 150; i++)
        printf("%d", s[i]);



    printf("\n");*/

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
        computationZ(n, hLength, h1TransCompact, sTemp, upc1);
        computationZ(n, hLength, h2TransCompact, sTemp, upc2);
        for (int j = 0; j < n; j++) {
            if (upc1[j] >= T)
                u[j] = u[j] ^ 1;
            if (upc2[j] >= T)
                v[j] = v[j] ^ 1;
        }

        for (int i = 0; i < n; i++) {
            sTemp[i] = c1[i] = c2[i] = 0;
        }
        computationF2(n, hLength, h1Compact, u, c1);
        computationF2(n, hLength, h2Compact, v, c2);
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
        printf("success\n");

    } else {
        printf("fail\n");
    }
    return 0;
}

int main3() {

    srand(time(0));

    unsigned short masks[elementSize];
    masks[0] = 1 << (elementSize - 1);
    for (int i = 1; i < elementSize; i++)
        masks[i] = masks[i - 1] >> 1;

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

    /**заполняем вектора данными*/
    generateSparseArray(h1Compact, n, hLength);
    generateSparseArray(h2Compact, n, hLength);
    generateSparseArray(e1Compact, n, eLength);
    generateSparseArray(e2Compact, n, eLength);


    for (int i = 0; i < hLength; i++) {
        h1TransCompact[i] = n - h1Compact[i];
        h2TransCompact[i] = n - h2Compact[i];
    }

    unsigned short e1[n];
    unsigned short e2[n];
    for (int i = 0; i < n; i++)
        e1[i] = e2[i] = 0;
    for (int i = 0; i < eLength; i++)
        e1[e1Compact[i]] = 1;
    for (int i = 0; i < eLength; i++)
        e2[e2Compact[i]] = 1;

    unsigned short e1Compressed[m];
    unsigned short e2Compressed[m];
    for(int i = 0;i<m;i++)
    {
        e1Compressed[i]=e2Compressed[i]=0;
    }

    createDenseArray(e1Compressed,m,e1,n);
    createDenseArray(e2Compressed,m,e2,n);

    unsigned short c1[m];
    unsigned short c2[m];
    unsigned short s[m];
    for (int i = 0; i < m; i++) {
        c1[i] = c2[i] = s[i] = 0;
    }

    //нужны для сравнения
    unsigned short c1_2[n];
    unsigned short c2_2[n];
    for (int i = 0; i < n; i++) {
        c1_2[i] = c2_2[i] = 0;
    }

    calculateSparseAndUsual2(n, hLength, m, h1Compact, e1Compressed, c1);
    computationF2(n, hLength, h1Compact, e1, c1_2);

    if (!compare(n, c1_2, m, c1)) {
        printf("not Equal\n");
        for (int i = 0; i < 10; i++) {
            printf("%s ", toBinary(c1[i], elementSize));

        }
        printf("\n");

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 16; j++)
                printf("%d", c1_2[i * 16 + j]);
            printf(" ");
        }
        printf("\n");
    }


    calculateSparseAndUsual2(n, hLength, m, h2Compact, e2Compressed, c2);
    computationF2(n, hLength, h2Compact, e2, c2_2);
    if (!compare(n, c2_2, m, c2)) {
        printf("not Equal\n");
    }

    for (int i = 0; i < n; i++) {
        s[i] = c1[i] ^ c2[i];
    }

    /**посчитали S*/
    /*for (int i = 0; i < 150; i++)
        printf("%d", s[i]);
    printf("\n");*/
    printf("s calculated\n");

    unsigned short u_2[n], v_2[n];
    unsigned short u[m], v[m];
    for (int i = 0; i < m; i++)
        u[i] = v[i] = 0;

    unsigned short sTemp[m];
    for (int i = 0; i < m; i++) {
        sTemp[i] = s[i];
    }

    for (int z = 0; z < num_it; z++) {
        unsigned short upc1[n], upc2[n];
        for (int i = 0; i < n; i++)
            upc1[i] = upc2[i] = 0;
        unsigned short sTempDecompressed[n];

        decompressArray(n, m, sTempDecompressed, sTemp);

        computationZ(n, hLength, h1TransCompact, sTempDecompressed, upc1);
        computationZ(n, hLength, h2TransCompact, sTempDecompressed, upc2);

        /**тут пойдут изменения для u и v*/
        for (int j = 0; j < n; j++) {
            if (upc1[j] >= T)
                u[j / elementSize] = u[j / elementSize] ^ masks[j % elementSize];
            if (upc2[j] >= T)
                v[j / elementSize] = v[j / elementSize] ^ masks[j % elementSize];
        }

        for (int i = 0; i < m; i++) {
            sTemp[i] = c1[i] = c2[i] = 0;
        }
        for (int i = 0; i < n; i++) {
            c1_2[i] = c2_2[i] = 0;
        }

        decompressArray(n, m, u_2, u);
        decompressArray(n, m, v_2, v);


        calculateSparseAndUsual2(n, hLength, m, h1Compact, u, c1);
        computationF2(n, hLength, h1Compact, u_2, c1_2);

        if (!compare(n, c1_2, m, c1)) {
            printf("not Equal\n");
        }

        calculateSparseAndUsual2(n, hLength, m, h2Compact, v, c2);
        computationF2(n, hLength, h2Compact, v_2, c2_2);

        if (!compare(n, c2_2, m, c2)) {
            printf("not Equal\n");
        }

        for (int i = 0; i < n; i++) {
            sTemp[i] = s[i] ^ c1[i] ^ c2[i];
        }

        /**проверка s` на ноль*/
        flag = 0;
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
        printf("success\n");
    } else {
        printf("fail\n");
    }
    return 0;
}

int mainSushkoTest() {
    srand(time(0));

    unsigned short masks[elementSize];
    masks[0] = 1 << (elementSize - 1);
    for (int i = 1; i < elementSize; i++)
        masks[i] = masks[i - 1] >> 1;

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

    unsigned short n = 4801;

    unsigned short m = (n + elementSize - 1) / elementSize;

    /**компактное хранение h1*/
    unsigned short h1Compact[hLength];


    /**компактное хранение h1 в обратном порядке*/
    unsigned short h1TransCompact[hLength];


    /**компактное хранение e1*/
    unsigned short e1Compact[eLength];


    /**заполняем вектора данными*/
    generateSparseArray(h1Compact, n, hLength);
    generateSparseArray(e1Compact, n, eLength);


    for (int i = 0; i < hLength; i++) {
        h1TransCompact[i] = n - h1Compact[i];
    }

    unsigned short e1[n];
    for (int i = 0; i < n; i++)
        e1[i] = 0;
    for (int i = 0; i < eLength; i++)
        e1[e1Compact[i]] = 1;

    unsigned short e1Compressed[m];

    createDenseArray(e1Compressed, m, e1, n);
    //decompressArray(n,m,e1,e1Compact);

    unsigned short c1[m];
    for (int i = 0; i < m; i++) {
        c1[i] = 0;
    }

    //нужны для сравнения
    unsigned short c1_2[n];
    unsigned short c2_2[n];
    for (int i = 0; i < n; i++) {
        c1_2[i] = c2_2[i] = 0;
    }

    /*for (int i = 0; i < m; i++) {
        printf("%s ", toBinary(e1Compressed[i], elementSize));
    }
    printf("\n");*/

    calculateSparseAndUsual2(n, hLength, m, h1Compact, e1Compressed, c1);
    computationF2(n, hLength, h1Compact, e1, c1_2);

    if (!compare(n, c1_2, m, c1)) {
        printf("not Equal\n");
        for (int i = 0; i < m; i++) {
            printf("%s ", toBinary(c1[i], elementSize));

        }
        printf("\n");

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < 16; j++)
                printf("%d", c1_2[i * 16 + j]);
            printf(" ");
        }
        printf("\n");
    }


    return 0;
}

int main() {


    srand(time(0));


    /*unsigned short M = 67;
    unsigned short m = (M + elementSize - 1) / elementSize;
    const unsigned short arrSize = 7;
    unsigned short inda[7] = {1,5,17,32,35,64,66};
    unsigned short b[M];
    unsigned short c2[M];
    for (int i = 0; i < M; i++) {
        c2[i] = 0;
    }
    unsigned short b2[m], c[m];
    for (int i = 0; i < m; i++)
        b2[i] = c[i] = 0;

    generateDenseArray(b, M);
    createDenseArray(b2, m, b, M);
    for (int i = 0; i < M; i++) {
        printf("%d", b[i]);
    }
    printf("\n");
    for (int i = 0; i < m; i++) {
        printf("%s ", toBinary(b2[i], elementSize));
    }
    printf("\n");
    printf("\n");

    calculateSparseAndUsual2(M, arrSize, m, inda, b2, c);
    computationF2(M, arrSize, inda, b, c2);

    for (int z = 0; z < m; z++) {
        printf("%s ", toBinary(c[z], elementSize));
    }
    printf("\n");

    for (int i = 0; i < M; i++) {
        if (i != 0 && (i % elementSize) == 0)
            printf(" ");
        printf("%d", c2[i]);
    }
    printf("\n");
    printf("equal = %d\n", compare(M, c2, m, c));*/

    main3();
    /*for (int i = 0; i < 1000; i++)
        main2();
    return 0;*/

    //mainSushkoTest();
}


void
computationF2(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < M; j++) {
            c[j] = (c[j] + b[(j - inda[i] + M) % M]) % 2;
        }
    }
}

void
computationZ(unsigned short M, int arrSize, const unsigned short *inda, const unsigned short *b, unsigned short *c) {
    for (int i = 0; i < arrSize; i++) {
        for (int j = 0; j < M; j++) {
            c[j] = (c[j] + b[(j - inda[i] + M) % M]);
        }
    }
}
