#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define A 0 // A
#define T 1
#define G 2
#define C 3
#define sair 14

#define maxSeq 10000
#define maxNumThreads 1000

char baseMapa[5] = {'A', 'T', 'G', 'C', '-'};

char seqMaior[maxSeq], seqMenor[maxSeq], alinhaMaior[maxSeq][maxSeq], alinhaMenor[maxSeq][maxSeq];

int matrizScores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior = 6, tamSeqMenor = 6, tamAlinha[maxSeq], penalGap = 0, grauMuta = 0, diagScore, linScore, colScore, k = 1, numThreads = 1;

int matrizPesos[4][4] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

int alinhamentoScores[maxSeq];

typedef struct
{
    int tbLin;
    int tbCol;
    int pos;
    int thread_id;
} traceback_data_t;

typedef struct
{
    int thread_id;
    int num_threads;
    int lin_start;
    int lin_end;
} thread_data_t;

// Declaração de funções
int leTamMaior(void);
int leTamMenor(void);
int lePenalidade(void);
void leMatrizPesos();
void mostraMatrizPesos(void);
int leGrauMutacao(void);
void leSequenciasDeArquivo();
void leSequencias();
void geraSequencias();
void mostraSequencias(void);
void mostraAlinhamentoGlobal(void);
void salvaMatrizScores(void);
void mostraMatrizScores();
void *preencheMatrizScores(void *threadarg);
void geraMatrizScores(void);
void *traceBackThread(void *threadarg);
void traceBack(int k);
int leNumeroDeAlinhamentos(void);
int leNumeroDeThreads(void);
int menuOpcao(void);
void trataOpcao(int op);

// Implementação das funções
int leTamMaior(void)
{
    printf("\nLeitura do Tamanho da Sequencia Maior:");
    do
    {
        printf("\nDigite 0 < valor < %d = ", maxSeq);
        scanf("%d", &tamSeqMaior);
    } while ((tamSeqMaior < 1) || (tamSeqMaior > maxSeq));
    return tamSeqMaior;
}

int leTamMenor(void)
{
    printf("\nLeitura do Tamanho da Sequencia Menor:");
    do
    {
        printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
        scanf("%d", &tamSeqMenor);
    } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
    return tamSeqMenor;
}

int lePenalidade(void)
{
    int penal;

    printf("\nLeitura da Penalidade de Gap:");
    do
    {
        printf("\nDigite valor >= 0 = ");
        scanf("%d", &penal);
    } while (penal < 0);

    return penal;
}

void leMatrizPesos()
{
    int i, j;

    printf("\nLeitura da Matriz de Pesos:\n");
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            printf("Digite valor %c x %c = ", baseMapa[i], baseMapa[j]);
            scanf("%d", &(matrizPesos[i][j]));
        }
        printf("\n");
    }
}

void mostraMatrizPesos(void)
{
    int i, j;

    printf("\nMatriz de Pesos Atual:");
    printf("\n%4c%4c%4c%4c%4c\n", ' ', 'A', 'T', 'G', 'C');
    for (i = 0; i < 4; i++)
    {
        printf("%4c", baseMapa[i]);
        for (j = 0; j < 4; j++)
            printf("%4d", matrizPesos[i][j]);
        printf("\n");
    }
}

int leGrauMutacao(void)
{
    int prob;

    printf("\nLeitura da Porcentagem Maxima de Mutacao Aleatoria:\n");
    do
    {
        printf("\nDigite 0 <= valor <= 100 = ");
        scanf("%d", &prob);
    } while ((prob < 0) || (prob > 100));

    return prob;
}

void leSequenciasDeArquivo()
{
    FILE *fileMaior, *fileMenor;
    char nomeArquivoMaior[100], nomeArquivoMenor[100];

    printf("\nDigite o nome do arquivo para a Sequencia Maior: ");
    scanf("%s", nomeArquivoMaior);
    fileMaior = fopen(nomeArquivoMaior, "r");
    if (fileMaior == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nomeArquivoMaior);
        exit(1);
    }

    printf("\nDigite o nome do arquivo para a Sequencia Menor: ");
    scanf("%s", nomeArquivoMenor);
    fileMenor = fopen(nomeArquivoMenor, "r");
    if (fileMenor == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nomeArquivoMenor);
        exit(1);
    }

    int oldTamMaior = tamSeqMaior;
    int oldTamMenor = tamSeqMenor;
    char oldMaior[maxSeq] = seqMaior;
    char oldMenor[maxSeq] = seqMenor;

    tamSeqMaior = fread(seqMaior, sizeof(char), maxSeq, fileMaior);
    tamSeqMenor = fread(seqMenor, sizeof(char), maxSeq, fileMenor);

    for (int i = 0; i < tamSeqMaior; i++)
    {
        if (seqMaior[i] == 'A')
        {
            seqMaior[i] = 0;
        }else if (seqMaior[i] == 'T')
        {
            seqMaior[i] = 1;
        }else if (seqMaior[i] == 'G')
        {
            seqMaior[i] = 2;
        }else if (seqMaior[i] == 'C'){
            seqMaior[i] = 3;

        }else 
        {
            tamSeqMaior = oldTamMaior;
            seqMaior = oldMaior;
            printf("Sequencia maior invalida, valor original restaurado.\n");
            exit(1);
        }
    }

    for (int i = 0; i < tamSeqMenor; i++)
    {
        if (seqMenor[i] == 'A')
        {
            seqMenor[i] = 0;
        }else if (seqMenor[i] == 'T')
        {
            seqMenor[i] = 1;
        }else if (seqMenor[i] == 'G')
        {
            seqMenor[i] = 2;
        }else if (seqMenor[i] == 'C'){
            seqMenor[i] = 3;

        }else 
        {
            tamSeqMenor = oldTamMenor;
            seqMenor = oldMenor;
            printf("Sequencia menor invalida, valor original restaurado. \n");
            exit(1);
        }
    }

    fclose(fileMaior);
    fclose(fileMenor);

    if (tamSeqMenor > tamSeqMaior)
    {
        printf("Erro: Sequencia menor é maior que a sequencia maior.\n");
        exit(1);
    }
}

void leSequencias()
{
    int i, erro, opcao;

    printf("\nSelecione o modo de entrada das sequencias:\n");
    printf("1. Manual\n");
    printf("2. Arquivo\n");
    printf("3. Aleatoria\n");
    scanf("%d", &opcao);

    if (opcao == 1)
    {
        printf("\nLeitura das Sequencias:\n");
        do
        {
            printf("\nPara a Sequencia Maior,");
            printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
            do
            {
                printf("\n> ");
                fgets(seqMaior, maxSeq, stdin);
                tamSeqMaior = strlen(seqMaior) - 1;
            } while (tamSeqMaior < 1);
            printf("\ntamSeqMaior = %d\n", tamSeqMaior);
            i = 0;
            erro = 0;
            do
            {
                switch (seqMaior[i])
                {
                case 'A':
                    seqMaior[i] = (char)A;
                    break;
                case 'T':
                    seqMaior[i] = (char)T;
                    break;
                case 'G':
                    seqMaior[i] = (char)G;
                    break;
                case 'C':
                    seqMaior[i] = (char)C;
                    break;
                default:
                    erro = 1;
                }
                i++;
            } while ((erro == 0) && (i < tamSeqMaior));
        } while (erro == 1);

        do
        {
            printf("\nPara a Sequencia Menor, ");
            printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
            do
            {
                printf("\n> ");
                fgets(seqMenor, maxSeq, stdin);
                tamSeqMenor = strlen(seqMenor) - 1;
            } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
            printf("\ntamSeqMenor = %d\n", tamSeqMenor);

            i = 0;
            erro = 0;
            do
            {
                switch (seqMenor[i])
                {
                case 'A':
                    seqMenor[i] = (char)A;
                    break;
                case 'T':
                    seqMenor[i] = (char)T;
                    break;
                case 'G':
                    seqMenor[i] = (char)G;
                    break;
                case 'C':
                    seqMenor[i] = (char)C;
                    break;
                default:
                    erro = 1;
                }
                i++;
            } while ((erro == 0) && (i < tamSeqMenor));
        } while (erro == 1);
    }
    else if (opcao == 2)
    {
        leSequenciasDeArquivo();
    }
    else if (opcao == 3)
    {
        leTamMaior();
        leTamMenor();
        grauMuta = leGrauMutacao();
        geraSequencias();
    }
    else
    {
        printf("Opcao invalida!\n");
        exit(1);
    }
}

void geraSequencias()
{
    int i, dif, probAux, ind, nTrocas;
    char base;

    srand(time(NULL));

    printf("\nGeracao Aleatoria das Sequencias:\n");

    for (i = 0; i < tamSeqMaior; i++)
    {
        base = (char)(rand() % 4);
        seqMaior[i] = base;
    }

    dif = tamSeqMaior - tamSeqMenor;
    ind = 0;
    if (dif > 0)
        ind = rand() % dif;

    for (i = 0; i < tamSeqMenor; i++)
        seqMenor[i] = seqMaior[ind + i];

    i = 0;
    nTrocas = 0;
    while ((i < tamSeqMenor) && (nTrocas < grauMuta))
    {
        probAux = rand() % 100;

        if (probAux < grauMuta)
        {
            seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
            nTrocas++;
        }
    }

    printf("\nSequencias Geradas, Dif = %d Ind = %d\n", dif, ind);
}

void mostraSequencias(void)
{
    int i;

    printf("\nSequencias Atuais:\n");
    printf("\nSequencia Maior, Tam = %d\n", tamSeqMaior);
    printf("%c", baseMapa[(int)seqMaior[0]]);
    for (i = 1; i < tamSeqMaior; i++)
        printf("%c", baseMapa[(int)seqMaior[i]]);
    printf("\n");

    printf("\nSequencia Menor, Tam = %d\n", tamSeqMenor);
    printf("%c", baseMapa[(int)seqMenor[0]]);
    for (i = 1; i < tamSeqMenor; i++)
        printf("%c", baseMapa[(int)seqMenor[i]]);
    printf("\n");
}

void mostraAlinhamentoGlobal(void)
{
    int i, j;
    int bestThread = 0;
    int bestScore = alinhamentoScores[0];

    printf("\nAlinhamentos Atuais - Mostrando %d Alinhamentos:\n", k);

    // Mostrar todos os alinhamentos gerados pelas threads
    for (j = 0; j < k; j++)
    {
        printf("\nThread %d -> Alinhamento %d (Score: %d):\n", j + 1, j + 1, alinhamentoScores[j]);
        printf("%c", baseMapa[(int)alinhaMaior[j][0]]);
        for (i = 1; i < tamAlinha[j]; i++)
            printf("%c", baseMapa[(int)alinhaMaior[j][i]]);
        printf("\n");

        printf("%c", baseMapa[(int)alinhaMenor[j][0]]);
        for (i = 1; i < tamAlinha[j]; i++)
            printf("%c", baseMapa[(int)alinhaMenor[j][i]]);
        printf("\n");

        // Verificar o melhor score
        if (alinhamentoScores[j] > bestScore)
        {
            bestScore = alinhamentoScores[j];
            bestThread = j;
        }
    }

    // Mostrar o melhor alinhamento
    printf("\nMelhor Alinhamento (Thread %d - Score: %d):\n", bestThread + 1, bestScore);
    printf("%c", baseMapa[(int)alinhaMaior[bestThread][0]]);
    for (i = 1; i < tamAlinha[bestThread]; i++)
        printf("%c", baseMapa[(int)alinhaMaior[bestThread][i]]);
    printf("\n");

    printf("%c", baseMapa[(int)alinhaMenor[bestThread][0]]);
    for (i = 1; i < tamAlinha[bestThread]; i++)
        printf("%c", baseMapa[(int)alinhaMenor[bestThread][i]]);
    printf("\n");
}

void salvaMatrizScores(void)
{
    FILE *file;
    char nomeArquivo[100];
    int i, lin, col;

    printf("\nDigite o nome do arquivo para salvar a Matriz de Scores: ");
    scanf("%s", nomeArquivo);
    file = fopen(nomeArquivo, "w");
    if (file == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nomeArquivo);
        exit(1);
    }

    fprintf(file, "Matriz de Scores:\n");
    fprintf(file, "%4c%4c", ' ', ' ');
    for (i = 0; i <= tamSeqMaior; i++)
        fprintf(file, "%4d", i);
    fprintf(file, "\n");

    fprintf(file, "%4c%4c%4c", ' ', ' ', '-');
    for (i = 0; i < tamSeqMaior; i++)
        fprintf(file, "%4c", baseMapa[(int)(seqMaior[i])]);
    fprintf(file, "\n");

    fprintf(file, "%4c%4c", '0', '-');
    for (col = 0; col <= tamSeqMaior; col++)
        fprintf(file, "%4d", matrizScores[0][col]);
    fprintf(file, "\n");

    for (lin = 1; lin <= tamSeqMenor; lin++)
    {
        fprintf(file, "%4d%4c", lin, baseMapa[(int)(seqMenor[lin - 1])]);
        for (col = 0; col <= tamSeqMaior; col++)
        {
            fprintf(file, "%4d", matrizScores[lin][col]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void mostraMatrizScores()
{
    int i, lin, col;

    printf("\nMatriz de Scores Atual:\n");

    printf("%4c%4c", ' ', ' ');
    for (i = 0; i <= tamSeqMaior; i++)
        printf("%4d", i);
    printf("\n");

    printf("%4c%4c%4c", ' ', ' ', '-');
    for (i = 0; i < tamSeqMaior; i++)
        printf("%4c", baseMapa[(int)(seqMaior[i])]);
    printf("\n");

    printf("%4c%4c", '0', '-');
    for (col = 0; col <= tamSeqMaior; col++)
        printf("%4d", matrizScores[0][col]);
    printf("\n");

    for (lin = 1; lin <= tamSeqMenor; lin++)
    {
        printf("%4d%4c", lin, baseMapa[(int)(seqMenor[lin - 1])]);
        for (col = 0; col <= tamSeqMaior; col++)
        {
            printf("%4d", matrizScores[lin][col]);
        }
        printf("\n");
    }
}

void *preencheMatrizScores(void *threadarg)
{
    int lin, col, peso;
    thread_data_t *data = (thread_data_t *)threadarg;
    int lin_start = data->lin_start;
    int lin_end = data->lin_end;

    for (lin = lin_start; lin < lin_end; lin++)
    {
        for (col = 1; col <= tamSeqMaior; col++)
        {
            peso = matrizPesos[(int)(seqMenor[lin - 1])][(int)(seqMaior[col - 1])];
            diagScore = matrizScores[lin - 1][col - 1] + peso;
            linScore = matrizScores[lin][col - 1] - penalGap;
            colScore = matrizScores[lin - 1][col] - penalGap;

            if ((diagScore >= linScore) && (diagScore >= colScore))
                matrizScores[lin][col] = diagScore;
            else if (linScore > colScore)
                matrizScores[lin][col] = linScore;
            else
                matrizScores[lin][col] = colScore;
        }
    }
    pthread_exit(NULL);
}

void geraMatrizScores(void)
{
    int lin, col, linMaior, colMaior, maior;
    pthread_t threads[numThreads];
    thread_data_t thread_data[numThreads];
    int rc;
    long t;
    void *status;

    printf("\nGeracao da Matriz de Scores:\n");

    for (col = 0; col <= tamSeqMaior; col++)
        matrizScores[0][col] = -1 * (col * penalGap);

    for (lin = 0; lin <= tamSeqMenor; lin++)
        matrizScores[lin][0] = -1 * (lin * penalGap);

    for (t = 0; t < numThreads; t++)
    {
        thread_data[t].thread_id = t;
        thread_data[t].num_threads = numThreads;
        thread_data[t].lin_start = 1 + t * (tamSeqMenor / numThreads);
        thread_data[t].lin_end = 1 + (t + 1) * (tamSeqMenor / numThreads);
        if (thread_data[t].lin_end > tamSeqMenor + 1)
            thread_data[t].lin_end = tamSeqMenor + 1;
        rc = pthread_create(&threads[t], NULL, preencheMatrizScores, (void *)&thread_data[t]);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    for (t = 0; t < numThreads; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
    }

    linMaior = 1;
    colMaior = 1;
    maior = matrizScores[1][1];
    for (lin = 1; lin <= tamSeqMenor; lin++)
    {
        for (col = 1; col <= tamSeqMaior; col++)
        {
            if (maior <= matrizScores[lin][col])
            {
                linMaior = lin;
                colMaior = col;
                maior = matrizScores[lin][col];
            }
        }
    }
    printf("\nMatriz de Scores Gerada.");
    printf("\nUltimo Maior Score = %d na celula [%d,%d]", maior, linMaior, colMaior);
}

void *traceBackThread(void *threadarg)
{
    traceback_data_t *data = (traceback_data_t *)threadarg;
    int tbLin = data->tbLin;
    int tbCol = data->tbCol;
    int pos = data->pos;
    int thread_id = data->thread_id;
    int peso;
    int score = 0;

    while (tbLin > 0 && tbCol > 0)
    {
        peso = matrizPesos[(int)(seqMenor[tbLin - 1])][(int)(seqMaior[tbCol - 1])];
        diagScore = matrizScores[tbLin - 1][tbCol - 1] + peso;
        linScore = matrizScores[tbLin][tbCol - 1] - penalGap;
        colScore = matrizScores[tbLin - 1][tbCol] - penalGap;

        if ((diagScore >= linScore) && (diagScore >= colScore))
        {
            alinhaMenor[thread_id][pos] = seqMenor[tbLin - 1];
            alinhaMaior[thread_id][pos] = seqMaior[tbCol - 1];
            tbLin--;
            tbCol--;
            pos++;
            score += peso;
        }
        else if (linScore > colScore)
        {
            alinhaMenor[thread_id][pos] = (char)4;
            alinhaMaior[thread_id][pos] = seqMaior[tbCol - 1];
            tbCol--;
            pos++;
            score -= penalGap;
        }
        else
        {
            alinhaMenor[thread_id][pos] = seqMenor[tbLin - 1];
            alinhaMaior[thread_id][pos] = (char)4;
            tbLin--;
            pos++;
            score -= penalGap;
        }
    }

    while (tbLin > 0)
    {
        alinhaMenor[thread_id][pos] = seqMenor[tbLin - 1];
        alinhaMaior[thread_id][pos] = (char)4;
        tbLin--;
        pos++;
        score -= penalGap;
    }

    while (tbCol > 0)
    {
        alinhaMenor[thread_id][pos] = (char)4;
        alinhaMaior[thread_id][pos] = seqMaior[tbCol - 1];
        tbCol--;
        pos++;
        score -= penalGap;
    }

    tamAlinha[thread_id] = pos;
    alinhamentoScores[thread_id] = score;
    pthread_exit(NULL);
}

void traceBack(int k)
{
    int tbLin = tamSeqMenor;
    int tbCol = tamSeqMaior;
    int pos = 0;
    pthread_t threads[k];
    traceback_data_t traceback_data[k];
    int rc;
    long t;
    void *status;

    for (t = 0; t < k; t++)
    {
        // Alterar o ponto inicial de tbLin e tbCol para gerar alinhamentos diferentes
        traceback_data[t].tbLin = tbLin - t;
        traceback_data[t].tbCol = tbCol - t;
        traceback_data[t].pos = pos;
        traceback_data[t].thread_id = t;
        rc = pthread_create(&threads[t], NULL, traceBackThread, (void *)&traceback_data[t]);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    for (t = 0; t < k; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
    }

    for (int j = 0; j < k; j++)
    {
        for (int i = 0; i < tamAlinha[j] / 2; i++)
        {
            char aux = alinhaMenor[j][i];
            alinhaMenor[j][i] = alinhaMenor[j][tamAlinha[j] - i - 1];
            alinhaMenor[j][tamAlinha[j] - i - 1] = aux;

            aux = alinhaMaior[j][i];
            alinhaMaior[j][i] = alinhaMaior[j][tamAlinha[j] - i - 1];
            alinhaMaior[j][tamAlinha[j] - i - 1] = aux;
        }
    }

    printf("\nAlinhamento Global Gerado.");
}

int leNumeroDeAlinhamentos(void)
{
    int k;
    printf("\nLeitura do Numero de Alinhamentos a Mostrar (k):\n");
    do
    {
        printf("\nDigite valor > 0: ");
        scanf("%d", &k);
    } while (k <= 0);
    return k;
}

int leNumeroDeThreads(void)
{
    int numThreads;
    printf("\nLeitura do Numero de Threads a Utilizar:\n");
    do
    {
        printf("\nDigite valor > 0 e < %d: ", maxNumThreads);
        scanf("%d", &numThreads);
    } while (numThreads <= 0 && numThreads > maxNumThreads);
    return numThreads;
}

int menuOpcao(void)
{
    int op;
    char enter;

    do
    {
        printf("\nMenu de Opcao:");
        printf("\n<01> Ler Matriz de Pesos");
        printf("\n<02> Mostrar Matriz de Pesos");
        printf("\n<03> Ler Penalidade de Gap");
        printf("\n<04> Mostrar Penalidade");
        printf("\n<05> Definir Sequencias Genomicas");
        printf("\n<06> Mostrar Sequencias");
        printf("\n<07> Gerar Matriz de Scores");
        printf("\n<08> Mostrar Matriz de Scores");
        printf("\n<09> Salvar Matriz de Scores");
        printf("\n<10> Gerar Alinhamento Global");
        printf("\n<11> Mostrar Alinhamento Global");
        printf("\n<12> Definir Numero de Threads");
        printf("\n<13> Sair");
        printf("\nDigite a opcao => ");
        scanf("%d", &op);
        scanf("%c", &enter);
    } while ((op < 1) || (op > 13));

    return (op);
}

void trataOpcao(int op)
{
    int resp;
    char enter;

    switch (op)
    {
    case 1:
        leMatrizPesos();
        break;
    case 2:
        mostraMatrizPesos();
        break;
    case 3:
        penalGap = lePenalidade();
        break;
    case 4:
        printf("\nPenalidade = %d", penalGap);
        break;
    case 5:
        printf("\nDeseja Definicao: <1>MANUAL, <2>ARQUIVO ou <3>ALEATORIA? = ");
        scanf("%d", &resp);
        scanf("%c", &enter); /* remove o enter */
        if (resp == 1)
        {
            leSequencias();
        }
        else if (resp == 2)
        {
            leSequenciasDeArquivo();
        }
        else if (resp == 3)
        {
            leTamMaior();
            leTamMenor();
            grauMuta = leGrauMutacao();
            geraSequencias();
        }
        break;
    case 6:
        mostraSequencias();
        break;
    case 7:
        geraMatrizScores();
        break;
    case 8:
        mostraMatrizScores();
        break;
    case 9:
        salvaMatrizScores();
        break;
    case 10:
        k = leNumeroDeAlinhamentos();
        traceBack(k);
        break;
    case 11:
        mostraAlinhamentoGlobal();
        break;
    case 12:
        numThreads = leNumeroDeThreads();
        break;
    case 13:
        exit(0);
    }
}

int main(void)
{
    int opcao;

    do
    {
        printf("\n\nPrograma Needleman-Wunsch Sequencial\n");
        opcao = menuOpcao();
        trataOpcao(opcao);
    } while (opcao != 14);

    return 0;
}
