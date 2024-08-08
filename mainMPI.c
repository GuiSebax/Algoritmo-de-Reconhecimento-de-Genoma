#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define A 0
#define T 1
#define G 2
#define C 3

#define maxSeq 10000

char baseMapa[5] = {'A', 'T', 'G', 'C', '-'};

char seqMaior[maxSeq], seqMenor[maxSeq];
int matrizScores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior = 6, tamSeqMenor = 6, penalGap = 0;
int matrizPesos[4][4] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

void leMatrizPesos();
void mostraMatrizPesos(void);
void leSequenciasDeArquivo();
void leSequencias();
void geraSequencias();
void mostraSequencias(void);
void mostraAlinhamentoGlobal(void);
void salvaMatrizScores(void);
void mostraMatrizScores();
void geraMatrizScores(int np, int rank, int blocoTamanho);
void traceBack(int rank);
int leTamMaior(void);
int leTamMenor(void);
int lePenalidade(void);
int leTamanhoBloco(void);
int menuOpcao(void);
void trataOpcao(int op, int np, int rank);

int main(int argc, char **argv)
{
    int opcao;
    int np, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    do
    {
        if (rank == 0)
        {
            printf("\n\nPrograma Needleman-Wunsch Paralelo com MPI\n");
            opcao = menuOpcao();
        }

        // Distribuindo a opção para todos os processos
        if (rank == 0)
        {
            for (int i = 1; i < np; i++)
            {
                MPI_Send(&opcao, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(&opcao, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        trataOpcao(opcao, np, rank);
    } while (opcao != 11);

    MPI_Finalize();
    return 0;
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

    tamSeqMaior = fread(seqMaior, sizeof(char), maxSeq, fileMaior);
    tamSeqMenor = fread(seqMenor, sizeof(char), maxSeq, fileMenor);

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
    int i, erro;

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
        i++;
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
    printf("\nAlinhamento Global Gerado.\n");
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

int leTamanhoBloco(void)
{
    int blocoTamanho;

    printf("\nDigite o tamanho do bloco para comunicação MPI: ");
    do
    {
        scanf("%d", &blocoTamanho);
    } while (blocoTamanho <= 0 || blocoTamanho > tamSeqMaior);

    return blocoTamanho;
}

void geraMatrizScores(int np, int rank, int blocoTamanho)
{
    int lin, col, peso;
    int lin_start = rank + 1;
    int localScore[maxSeq + 1];

    if (rank == 0)
    {
        printf("\nGeracao da Matriz de Scores:\n");
        for (col = 0; col <= tamSeqMaior; col++)
            matrizScores[0][col] = -1 * (col * penalGap);
        for (lin = 0; lin <= tamSeqMenor; lin++)
            matrizScores[lin][0] = -1 * (lin * penalGap);
    }

    while (lin_start <= tamSeqMenor)
    {
        for (col = 1; col <= tamSeqMaior; col++)
        {
            peso = matrizPesos[(int)(seqMenor[lin_start - 1])][(int)(seqMaior[col - 1])];
            int diagScore = matrizScores[lin_start - 1][col - 1] + peso;
            int linScore = matrizScores[lin_start][col - 1] - penalGap;
            int colScore = matrizScores[lin_start - 1][col] - penalGap;

            if ((diagScore >= linScore) && (diagScore >= colScore))
                localScore[col] = diagScore;
            else if (linScore > colScore)
                localScore[col] = linScore;
            else
                localScore[col] = colScore;
        }

        for (col = 1; col <= tamSeqMaior; col += blocoTamanho)
        {
            int tamanhoEnvio = (col + blocoTamanho > tamSeqMaior) ? tamSeqMaior - col + 1 : blocoTamanho;
            MPI_Send(&localScore[col], tamanhoEnvio, MPI_INT, 0, lin_start, MPI_COMM_WORLD);
        }

        lin_start += np;
    }

    if (rank == 0)
    {
        for (lin = 1; lin <= tamSeqMenor; lin++)
        {
            for (col = 1; col <= tamSeqMaior; col += blocoTamanho)
            {
                int tamanhoRecebimento = (col + blocoTamanho > tamSeqMaior) ? tamSeqMaior - col + 1 : blocoTamanho;
                MPI_Recv(&matrizScores[lin][col], tamanhoRecebimento, MPI_INT, MPI_ANY_SOURCE, lin, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        printf("\nMatriz de Scores Gerada.\n");
    }
}

void traceBack(int rank)
{
    if (rank == 0)
    {
        int tbLin = tamSeqMenor;
        int tbCol = tamSeqMaior;

        while (tbLin > 0 && tbCol > 0)
        {
            int peso = matrizPesos[(int)(seqMenor[tbLin - 1])][(int)(seqMaior[tbCol - 1])];
            int diagScore = matrizScores[tbLin - 1][tbCol - 1] + peso;
            int linScore = matrizScores[tbLin][tbCol - 1] - penalGap;
            int colScore = matrizScores[tbLin - 1][tbCol] - penalGap;

            if ((diagScore >= linScore) && (diagScore >= colScore))
            {
                tbLin--;
                tbCol--;
            }
            else if (linScore > colScore)
            {
                tbCol--;
            }
            else
            {
                tbLin--;
            }
        }

        printf("\nAlinhamento Global Gerado.\n");
    }
}

int menuOpcao(void)
{
    int op;

    printf("\nMenu de Opcao:");
    printf("\n<01> Ler Matriz de Pesos");
    printf("\n<02> Mostrar Matriz de Pesos");
    printf("\n<03> Ler Penalidade de Gap");
    printf("\n<04> Definir Sequencias Genomicas");
    printf("\n<05> Mostrar Sequencias");
    printf("\n<06> Gerar Matriz de Scores");
    printf("\n<07> Mostrar Matriz de Scores");
    printf("\n<08> Salvar Matriz de Scores");
    printf("\n<09> Gerar Alinhamento Global");
    printf("\n<10> Mostrar Alinhamento Global");
    printf("\n<11> Sair");
    printf("\nDigite a opcao => ");
    scanf("%d", &op);

    return (op);
}

void trataOpcao(int op, int np, int rank)
{
    int blocoTamanho;

    switch (op)
    {
    case 1:
        leMatrizPesos();
        break;
    case 2:
        if (rank == 0)
            mostraMatrizPesos();
        break;
    case 3:
        penalGap = lePenalidade();
        break;
    case 4:
        leSequencias();
        break;
    case 5:
        if (rank == 0)
            mostraSequencias();
        break;
    case 6:
        if (rank == 0)
        {
            blocoTamanho = leTamanhoBloco();
            for (int i = 1; i < np; i++)
            {
                MPI_Send(&blocoTamanho, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(&blocoTamanho, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        geraMatrizScores(np, rank, blocoTamanho);
        break;
    case 7:
        if (rank == 0)
            mostraMatrizScores();
        break;
    case 8:
        if (rank == 0)
            salvaMatrizScores();
        break;
    case 9:
        traceBack(rank);
        break;
    case 10:
        if (rank == 0)
            mostraAlinhamentoGlobal();
        break;
    case 11:
        exit(0);
    }
}
