/* * * * * * * * * * * * * * * * * * * * */
/* @author Rodrigo Kenji Asato Kobayashi */
/* @RGA 202319040134                     */
/* * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <time.h>
#include <string.h>

#define EPSILON 0.00001

#ifdef DEBUG
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...) 
#endif

typedef struct _celula{
    int conectado;
    float custo;
}celula;

typedef struct _tGrafo{
    int vertices; // numero de vertices = n
    int arestas; // numero de arestas = m
    celula **MatrizAdj; // matriz de adjacencia
}tGrafo;

typedef struct _t_arco{
    int origem;
    int destino;
}t_arco;

int  carga_instancia(FILE* fin, tGrafo *grafo);
void free_graph_matrix(tGrafo *grafo);
int  carga_lp(glp_prob **lp,int s, int t, tGrafo* grafo);
int  busca_prox(t_arco *arcos, int tam, int s);

int carga_instancia(FILE* fin, tGrafo *grafo)
{
    int q, i, j;
    float custo;

    if (fscanf(fin, "P %d %d\n", &(grafo->vertices), &(grafo->arestas)) != 2) {
        return 0;
    }
    PRINTF("arestas = %d  and vertices = %d\n", grafo->vertices,grafo->arestas);

    // Aloca memoria para a matriz de adjacência
    grafo->MatrizAdj = (celula **)malloc(grafo->vertices * sizeof(celula *));

    for(i = 0; i < grafo->vertices; i++){
        grafo->MatrizAdj[i] = (celula*)malloc((grafo->vertices)*sizeof(celula));
    }

    for(i = 0; i < grafo->vertices ;i++){
        for(j = 0; j < grafo->vertices; j++){
            (grafo->MatrizAdj[i][j]).conectado = 0;
            (grafo->MatrizAdj[i][j]).custo = 0;
        }
    }
    
    // Le todas os arcos e seus custos
    q = 0;
    while (!feof(fin)) {
        if (fscanf(fin, "A %d %d %f\n", &i, &j, &custo) == 3) {
            PRINTF("aresta %d_%d custo = %.3f\n", i, j, custo);
            q++;
            grafo->MatrizAdj[i][j].conectado = 1;
            grafo->MatrizAdj[i][j].custo = custo;
        } else {
            break;
        }
    }
    PRINTF("\n");

    // Caso que numero de arestas lidos eh diferente do informado
    if(q != grafo->arestas){ 
        free_graph_matrix(grafo);
        return 0;
    }
    return 1;
}

void free_graph_matrix(tGrafo *grafo){
    for(int i = 0; i < grafo->vertices; i++){
        free(grafo->MatrizAdj[i]);
    }
    free(grafo->MatrizAdj);
}

int carga_lp(glp_prob **lp,int s, int t, tGrafo* grafo)
{
    int *ia, *ja, rows, cols, i, j, k, nz;
    double *ar;
    char name[80];
    t_arco *arcos;

    rows = grafo->vertices;
    cols = grafo->arestas; 

    // Aloca matriz de coeficientes
    ia = (int*)malloc(sizeof(int)*(2*(grafo->arestas)+1));
    ja = (int*)malloc(sizeof(int)*(2*(grafo->arestas)+1));
    ar = (double*)malloc(sizeof(double)*(2*(grafo->arestas)+1));

    // Aloca vetor de arcos, que guardara os vetores lidos em ordem lexicografica
    arcos = (t_arco*)malloc(sizeof(t_arco)*grafo->arestas);
      
    // Cria problema de PL
    *lp = glp_create_prob();
    glp_set_prob_name(*lp, "caminho_min");
    glp_set_obj_dir(*lp, GLP_MIN);

    // Configura restricoes
    glp_add_rows(*lp, rows);
    for(i = 0; i < rows; i++){
        name[0]='\0';
        sprintf(name,"v_%d", i);
        glp_set_row_name(*lp, i+1, name);
        if(i == s){
            glp_set_row_bnds(*lp, i+1, GLP_FX, 1.0, 1.0);
        }else if(i == t){
            glp_set_row_bnds(*lp, i+1, GLP_FX, -1.0, -1.0);
        }else{
            glp_set_row_bnds(*lp, i+1, GLP_FX, 0.0, 0.0);
        }
    }

    // Configura variaveis
    glp_add_cols(*lp, cols);
    k = 1;
    for(i = 0; i < grafo->vertices; i++){
        for(j = 0; j < grafo->vertices; j++){
            if((grafo->MatrizAdj[i][j]).conectado){
                name[0]='\0';
                sprintf(name,"x_%d_%d", i, j);
                glp_set_col_name(*lp, k, name);
                glp_set_col_bnds(*lp, k, GLP_DB, 0.0, 1.0); // 0<=x_i_j<=1
                glp_set_obj_coef(*lp, k, (grafo->MatrizAdj[i][j]).custo);
                (arcos[k-1]).origem = i; (arcos[k-1]).destino = j; // 
                k++;
            }
        }
    }

    // Configura matriz de coeficientes // vertices x arestas 
    nz = 1;
    for(i = 1; i <= grafo->vertices; i++){
        k = 0;
        for(j = 1; j <= grafo->arestas; j++){
            if(i-1 == arcos[k].origem){
                ia[nz] = i; ja[nz] = j; ar[nz++] =  1.0;
            }else if(i-1 == arcos[k].destino){
                ia[nz] = i; ja[nz] = j; ar[nz++] =  -1.0;
            }
            k++;
        }
    }

    // Carrega PL
    glp_load_matrix(*lp, nz-1, ia, ja, ar);

    // libera memoria
    free(ia); free(ja); free(ar); free(arcos);
    return 1;
}

int busca_prox(t_arco *arcos, int tam, int s){ // s = source, return t = target
    int t = -1, found = 0;
    for(int i = 0; i < tam && !found; i++){
        if(arcos[i].origem == s){
            t = (arcos)[i].destino;
            found = 1;
        }
    }
    return t;
}

int main(int argc, char **argv)
{
    glp_prob *lp;
    double z, valor, tempo;
    FILE *fin;
    int i, j, k, s, t, tam;
    tGrafo grafo;
    clock_t antes, agora;
    char nome[80];
    t_arco *arcos;

    glp_smcp param_lp;

    if(argc < 4){
        PRINTF("Sintaxe: ./caminho <grafo> s t\n");
        printf("E\n");
        exit(1);
    }

    fin = fopen(argv[1], "r");
    if (!fin){
        printf("problema na abertura do arquivo: %s", argv[1]);
        exit(1);
    }

    s = atoi(argv[2]);
    t = atoi(argv[3]);

    if(!carga_instancia(fin, &grafo)){
        printf("E\n");
        exit(1);
    }
    if(!(s < grafo.vertices && t < grafo.vertices) || s == t){ // caso em que s ou t são >= a grafo.vertices (assume-se que os vertices estão nomeados de 0 a n-1)
        printf("E\n");                                         // caso s == t
        exit(1);
    }

    // desabilita saidas do GLPK no terminal
    glp_term_out(GLP_OFF);

    // carga do lp
    carga_lp(&lp, s, t, &grafo); 

    // configura simplex
    glp_init_smcp(&param_lp);
    param_lp.msg_lev = GLP_MSG_OFF;

    antes = clock();
    // Executa Solver de PL
    glp_simplex(lp, &param_lp);
    agora = clock();

    // aloca vetor de arcos para no final imprimir o caminho utilizando busca_prox
    arcos = (t_arco*)malloc(grafo.arestas*sizeof(t_arco));
    tam = 0;

    // Recupera solucao 
    if(glp_get_status(lp) == GLP_OPT){
        z = glp_get_obj_val(lp);
        printf("C %.3f\n", z);
        for(k = 0; k < grafo.arestas; k++){
            valor=glp_get_col_prim(lp, k+1);
            strcpy(nome, glp_get_col_name(lp, k+1));
            sscanf(nome, "x_%d_%d", &i, &j);
            if(valor > EPSILON) {
                printf("V %d %d %.3f\n", i, j, valor); // imprime a variavel não nula
                arcos[tam].origem = i; arcos[tam].destino = j; tam++; // guarda tal arco em arcos
            }
        }

        // Imprime Caminho comecando em s ate t
        printf("P ");
        k = s;
        do{
            printf("%d ", k);
            k = busca_prox(arcos, tam, k);
        }while(k != t && k >= 0);
        printf("%d\n", k);

    }else if(glp_get_status(lp) == GLP_NOFEAS){ // nao tem solucao
        printf("I\n");
    }

#ifdef DEBUG	
      // Grava solucao e PL
      PRINTF("\n---LP gravado em caminho.lp e solucao em caminho.sol");
      glp_write_lp(lp, NULL,"caminho.lp");
      glp_print_sol(lp, "caminho.sol");
#endif
    tempo = ((double)agora-antes)/CLOCKS_PER_SEC;
    PRINTF("\n\n**** z =%.2lf tempo=%.2lf\n", z, tempo);

    // libera memoria alocada
    free_graph_matrix(&grafo); free(arcos);

    // Destroi problema
    glp_delete_prob(lp);
    return 0;
}

/* eof */
