#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//---------------------------------------------------------------------------------------------------------------------------------------------
/* struct vert, guarda as duas coordenadas de um vertice*/
typedef struct _vertice{
    double x;
    double y;
}vert;
//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
/* struct pol, guarda o numero de vertices e o endereco de onde teremos um vetor de vertices*/
typedef struct _poligono{
    int n; //numero de vertices
    vert *ponto; //ponto[0].x  ponto[0].y - guarda o endereco de algum lugar que vai ter o vetor de pontos
}pol;
//---------------------------------------------------------------------------------------------------------------------------------------------
typedef struct _nfp{
    int partes; //numero de partes do nfp
    pol *polig; //cada parte eh um poligono
}nfp;
//---------------------------------------------------------------------------------------------------------------------------------------------
/* struct item, guarda todas as informacoes necessarias de cada item*/
typedef struct _item{
    int partes; //numero de partes do item
    pol *polig; //cada parte do item eh um poligono
    double area; //area do item
    double esquerda; //ponto mais a esquerda
    double direita; //ponto mais a direita
    double alto; //ponto mais alto
    double baixo; //ponto mais baixo
    int copias; //demanda do item
    nfp *nfps; //guarda o endereco de um lugar que vai ter uma struct tipo pol (nfp desse item com cada item)
    vert *pfinal; //guarda o endereco de um lugar que vai ter uma struct tipo vert (onde cada copia do item foi alocada)
    int cop_alocad; //quantidade de quantas copias ja foram alocadas
    int num_rotacoes;
    double *rotacoes;
}item;  
//---------------------------------------------------------------------------------------------------------------------------------------------
//=============================================================================================================================================
/* funcao para elevar um valor ao quadrado*/
double sq(double x){
    return (x * x);
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao que abre todos os arquivos de NFPs e aloca as informacoes recebidas*/
int LerNFPs(item *itens, int numitens){
    char lixo;
    char nome[35];
    int numvert, i, j, k, m;

    FILE *lenfp;

    for(i = 0; i<numitens; i++){        
        itens[i].nfps = (nfp*)malloc((numitens)*sizeof(nfp));
        if(itens[i].nfps == NULL){
            printf("erro no malloc dos nfps do item %d", i);
            return -1; //caso malloc vazio, retorna -1
        }

        for(j = 0; j<numitens; j++){
            sprintf(nome, "NFP_%d_0_%d_0.txt", j+1, i+1); //i eh o NFP girando em torno de j
            //note que, como os itens sao dados a partir de 1 e os dados contamos de 0 em diante, vou olhar o nfp com os numeros em 1 a mais
            lenfp = fopen(nome, "r");

            do{
                fscanf(lenfp, "%c", &lixo);
            }while(lixo != '#');

            fscanf(lenfp, "%d", &itens[i].nfps[j].partes);

            itens[i].nfps[j].polig = (pol*)malloc((itens[i].nfps[j].partes*sizeof(pol)));
            if(itens[i].nfps[j].polig == NULL){
                printf("erro no malloc do vertices das partes do nfp %d do item %d",j, i);
                return -1;
            }
            
            for(k = 0; k < itens[i].nfps[j].partes; k++){
                do{
                    fscanf(lenfp, "%c", &lixo);
                }while(lixo != '#');

                fscanf(lenfp, "%d", &itens[i].nfps[j].polig[k].n);
                itens[i].nfps[j].polig[k].ponto = (vert*)malloc((itens[i].nfps[j].polig[k].n)*sizeof(vert));
                if(itens[i].nfps[j].polig[k].ponto == NULL){
                    printf("erro no malloc do vertices da parte %d do nfp %d do item %d", k, j, i);
                    return -1;
                }

                do{
                    fscanf(lenfp, "%c", &lixo);
                }while(lixo != '#');

                for(m = 0; m<(itens[i].nfps[j].polig[k].n); m++){ //le as coordenadas de cada vertice do NFP
                    fscanf(lenfp, "%lf %lf", &itens[i].nfps[j].polig[k].ponto[m].x, &itens[i].nfps[j].polig[k].ponto[m].y);
                }
            }    
            fclose(lenfp);
        }    
    }
    return 0;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao que abre arquivos para ler as informacoes gerais (altura, numero de itens, 
demanda de cada item) e aloca-las no vetor de itens*/
item *LerDados(int *numitens, double *altura){
    char lixo, nome[21];
    int i, j, k, numvertices, checarNFP, cortes, defeitos, centro[2];
    int it; //utilizar outro contador quando estiver checando qual item estou vendo no arquivo, para caso esteja fora de ordem
    FILE *infos, *listait;
    double largura;

    infos = fopen("general_data.txt", "r");
        do{
            fscanf(infos, "%c", &lixo);
        }while(lixo != '#');

        fscanf(infos, "%lf", &largura);

        do{
            fscanf(infos, "%c", &lixo);
        }while(lixo != '#');

        fscanf(infos, "%lf", altura);

        do{
            fscanf(infos, "%c", &lixo);
        }while(lixo != '#');

        fscanf(infos, "%d", &cortes);

        do{
            fscanf(infos, "%c", &lixo);
        }while(lixo != '#');

        fscanf(infos, "%d", &defeitos);

        do{
            fscanf(infos, "%c", &lixo);
        }while(lixo != '#');

        fscanf(infos, "%d", numitens);

        item *itens;
        itens = (item*) malloc((*numitens)*sizeof(item));
        if(itens == NULL){
            printf("erro do malloc dos itens");
            return NULL;
        }
        
        for(i = 0; i<*numitens; i++){
            itens[i].copias = 0;
            do{
                fscanf(infos, "%c", &lixo);
            }while(lixo != '#');

            fscanf(infos, "%d", &it); //variavel de qual eh o item - caso fora de ordem
            it--; // -------> Como os itens comecam de 1 e o vetor comeca de 0, precisa subtrair 1 // os itens estarao no espaco 0 ate m-1
            do{
                fscanf(infos, "%c", &lixo);
            }while(lixo != '#');

            fscanf(infos, "%d", &itens[it].copias);

            do{
                fscanf(infos, "%c", &lixo);
            }while(lixo != '#');

            fscanf(infos, "%d", &itens[i].num_rotacoes);

            itens[i].rotacoes = (double*)malloc(itens[i].num_rotacoes*sizeof(double));
            if(itens[i].rotacoes==NULL){
                printf("erro no malloc das rotacoes do item %d", i);
                return NULL;
            }

            for(j = 0; j<itens[i].num_rotacoes; j++){
                fscanf(infos, "%lf", &itens[i].rotacoes[j]);
            }
            
            do{
                fscanf(infos, "%c", &lixo);
            }while(lixo != '#');
        }
    fclose(infos);

    for(i = 0; i<*numitens; i++){
        sprintf(nome, "piece%d.txt", i+1); //note que, como os itens sao dados a partir de 1 e os dados contamos de 0 em diante, vou olhar o nfp com os numeros em 1 a mais
        listait = fopen(nome, "r");
            do{
                fscanf(listait, "%c", &lixo);
            }while(lixo != '#');

            fscanf(listait, "%d %d", &centro[0], &centro[1]);

            do{
                fscanf(listait, "%c", &lixo);
            }while(lixo != '#');
            
            fscanf(listait, "%d", &itens[i].partes);//li o numero de partes do item
            itens[i].polig = (pol*)malloc((itens[i].partes)*sizeof(pol));
            if(itens[i].polig == NULL){
                printf("erro no malloc das partes do item %d", it);
                return NULL; //caso malloc vazio, retorna -1
            }
            
            for(j = 0; j < itens[i].partes; j++){
                do{
                    fscanf(listait, "%c", &lixo);
                }while(lixo != '#');
                
                
                fscanf(listait, "%d", &itens[i].polig[j].n);
                numvertices = itens[i].polig[j].n;
                

                itens[i].polig[j].ponto = (vert*) malloc(numvertices*sizeof(vert));
                if(itens[i].polig[j].ponto == NULL){
                    printf("erro do malloc dos vertices do item %d", it);
                    return NULL;
                }
                

                do{
                    fscanf(listait, "%c", &lixo);
                }while(lixo != '#');
                    
                for(k = 0; k<numvertices; k++){
                    fscanf(listait, "%lf %lf", &itens[i].polig[j].ponto[k].x, &itens[i].polig[j].ponto[k].y);
                }
            }
        fclose(listait);
    }

    checarNFP = LerNFPs(itens, *numitens);
    if(checarNFP != 0){
        printf("algo deu errado checando nfps...\n");
    }

    return itens;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para calculo da area de um poligono, a partir da formula de determinates com vertices, que precisam estar ordenados. 
se os vertices estiverem no sentido horario, o resultado eh o oposto, entao realizamos essa correcao*/
double Area(item it_atual){
    double area;
    int i, j;
    area = 0;

    for(i = 0; i<it_atual.partes; i++){
        for(j = 0; j<(it_atual.polig[i].n); j++){
            area += (it_atual.polig[i].ponto[j].x*it_atual.polig[i].ponto[(j+1)%it_atual.polig[i].n].y) - (it_atual.polig[i].ponto[(j+1)%it_atual.polig[i].n].x*it_atual.polig[i].ponto[j].y);
        }
    }
    area = area/2;

    if(area < 0 ){
        return -area;
    }
    return area;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para ordenar a area dos itens de forma descrescente, retornando apenas os indices 
ordenados num vetor, de forma que nao embaralhe as areas dos itens sem querer*/
int* OrdenarAreaDecrescente(item *itens, int numtotal, int numitens){
    int i, j, k, auxlocal;

    int *indice;
    indice = (int*) malloc((numtotal)*sizeof(int));
    if(indice == NULL){
        printf("erro do malloc do indice");
        return NULL;
    }

    k = 0;
    for(i = 0; i<numitens; i++){
        for(j=0; j<itens[i].copias; j++){
            indice[k] = i; //cada espaco do indice recebe o numero do item associado a ele (tem o mesmo numero de copias de um item com seu valor no vetor)
            k++;
        }
    }

    for(int i = 0; i<(numtotal-1); i++){
        for(int j = i+1; j<numtotal; j++){
          if(itens[indice[j]].area>itens[indice[i]].area){
            auxlocal = indice[i];
            indice[i] = indice[j];
            indice[j] = auxlocal;
          }
        }
      }
    return indice;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para calcular o ponto mais a esquerda de um poligono*/
double Esquerda(item it_atual){
    double esquerda;
    int i, j;

    esquerda = it_atual.polig[0].ponto[0].x;
    for(i = 0; i<it_atual.partes; i++){
        for(j = 0; j<it_atual.polig[i].n; j++){
            if(it_atual.polig[i].ponto[j].x < esquerda){
                esquerda = it_atual.polig[i].ponto[j].x;
            }
        }
    }
    return esquerda;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para calcular o ponto mais a direita de um poligono*/
double Direita(item it_atual){
    double direita;
    int i, j;

    direita = it_atual.polig[0].ponto[0].x;
    for(i = 0; i<it_atual.partes; i++){
        for(j = 0; j<it_atual.polig[i].n; j++){
            if(it_atual.polig[i].ponto[j].x > direita){
                direita = it_atual.polig[i].ponto[j].x;
            }
        }
    }
    return direita;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para calcular o ponto mais acima de um poligono*/
double Alto(item it_atual){
    double alto;
    int i, j;

    alto = it_atual.polig[0].ponto[0].y;
    for(i = 0; i<it_atual.partes; i++){
        for(j = 0; j<it_atual.polig[i].n; j++){
            if(it_atual.polig[i].ponto[j].y > alto){
                alto = it_atual.polig[i].ponto[j].y;
            }
        }
    }
    return alto;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para calcular o ponto mais abaixo de um poligono*/
double Baixo(item it_atual){
    double baixo;
    int i, j;

    baixo = it_atual.polig[0].ponto[0].y;
    for(i = 0; i<it_atual.partes; i++){
        for(j = 0; j<it_atual.polig[i].n; j++){
            if(it_atual.polig[i].ponto[j].y < baixo){
                baixo = it_atual.polig[i].ponto[j].y;
            }
        }
    }
    return baixo;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para calcular a maior largura possivel do recipiente*/
double MuitoGrande(item *itens, int numitens){
    double mg = 0;
    int i;

    for(i = 0; i<numitens; i++){
        mg += (itens[i].direita - itens[i].esquerda)*itens[i].copias;
    }
    return mg+1;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para calcular a interseccao entre dois segmentos.
recebe 4 verts p - lista dos segmentos e o espaco de 2 vert atual - que sera alterado com os valores recebidos*/
int Intersec(vert* p, vert*atual){

    double num1, num2, t, u, denom, maiord, dist;
    int i, j, pontoa, pontod;
    
    denom = ((p[3].x-p[2].x)*(p[1].y-p[0].y))-((p[3].y-p[2].y)*(p[1].x-p[0].x));
    num1 = ((p[3].x-p[2].x)*(p[2].y-p[0].y))-((p[3].y-p[2].y)*(p[2].x-p[0].x));
    num2 = ((p[1].x-p[0].x)*(p[2].y-p[0].y))-((p[1].y-p[0].y)*(p[2].x-p[0].x));

    t = num1/denom;
    u = num2/denom;

    maiord = -1;

    if((fabs(denom)) <= 1.0e-5 && fabs(num1) <= 1.0e-5){ //caso onde os 4 pontos sao colineares
        if((((p[0].x<=p[2].x)&&(p[2].x<=p[1].x))&&((p[0].y<=p[2].y)&&(p[2].y<=p[1].y)))||(((p[0].x<=p[3].x)&&(p[3].x<=p[1].x))&&((p[0].y<=p[3].y)&&(p[3].y<=p[1].y)))||(((p[1].x<=p[2].x)&&(p[2].x<=p[0].x))&&((p[1].y<=p[2].y)&&(p[2].y<=p[0].y)))||(((p[1].x<=p[3].x)&&(p[3].x<=p[0].x))&&((p[1].y<=p[3].y)&&(p[3].y<=p[0].y)))){//checar se se interceptam
            for(i=0; i<3; i++){
                for(j=1; j<4; j++){
                    dist = sqrt(sq(p[j].x-p[i].x) + sq(p[j].y-p[i].y));
                    if(dist>maiord){
                        maiord = dist;
                        pontoa = i;
                        pontod = j;
                    }
                }
            }
            for(i=0; i<4; i++){
                if((i!=pontoa)&&(i!=pontod)){
                    atual[0] = p[i];
                    for(j = i+1; j<4; j++){
                        if((j!=pontoa)&&(j!=pontod)){
                            atual[1] = p[j];
                            return 2; //calculou 2 interseccoes
                        }
                    }
                }
            }
        }
    }else if (denom != 0){ //caso os segmentos nao sejam paralelos
        if((-1.0e-5<=t)&&(t<=(1+1.0e-5))&&(-1.0e-5<=u)&&(u<=(1+1.0e-5))){
            atual[0].x = p[0].x + t*(p[1].x - p[0].x);
            atual[0].y = p[0].y + t*(p[1].y - p[0].y);
            return 1; //calculou 1 interseccao
        }
    }
     return 0; //calculou 0 interseccoes
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao que checa se o vert p ta fora do pol nfp deslocado (para a posicao final do item)*/
int Fora(vert p, nfp nfps, vert desloca){
    int i, j, k, fora = 0;
    double valor;
    nfp nfp_desloc;

    nfp_desloc.polig = (pol*)malloc(nfps.partes*sizeof(pol));
    if(nfp_desloc.polig == NULL){
        printf("deu erro na criação das partes do nfp deslocado");
        return -1;
    }

    for(i = 0; i < nfps.partes; i++){
        nfp_desloc.polig[i].ponto = (vert*)malloc((nfps.polig[i].n)*sizeof(vert));
        if(nfp_desloc.polig[i].ponto == NULL){
            printf("deu erro na parte %d do nfp deslocado", i);
            return -1;
        }
    }

    for(i = 0; i<nfps.partes; i++){
        for(j = 0; j<nfps.polig[i].n; j++){
            nfp_desloc.polig[i].ponto[j].x = nfps.polig[i].ponto[j].x + desloca.x;
            nfp_desloc.polig[i].ponto[j].y = nfps.polig[i].ponto[j].y + desloca.y;
        }
    }

    for(i = 0; i<nfps.partes; i++){
        for(j = 0; j<nfps.polig[i].n; j++){
            valor = ((nfp_desloc.polig[i].ponto[j].x-nfp_desloc.polig[i].ponto[(j+1)%nfps.polig[i].n].x)*(nfp_desloc.polig[i].ponto[j].y-p.y))-((nfp_desloc.polig[i].ponto[j].y-nfp_desloc.polig[i].ponto[(j+1)%nfps.polig[i].n].y)*(nfp_desloc.polig[i].ponto[j].x-p.x));
            if(valor <=1.0e-5){
                fora++;
                break; ////////////////////////////DUVIDA
            }
        }
    }


    if(fora>=nfps.partes){
        for(k = 0; k < nfps.partes; k++){
            free(nfp_desloc.polig[k].ponto);
        }
        free(nfp_desloc.polig);
        return 1; //verdadeiro - esta fora
    }
    
    for(k = 0; k < nfps.partes; k++){
        free(nfp_desloc.polig);
        free(nfp_desloc.polig[k].ponto);
    }
    return 0; //falso - esta dentro
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao que verifica se um ponto eh valido para alocar um item, checando se esta fora do ifp ou dentro de algum nfp*/
int Valido(vert pt_atual, int it_atual, item *itens, int numitens, double altura){
    int i, j;

    if(it_atual == 3 && itens[it_atual].cop_alocad == 2 && pt_atual.x > 23 && pt_atual.x < 23.3 && pt_atual.y > 4.6 && pt_atual.y < 5){
                printf("%.2lf %.2lf\n", pt_atual.x, pt_atual.y);
            }
    
    if((pt_atual.x<(-itens[it_atual].esquerda))||(pt_atual.y>(altura - itens[it_atual].alto))||(pt_atual.y<(-itens[it_atual].baixo))){
        return 0; //invalido - fora do ifp
    }

    for(i=0; i<numitens; i++){
        for(j = 0; j<itens[i].cop_alocad; j++){
            /*if(it_atual == 3 && itens[it_atual].cop_alocad == 2 && pt_atual.x == 23.15 && pt_atual.y == 4.8){
                printf(" %.2lf %.2lf  ", pt_atual.x, pt_atual.y);
            }*/
            if(Fora(pt_atual, itens[it_atual].nfps[i], itens[i].pfinal[j]) == 0){
                return 0; //invalido - dentro de algum nfp
            }
        }
    }
    return 1; //valido
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao que aloca uma copia de um item, analisando cada interseccao possivel e, caso um ponto seja melhor que
o candidato anterior, vendo se eh um ponto valido para alocacao. caso seja, vira o novo candidato*/
vert AlocaUm(int it_atual, item *itens, int numitens, double mtgd, double altura){
    int i, j, k, m, n, n_int, f, g, h, alfa, beta;
    vert candidato, ifp[4];
    candidato.x = mtgd;
    candidato.y = 0;

///checando o ponto inferior esquerdo do ifp primeiro - melhor alocacao possivel
    ifp[1].x = -itens[it_atual].esquerda;
    ifp[1].y = -itens[it_atual].baixo;
    if(Valido(ifp[1], it_atual, itens, numitens, altura) == 1){
        return ifp[1];
    }

///checando as interseccoes entre ifps e nfps
    ifp[0].x = mtgd;
    ifp[0].y = -itens[it_atual].baixo;
    
    ifp[2].x = -itens[it_atual].esquerda;;
    ifp[2].y = altura - itens[it_atual].alto;
    
    ifp[3].x = mtgd;
    ifp[3].y = altura - itens[it_atual].alto;

    vert* p = malloc(4*sizeof(vert));
    if(p == NULL){
        printf("erro no malloc dos quatro pontos p - utilizados para verificar as interseccoes entre dois segmentos");
        return candidato;
    }

    vert* atual = malloc(2*sizeof(vert));
    if(atual == NULL){
        printf("erro no malloc dos vertices atuais");
        return candidato;
    }

    for(i = 0; i<numitens; i++){ //para cada item i
        for(j = 0; j<itens[i].cop_alocad; j++){ // ve quantas copias foram alocadas
            for(alfa = 0; alfa < itens[it_atual].nfps[i].partes; alfa++){ //alfa olha cada parte do nfp
                for(k = 0; k<itens[it_atual].nfps[i].polig[alfa].n; k++){ // k passa por cada aresta da parte alfa do nfp
                    p[0].x = itens[it_atual].nfps[i].polig[alfa].ponto[k].x + itens[i].pfinal[j].x; //os dois primeiros pontos sao nfps deslocados para a posicao da copia j
                    p[0].y = itens[it_atual].nfps[i].polig[alfa].ponto[k].y + itens[i].pfinal[j].y;
                    p[1].x = itens[it_atual].nfps[i].polig[alfa].ponto[(k+1)%itens[it_atual].nfps[i].polig[alfa].n].x + itens[i].pfinal[j].x;
                    p[1].y = itens[it_atual].nfps[i].polig[alfa].ponto[(k+1)%itens[it_atual].nfps[i].polig[alfa].n].y + itens[i].pfinal[j].y;
                    for(m = 0; m<3; m++){ //os outros dois pontos sao as arestas do ifp - tirando a aresta vertical da direita (infinita)
                        p[2] = ifp[m]; 
                        p[3] = ifp[m+1];
                        n_int = Intersec(p, atual); //checa as interseccoes e devolve quantas sao
                        for(n = 0; n<n_int; n++){ //passa por todas as intersec pra ver se sao validas
                            if((atual[n].x <candidato.x)||(atual[n].x == candidato.x && atual[n].y<candidato.y)){
                                if(Valido(atual[n], it_atual, itens, numitens, altura) == 1){
                                    candidato = atual[n];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
///checando as interseccoes entre nfps
    for(i = 0; i<numitens; i++){ //para cada item
        for(j = 0; j<itens[i].cop_alocad; j++){ //ve quantas copias foram alocadas
            for(alfa = 0; alfa<itens[it_atual].nfps[i].partes; alfa++){ //alfa olha cada parte do nfp
                for(k = 0; k<itens[it_atual].nfps[i].polig[alfa].n; k++){ // k passa por cada aresta da parte alfa do nfp
                    p[0].x = itens[it_atual].nfps[i].polig[alfa].ponto[k].x + itens[i].pfinal[j].x; //os dois primeiros pontos sao nfps deslocados para a posicao da copia j
                    p[0].y = itens[it_atual].nfps[i].polig[alfa].ponto[k].y + itens[i].pfinal[j].y;
                    p[1].x = itens[it_atual].nfps[i].polig[alfa].ponto[(k+1)%itens[it_atual].nfps[i].polig[alfa].n].x + itens[i].pfinal[j].x;
                    p[1].y = itens[it_atual].nfps[i].polig[alfa].ponto[(k+1)%itens[it_atual].nfps[i].polig[alfa].n].y + itens[i].pfinal[j].y;
                    for(g = j; g<itens[i].cop_alocad; g++){ //checa contra todas as >outras< copias ja alocadas desse mesmo item
                        for(beta = 0; beta<itens[it_atual].nfps[i].partes; beta++){ //passa por todas as partes de cada nfp
                            for(h = 0; h<itens[it_atual].nfps[i].polig[beta].n; h++){ 
                                p[2].x = itens[it_atual].nfps[i].polig[beta].ponto[h].x + itens[i].pfinal[g].x;
                                p[2].y = itens[it_atual].nfps[i].polig[beta].ponto[h].y + itens[i].pfinal[g].y;
                                p[3].x = itens[it_atual].nfps[i].polig[beta].ponto[(h+1)%itens[it_atual].nfps[i].polig[beta].n].x + itens[i].pfinal[g].x;
                                p[3].y = itens[it_atual].nfps[i].polig[beta].ponto[(h+1)%itens[it_atual].nfps[i].polig[beta].n].y + itens[i].pfinal[g].y;
                                n_int = Intersec(p, atual); //checa as interseccoes e devolve quantas sao
                                for(n = 0; n<n_int; n++){ //passa por todas as intersec pra ver se sao validas
                                    if((atual[n].x <candidato.x)||(atual[n].x == candidato.x && atual[n].y<candidato.y)){
                                        if(Valido(atual[n], it_atual, itens, numitens, altura) == 1){
                                            candidato = atual[n];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    for(f = i+1; f<numitens; f++){ //os outros dois pontos sao dos nfps dos >itens seguintes<
                        for(g = j; g<itens[f].cop_alocad; g++){ //checa contra >todas< as copias ja alocadas do item
                            for(beta = 0; beta < itens[it_atual].nfps[f].partes; beta++){
                                for(h = 0; h<itens[it_atual].nfps[f].polig[beta].n; h++){ 
                                    p[2].x = itens[it_atual].nfps[f].polig[beta].ponto[h].x + itens[f].pfinal[g].x;
                                    p[2].y = itens[it_atual].nfps[f].polig[beta].ponto[h].y + itens[f].pfinal[g].y;
                                    p[3].x = itens[it_atual].nfps[f].polig[beta].ponto[(h+1)%itens[it_atual].nfps[f].polig[beta].n].x + itens[f].pfinal[g].x;
                                    p[3].y = itens[it_atual].nfps[f].polig[beta].ponto[(h+1)%itens[it_atual].nfps[f].polig[beta].n].y + itens[f].pfinal[g].y;
                                    n_int = Intersec(p, atual); //checa as interseccoes e devolve quantas sao
                                    for(n = 0; n<n_int; n++){ //passa por todas as intersec pra ver se sao validas
                                        if((atual[n].x <candidato.x)||(atual[n].x == candidato.x && atual[n].y<candidato.y)){
                                            if(Valido(atual[n], it_atual, itens, numitens, altura) == 1){
                                                candidato = atual[n];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    free(atual);
    free(p);

    return candidato;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao que calcula o numero total de itens, o maior valor possivel do comprimento (y), ordena o indice e aloca todos os itens*/
int AlocaTodos(item *itens, int numitens, double altura){
    int i, j, numtotal;
    int *indice;
    double mtgd;
    vert teste;
    
    numtotal = 0;
    for(i = 0; i<numitens; i++){
        numtotal += itens[i].copias;
    }

    mtgd = MuitoGrande(itens, numitens);

    indice = OrdenarAreaDecrescente(itens, numtotal, numitens);

    
    for(i = 0; i<numtotal; i++){
        itens[indice[i]].pfinal[itens[indice[i]].cop_alocad] = AlocaUm(indice[i], itens, numitens, mtgd, altura);
        itens[indice[i]].cop_alocad++;
    }

    free(indice);

    return 0;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao que imprime as posicoes finais da alocacao em forma de leiaute com linguagem TikZ, para LateX*/
int Imprimir(item *itens, int numitens, double altura){
    int i, j, k, alfa, aux = 0;
    double comp_total = 0;
    int colors[23];

    for(i = 0; i<23 ; i++){
        colors[i] = 65+i;
    }

    for(i = 0; i<numitens; i++){
        for(j = 0; j<itens[i].copias; j++){
            for(alfa = 0; alfa < itens[i].partes; alfa++){
                //printf("\n %% item %d copia %d \n\\filldraw [color=black, fill=red!60] ", i, j);
                printf("\n %% item %d copia %d \n\\filldraw [fill= \\cor%c] ", i, j, (char)(colors[(aux)%(23)]));
                for(k = 0; k<itens[i].polig[alfa].n; k++){
                    printf(" (%.3lf,%.3lf) --", (itens[i].polig[alfa].ponto[k].x + itens[i].pfinal[j].x), (itens[i].polig[alfa].ponto[k].y + itens[i].pfinal[j].y));
                }
                if(comp_total < (itens[i].direita + itens[i].pfinal[j].x)){
                    comp_total = itens[i].direita + itens[i].pfinal[j].x;
                }
                printf(" cycle;\n");
            }
            aux++;
        }
    }
    printf("\n\\draw[blue, very thick] (%.2lf, 0) -- (0,0) -- (0,%.2lf) -- (%.2lf, %.2lf);\n", 1.2*comp_total, altura, 1.2*comp_total, altura);
    printf("\n \\draw[dotted, very thick] (%.2lf, 0) -- (%.2lf,%.2lf);\n", comp_total, comp_total, altura);
    return 0;
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao para limpar todos os mallocs feitos ao longo do programa*/
void Limpa(item *itens, int numitens){
    int i, j, k;
    for(i = 0; i<numitens; i++){
        for(j = 0; j<numitens; j++){
            for(k = 0; k<itens[i].nfps[j].partes; k++){
                free(itens[i].nfps[j].polig[k].ponto);
            }
            free(itens[i].nfps[j].polig);
        }

        for(j = 0; j<itens[i].partes; j++){
            free(itens[i].polig[j].ponto);
        }
        free(itens[i].polig);
        free(itens[i].rotacoes);
        free(itens[i].nfps);
        free(itens[i].pfinal);
    }
    free(itens);
}
//=============================================================================================================================================
//=============================================================================================================================================
/* funcao principal, chama LerDados e, para cada item, chama os calculos da area, pontos extremos,
coloca o numero de copias ja alocadas em 0 e chama AlocaTodos e Imprimir*/
int main() {
    int numitens; //numero de ITENS (ITEM = poligono que queremos alocar)
    int i, teste, teste2;
    double altura; // altura do recipiente

    item *itens;
    itens = LerDados(&numitens,&altura);

    for(i = 0; i<numitens; i++){
        itens[i].area = Area(itens[i]);
        itens[i].esquerda = Esquerda(itens[i]);
        itens[i].alto = Alto(itens[i]);
        itens[i].baixo = Baixo(itens[i]);
        itens[i].direita = Direita(itens[i]);
        itens[i].cop_alocad = 0;
        itens[i].pfinal = (vert*)malloc((itens[i].copias)*sizeof(vert));
        if(itens[i].pfinal == NULL){
            printf("erro do malloc do vetor de pontos finais do item %d\n", i);
            return 0;
        } 
    }


    teste = AlocaTodos(itens, numitens, altura);
    teste2 = Imprimir(itens, numitens, altura);

    Limpa(itens, numitens);
    
    return 0;
}
