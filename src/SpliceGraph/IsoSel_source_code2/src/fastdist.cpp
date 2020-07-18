#include "fastdist.h"
/*;;;;;;;;;;;;;;;;;;;;;;;;;; Compilation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; No OpenMP library:                                                        ;
; g++ -I eigen -O -o fastdist fastdist.cpp                                  ;
;                                                                           ;
; If OpenMP:                                                                ;
; g++ -I eigen -O -o fastdist fastdist.cpp -fopenmp                         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void   Initialize(float **, FILE *, int, POINTERS *);
void   Compute_sums_Sx(float **, int);
void   Best_pair(float **, int, int *, int *, int);
void   Finish(float **, int, POINTERS *, FILE *);
void   Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post);
void   write_trees(int, POINTERS *, FILE *);
void   bootstrap_weights(allseq *, unsigned);
void   exit_error(char *);
void   *my_calloc(int, size_t);
void   check_symmetry(const MatrixXd&);
void   set_exchanges(int, MatrixXd&);
void   set_equilibrium(int, int, SEA_VIEW, VectorXd&);
void   write_dst(FILE *, matrix *);
void   test_duplicated_names(allseq *);
void   dist_miss_ultra(matrix *, double);
void   dist_miss_add(matrix *, double);
void   Init_Mat(matrix *, allseq *);
void   majuscules(char *);
void   Free(void *);
void   Free_Mat(matrix *);
void   compute_dst(FILE *, allseq *, SEA_VIEW, int, int, int, double, int, unsigned);
void   compute_trees(FILE *, FILE *);
void   root_trees(FILE *, FILE *);
void   split_trees(FILE *, char *, int);

int    Emptied(int, float **);
int    Symmetrize(float **, int);
int    read_fasta_align(FILE *, char ***, char ***,  char ***, char **,
       char **, int);
int    Is_Ambigu(char *, int, int);
int    my_strncmp(const char*, const char*, size_t);
int    calc_tree_count(char *);

float  Distance(int, int, float **);
float  Variance(int, int, float **);
float  Sum_S(int, float **);
float  Agglomerative_criterion(int, int, float **, int);
float  Branch_length(int, int, float **, int);
float  Reduction4(int, float, int, float, int, float, float **);
float  Reduction10(int, int, int, float, float, float **);
float  Lamda(int, int, float, float **, int, int);
float  Finish_branch_length(int, int, int, float **);

double next_random(unsigned);
double calc_dist_model(double, const VectorXd&, const VectorXd&,
    const MatrixXd&, const MatrixXd&);
       
char *get_name(char *);
char *nextpar(char *);
char **strsplit(const char *, const char *, size_t *);

matrix *Make_Mat(int);
matrix *JC69_Dist(allseq *, model *);
matrix *Obs_Dist(allseq *, model *);
matrix *Obs_Dist_Gaps(allseq *, model *);
matrix *Kimura_p_Dist(allseq *);
matrix *Gamma_p_Dist(allseq *, double);
matrix *calc_dist_matrix(allseq *, int, int, double, char **);
allseq *view_to_allseq(SEA_VIEW *, int);

/*;;;;;;;;;;;;;;;;;;;;;;;;;; Initialize ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function reads an input file and return the            ;
;               delta matrix and trees: the list of taxa.                   ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              FILE *input    : pointer to input file                       ;
;              int n          : number of taxa                              ;
;              char **trees   : list of taxa                                ;
;                                                                           ;
; return value:                                                             ;
;              float **delta : delta matrix                                 ;
;              char *trees    : list of taxa                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Initialize(float **delta, FILE *input, int n, POINTERS *trees)
{
    int lig;                                          /* matrix line       */
    int col;                                          /* matrix column     */
    float distance;
    char name_taxon[LEN];                             /* taxon’s name      */
    WORD *name;
  
    for (lig = 1; lig <= n; lig++)
    {
        fscanf(input, "%s", name_taxon);               /* read taxon’s name */
        name = (WORD *)my_calloc(1, sizeof(WORD));     /* taxon’s name is   */
        strcpy(name->name, name_taxon);
        name->suiv = NULL;
        trees[lig].head = name;
        trees[lig].tail = name;
        for (col = 1; col <= n; col++)
        {
            fscanf(input, "%f", &distance);      /* read the distance  */
            delta[lig][col] = distance;
        }
    }
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; write_trees;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function prints out the tree in the output file.       ;
;                                                                           ;
; input       :                                                             ;
;              POINTERS *trees : pointer to the subtrees.                   ;
;              int i          : indicate the subtree i to be printed.       ;
:              FILE *output   : pointer to the output file.                 ;
;                                                                           ;
; return value: The phylogenetic tree in the output file.                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void write_trees(int i, POINTERS *trees, FILE *output)
{
    WORD *parcour;
    
    parcour = trees[i].head;
    while (parcour != NULL)
    {
        fprintf(output, "%s", parcour->name);
        parcour = parcour->suiv;
    }
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Symmetrize  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifies if the delta matrix is symmetric;    ;
;               if not the matrix is made symmetric.                        ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              int n          : number of taxa                              ;
;                                                                           ;
; return value:                                                             ;
;              int symmetric  : indicate if the matrix has been made        ;
;                               symmetric or not                            ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Symmetrize(float **delta, int n)
{
    int lig;                                 /* matrix line        */
    int col;                                 /* matrix column      */
    float value;                             /* symmetrized value  */
    int symmetric;
  
    symmetric = 1;
    for (lig = 1; lig  <=  n; lig++)
    {
        for (col = 1; col < lig; col++)
        {
            if (delta[lig][col] != delta[col][lig])
            {
                value = (delta[lig][col] + delta[col][lig])/2;
                delta[lig][col] = value;
                delta[col][lig] = value;
                symmetric = 0;
            }
        }
    }
    return(symmetric);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Concatenate ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function concatenates a string to another.             ;
;                                                                           ;
; input       :                                                             ;
;      char *chain1    : the string to be concatenated.                     ;
;      int ind         : indicate the subtree to which concatenate the      ;
;                        string                                             ;
;      POINTERS *trees  : pointer to subtrees.                              ;
;      int post        : position to which concatenate (front (0) or        ;
;                        end (1))                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post)
{
    WORD *bran;
  
    bran = (WORD *)my_calloc(1, sizeof(WORD));
    strcpy(bran->name, chain1);
    bran->suiv = NULL;
    if (post == 0)
    {
        bran->suiv = trees[ind].head;
        trees[ind].head = bran;
    }
    else
    {
        trees[ind].tail->suiv = bran;
        trees[ind].tail = trees[ind].tail->suiv;
    }
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Distance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve ant return de distance between taxa  ;
;               i and j from the delta matrix.                              ;
;                                                                           ;
; input       :                                                             ;
;               int i          : taxon i                                    ;
;               int j          : taxon j                                    ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               float distance : dissimilarity between the two taxa         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Distance(int i, int j, float **delta)
{
    if (i > j)
        return(delta[i][j]);
    else
        return(delta[j][i]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Variance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve and return the variance of the       ;
;               distance between i and j, from the delta matrix.            ;
;                                                                           ;
; input       :                                                             ;
;               int i           : taxon i                                   ;
;               int j           : taxon j                                   ;
;               float **delta  : the delta matrix                           ;
;                                                                           ;
; return value:                                                             ;
;               float distance : the variance of  Dij                       ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Variance(int i, int j, float **delta)
{
    if (i > j)
        return(delta[j][i]);
    else
        return(delta[i][j]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Emptied ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifie if a line is emptied or not.          ;
;                                                                           ;
; input       :                                                             ;
;               int i          : subtree (or line) i                        ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               0              : if not emptied.                            ;
;               1              : if emptied.                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Emptied(int i, float **delta)      /* test if the ith line is emptied */
{
    return((int)delta[i][0]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Sum_S;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function retrieves the sum Sx from the diagonal       ;
;                of the delta matrix.                                       ;
;                                                                           ;
;  input       :                                                            ;
;               int i          : subtree i                                  ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
;  return value:                                                            ;
;                float delta[i][i] : sum Si                                 ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Sum_S(int i, float **delta)          /* get sum Si form the diagonal */
{
    return(delta[i][i]);
}

/*;;;;;;;;;;;;;;;;;;;;;;;Compute_sums_Sx;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function computes the sums Sx and store them in the    ;
;               diagonal the delta matrix.                                  ;
;                                                                           ;
; input       :                                                             ;
;                  float **delta : the delta matrix.                        ;
;                  int n          : the number of taxa                      ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Compute_sums_Sx(float **delta, int n)
{
    float sum;
    int i, j;
  
    for (i = 1; i <= n ; i++)
    {
        if (!Emptied(i, delta))
        {
          sum = 0;
          for (j = 1; j <= n; j++)
          {
              if (i != j && !Emptied(j, delta))     /* compute the sum Si */
                  sum = sum + Distance(i, j, delta);
          }
        }
        delta[i][i] = sum;   /* store the sum Si in delta’s diagonal */
    }
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Best_pair;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function finds the best pair to be agglomerated by    ;
;                minimizing the agglomerative criterion (1).                ;
;                                                                           ;
;  input       :                                                            ;
;                float **delta : the delta matrix                           ;
;                int r          : number of subtrees                        ;
;                int *a         : contain the first taxon of the pair       ;
;                int *b         : contain the second taxon of the pair      ;
;                int n          : number of taxa                            ;
;                                                                           ;
;  return value:                                                            ;
;                int *a         : the first taxon of the pair               ;
;                int *b         : the second taxon of the pair              ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Best_pair(float **delta, int r, int *a, int *b, int n)
{
    float Qxy;                     /* value of the criterion calculated */
    int x, y;                      /* the pair which is tested          */
    float Qmin;                    /* current minimun of the criterion  */
  
    Qmin = 1.0e300;
    for (x = 1; x <= n; x++)
    {
        if (!Emptied(x, delta))
        {
            for (y = 1; y < x; y++)
            {
                if (!Emptied(y, delta))
                {
                    Qxy = Agglomerative_criterion(x, y, delta, r);
                    if (Qxy < Qmin - 0.000001)
                    {
                        Qmin = Qxy;
                        *a = x;          
                        *b = y;
                    }
                }  
            }
        }
    }
}

/*;;;;;;;;;;;;;;;;;;;;;;Finish_branch_length;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description :  Compute the length of the branch attached                 ;
;                 to the subtree i, during the final step                   ;
;                                                                           ;
;  input       :                                                            ;
;                int i          : position of subtree i                     ;
;                int j          : position of subtree j                     ;
;                int k          : position of subtree k                     ;
;                float **delta :                                            ;
;                                                                           ;
;  return value:                                                            ;
;                float length  : The length of the branch                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Finish_branch_length(int i, int j, int k, float **delta)
{
    float length;
    length = 0.5*(Distance(i, j, delta) + Distance(i, k, delta)
          - Distance(j, k, delta));
    return(length);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Finish;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function compute the length of the lasts three        ;
;                subtrees and write the tree in the output file.            ;
;                                                                           ;
;  input       :                                                            ;
;                float **delta  : the delta matrix                          ;
;                int n           : the number of taxa                       ;
;                WORD *trees   : list of subtrees                           ;
;                                                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Finish(float **delta, int n, POINTERS *trees, FILE *output)
{
    int l = 1;
    int i = 0;
    float length;
    char *str;
    WORD *bidon;
    WORD *ele;
    int last[3];                       /* the last three subtrees     */
  
    str = (char *)my_calloc(LEN, sizeof(char));
    while (l <= n)
    { /* find the last tree subtree  */
        if (!Emptied(l, delta))
        {
            last[i] = l;
            i++;
        }
        l++;
    }
    length = Finish_branch_length(last[0], last[1], last[2], delta);
    fprintf(output,"(");
    write_trees(last[0], trees, output);
    fprintf(output, ":");
    fprintf(output,"%f,", length);
    length = Finish_branch_length(last[1], last[0], last[2], delta);
    write_trees(last[1], trees, output);
    fprintf(output,":");
    fprintf(output,"%f,", length); 
    length = Finish_branch_length(last[2], last[1], last[0], delta);
    write_trees(last[2], trees, output);
    fprintf(output,":");
    fprintf(output,"%f", length);
    fprintf(output,");");
    fprintf(output,"\n");
//     for (i = 0; i < 3; i++)
//     {
//         bidon = trees[last[i]].head;
//         ele = bidon;
//         while (bidon != NULL)
//         {
//             ele = ele->suiv;
//             free(bidon);
//             bidon = ele;
//         }
//     }
    free(str);
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                          Formulae                                         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Agglomerative_criterion(int i, int j, float **delta, int r)
{
    float Qij;
    
    Qij = (r - 2)*Distance(i, j, delta) - Sum_S(i, delta) - Sum_S(j, delta);
    return(Qij);                       
}

float Branch_length(int a, int b, float **delta, int r)
{
    float length;
  
    length = 0.5*(Distance(a, b, delta) + (Sum_S(a, delta)
        - Sum_S(b, delta))/(r - 2));
    return(length);
}

float Reduction4(int a, float la, int b, float lb, int i, float lamda,
         float **delta)
{
    float Dui;
    
    Dui = lamda*(Distance(a, i, delta) - la)
        + (1 - lamda)*(Distance(b, i, delta) - lb);
    return(Dui);
}

float Reduction10(int a, int b, int i, float lamda, float vab,
    float **delta)
{
    float Vci;
  
    Vci = lamda*Variance(a, i, delta) + (1 - lamda)*Variance(b, i, delta)
        - lamda*(1 - lamda)*vab;
    return(Vci);
}

float Lamda(int a, int b, float vab, float **delta, int n, int r)
{
    float lamda = 0.0;
    int i;
  
    if (vab==0.0)
        lamda=0.5;
    else
    {
        for (i = 1; i <= n ; i++)
        {
            if (a != i && b != i && !Emptied(i, delta))
                lamda =lamda + (Variance(b, i, delta) - Variance(a, i, delta));
        }
        lamda = 0.5 + lamda/(2*(r - 2)*vab);
    }
    if (lamda > 1.0)
        lamda = 1.0;
    if (lamda < 0.0)
        lamda = 0.0;
    return(lamda);
}

double next_random(unsigned seed)
{
    static int first = TRUE;
    
    if (first)
    {
        first = FALSE;
        #ifdef WIN32
        SYSTEMTIME timing;
        GetSystemTime(&timing);
        if (seed == 0)
            seed = (unsigned)(timing.wMilliseconds + 1000*(timing.wSecond +
            60*(timing.wMinute + 60*(timing.wHour + 24*timing.wDay))));
        srand(seed);
        #else
        struct tms timing;
        if (seed == 0)
            seed = (unsigned)times(&timing);
        srandom(seed);
        #endif
    }
    #ifdef WIN32
    return rand() / (RAND_MAX + (double)1.0);
    #else
    return random() / (RAND_MAX + (double)1.0);
    #endif
}

void bootstrap_weights(allseq *seqs, unsigned seed)
{
    phydbl buff;
    int j, position;
    
    memset(seqs->wght, 0, seqs->clean_len*sizeof(int));
    for (j = 0; j < seqs->clean_len; j++)
    {
        buff = next_random(seed);
        buff *= seqs->clean_len;
        position = (int)floor(buff);
        seqs->wght[position] += 1;
    }
}

int read_fasta_align(FILE *in, char ***pseq, char ***pseqname, 
    char ***pcomments, char **pheader, char **err_message, int spaces_in_names)
{
    int totseqs, lseq, l2, l, lenseqs;
    char line[300], *p, *i, c, *q, *r;
    static char ret_message[200];
    char **seq, **seqname, **comments, *tmpseq = NULL;

    *ret_message = 0;
    *err_message = ret_message;
/*
    Number of sequences in file
*/
    totseqs = 0;
    while (fgets(line, sizeof(line), in) != NULL)
    {
        if (*line == '>')
            totseqs++;
    }
    rewind(in);
    seq = (char **)malloc(totseqs * sizeof(char *));
    if (seq == NULL)
        goto nomem;
    comments = (char **)malloc(totseqs * sizeof(char *));
    if (comments == NULL)
        goto nomem;
    seqname = (char **)malloc(totseqs * sizeof(char *));
    if (seqname == NULL)
        goto nomem;
    *pseq = seq;
    *pcomments = comments;
    *pseqname = seqname;
    lenseqs = MAXLENSEQ;
    tmpseq = (char *)malloc(lenseqs + 1);
    if (tmpseq == NULL)
        goto nomem;
    totseqs = -1;
    i = fgets(line, sizeof(line), in);
    if (line[0] != '>')
    {
        strcpy(ret_message, "Error: file not in Fasta format");
        totseqs = -1;
        goto fini;
    }
    while (i != NULL)
    {
/*
        Finish reading very long title line
*/
        c = line[strlen(line) - 1];
        while (c != '\n' && c != '\r' && c != EOF)
            c = getc(in);
        q = line + strlen(line) - 1;
        while (q > line + 1 && (*q == '\n' || *q == '\r'))
            *(q--) = 0;
        totseqs++;
        p = line + 1; 
        while (*p == ' ')
            p++;
        if (spaces_in_names)
        {
          while (*p != '\n')
              p++;
          while (*(p-1) == ' ')
              p--;
        }
        else
        {
            while (*p != ' ' && *p != '\n')
                p++;
        }
        r = line + 1;
        while (*r == ' ')
            r++;
        l = p - r;
        if ((seqname[totseqs] = (char *)malloc(l+1)) == NULL)
            goto nomem;
        memcpy(seqname[totseqs], r, l);
        seqname[totseqs][l] = 0;
/*
        Use rest of title line, if any, as comment
*/
        while (*p == ' ')
            p++;
        l = q - p + 1;
        if (l > 0)
        {
            comments[totseqs] = (char *)malloc(l + 3);
            if (comments[totseqs] != NULL)
            {
                strcpy(comments[totseqs], ";");
                strcpy(comments[totseqs] + 1, p);
                strcpy(comments[totseqs] + l + 1, "\n");
            }
        }
        else
            comments[totseqs] = NULL;
        lseq = 0;
        while ((i = fgets(line, sizeof(line), in)) != NULL && *i != '>' )
        {
            l2 = strlen(line);
            if (line[l2 - 1] == '\n' )
                l2--;
            while (l2 > 0 && line[l2-1] == ' ')
                l2--;
            if (lseq + l2 > lenseqs)
            {
              lenseqs += MAXLENSEQ;
              tmpseq = (char *)realloc(tmpseq, lenseqs + 1);
              if (tmpseq == NULL)
                  goto nomem;
            }
/*
            Copy seq data excluding spaces (because of Gblocks)
*/
            p = tmpseq + lseq;
            q = line;
            while (q < line + l2)
            {
                if (*q != ' ')
                    *(p++) = *q;
                q++;
            }
            lseq += p - (tmpseq+lseq);
        }
        tmpseq[lseq] = '\0';
        seq[totseqs] = (char *)malloc(lseq + 1);
        if (seq[totseqs] == NULL)
            goto nomem;
        memcpy(seq[totseqs], tmpseq, lseq + 1);
    }
fini:
    fclose(in);
    if (tmpseq != NULL)
        free(tmpseq);
    *pheader = NULL;
    return totseqs + 1;
nomem:
    sprintf(ret_message, "Error: not enough memory");
    totseqs = -1;
    goto fini;
}

void exit_error(char *s)
{
    fprintf(stderr, "%s\n", s);
    fflush(NULL);
    exit(1);
}

void *my_calloc(int nb, size_t size)
{
    void *allocated;
    char err_message[] = "Error: not enough memory";
    
    if ((allocated = calloc((size_t)nb, (size_t)size)) != NULL)
        return allocated;
    else
        exit_error(err_message);
    return NULL;
}

matrix *Make_Mat(int n_otu)
{
    matrix *mat;
    int i;
    
    mat = (matrix *)my_calloc(1, sizeof(matrix));
    mat->n_otu = n_otu;
    mat->P = (phydbl **)my_calloc(n_otu, sizeof(phydbl *));
    mat->Q = (phydbl **)my_calloc(n_otu, sizeof(phydbl *));
    mat->dist = (phydbl **)my_calloc(n_otu, sizeof(phydbl *));
    mat->on_off = (int *)my_calloc(n_otu, sizeof(int));
    mat->name = (char **)my_calloc(n_otu, sizeof(char *));
    mat->tip_node = (node **)my_calloc(n_otu, sizeof(node *));  
    for (i = 0; i < n_otu; i++)
    {
        mat->P[i] = (phydbl *)my_calloc(n_otu, sizeof(phydbl));
        mat->Q[i] = (phydbl *)my_calloc(n_otu, sizeof(phydbl));
        mat->dist[i] = (phydbl *)my_calloc(n_otu, sizeof(phydbl));
    }
    return mat;
}

void Init_Mat(matrix *mat, allseq *data)
{
    int i;
    
    mat->n_otu = data->n_otu;
    mat->r = mat->n_otu;
    mat->curr_int = mat->n_otu;
    mat->method = 1;
    for(i = 0; i < data->n_otu; i++)
    {
        mat->name[i] = strdup(data->c_seq[i]->name);
        mat->on_off[i] = 1;
    }
}

int Is_Ambigu(char *state, int datatype, int stepsize)
{
    int i;
    
    if (datatype == NT)
    {
        for (i = 0; i < stepsize; i++)
        {
            if (strchr("MRWSYKBDHVNXO?-.", state[i]))
                return 1;
        }
    }
    else
    {
        if (strchr("X?-.", state[0]))
            return 1;
    }
    return 0;
}

int my_strncmp(const char* s1, const char* s2, size_t n)
{
    while(n--)
    {
        if (s1[n] != s2[n])
            return 1;
        else
        {
            if ((s1[n] != 'X') )
                return 2;
        }
    }
    return 0;
}

void Free(void *p)
{
    if (p != NULL)
        free(p);
}

void Free_Mat(matrix *mat)
{
    int i;
    
    for (i = 0; i < mat->n_otu; i++)
    {
        Free(mat->P[i]);
        Free(mat->Q[i]);
        Free(mat->dist[i]);
        Free(mat->name[i]);
    }
    Free(mat->P);
    Free(mat->Q);
    Free(mat->dist);
    Free(mat->name);
    Free(mat->tip_node);
    Free(mat->on_off);
    Free(mat);
}

matrix *JC69_Dist(allseq *data, model *mod)
{
    int site, i, j, k;
    phydbl x;
    matrix *mat;
    phydbl **len;
    
    len = (phydbl **)my_calloc(data->n_otu,sizeof(phydbl *));
    for (i = 0; i < data->n_otu; i++)
        len[i] = (phydbl *)my_calloc(data->n_otu, sizeof(phydbl));
    mat = Make_Mat(data->n_otu);
    Init_Mat(mat, data);
    Fors (site, data->c_seq[0]->len, mod->stepsize)
    {
        for (j = 0; j < data->n_otu - 1; j++)
        {
            for (k = j + 1; k < data->n_otu; k++)
            {
                if ((!Is_Ambigu(data->c_seq[j]->state+site, mod->datatype,
                mod->stepsize)) && (!Is_Ambigu(data->c_seq[k]->state+site,
                mod->datatype,mod->stepsize)))
                {
                    len[j][k] += data->wght[site];
                    len[k][j] = len[j][k];
                    if (strncmp(data->c_seq[j]->state+site,
                    data->c_seq[k]->state+site, mod->stepsize))
                        mat->P[j][k] += data->wght[site];
                }
            }
        }
    }
    for (i = 0; i < data->n_otu - 1; i++)
    {
        for (j = i + 1; j < data->n_otu; j++)
        {
            if (len[i][j])
                mat->P[i][j] /= len[i][j];
            else
                mat->P[i][j] = 1.;
            mat->P[j][i] = mat->P[i][j];
            x = 1.0 - (mod->ns/(mod->ns - 1.0))*mat->P[i][j];
            if (x <= 0.0)
            {
                mat->dist[i][j] = NAN;
                mat->dist[j][i] = mat->dist[i][j];
            }
            else
            {
                mat->dist[i][j] = -((mod->ns-1.)/mod->ns)*(phydbl)log(x);
                mat->dist[j][i] = mat->dist[i][j];
            }
        }
    } 
out:
    for (i = 0; i < data->n_otu; i++)
        free(len[i]);
    free(len);
    return mat;
}

matrix *Obs_Dist(allseq *data, model *mod)
{
    int site, i, j, k;
    phydbl unc_len;
    matrix *mat;
    phydbl **len;
    
    len = (phydbl **)my_calloc(data->n_otu, sizeof(phydbl *));
    For (i, data->n_otu)
        len[i] = (phydbl *)my_calloc(data->n_otu, sizeof(phydbl));
    unc_len = .0;
    mat = Make_Mat(data->n_otu);
    Init_Mat(mat, data);
    Fors (site, data->c_seq[0]->len, mod->stepsize)
    {
        for (j = 0; j < data->n_otu - 1; j++)
        {
            for (k = j + 1; k < data->n_otu; k++)
            {
                if ((!Is_Ambigu(data->c_seq[j]->state+site, mod->datatype,
                mod->stepsize)) && (!Is_Ambigu(data->c_seq[k]->state + site,
                mod->datatype, mod->stepsize)))
                {
                    len[j][k] += data->wght[site];
                    len[k][j] = len[j][k];
                    if (strncmp(data->c_seq[j]->state+site,
                            data->c_seq[k]->state+site,
                            mod->stepsize))
                        mat->P[j][k] += data->wght[site];
                }
            }
        }
    }
    for (i = 0; i < data->n_otu - 1; i++)
    {
        for (j = i + 1; j < data->n_otu; j++)
        {
            if (len[i][j])
                mat->P[i][j] /= len[i][j];
            else
                mat->P[i][j] = 1.0;
            mat->P[j][i] = mat->P[i][j];
            mat->dist[i][j] = mat->P[i][j];
            mat->dist[j][i] = mat->dist[i][j];
        }
    }
    for (i = 0;  i < data->n_otu; i++)
        free(len[i]);
    free(len);
    return mat;
}

matrix *Obs_Dist_Gaps(allseq *data, model *mod)
{
    int site, i, j, k;
    int lengthJK = 0;
    phydbl unc_len;
    matrix *mat;
    phydbl **len;
    
    len = (phydbl **)my_calloc(data->n_otu, sizeof(phydbl *));
    for (i = 0; i < data->n_otu; i++)
        len[i] = (phydbl *)my_calloc(data->n_otu, sizeof(phydbl));
    unc_len = .0;
    mat = Make_Mat(data->n_otu);
    Init_Mat(mat, data);

    for (j = 0; j < data->n_otu - 1; j++)
    {
        for (k = j + 1; k < data->n_otu; k++)
        {
            lengthJK = 0;
            Fors (site, data->c_seq[0]->len, mod->stepsize)
            {
                if ((!Is_Ambigu(data->c_seq[j]->state+site, mod->datatype, mod->stepsize))
                    && (!Is_Ambigu(data->c_seq[k]->state + site, mod->datatype, mod->stepsize)))
                {
                    len[j][k] += data->wght[site];
                    len[k][j] = len[j][k];
                    if ((my_strncmp(data->c_seq[j]->state+site, data->c_seq[k]->state + site,
                        mod->stepsize)) != 0 ) // it is not 2 gaps (X-X)
                        lengthJK+=1;
                    if (my_strncmp(data->c_seq[j]->state+site, data->c_seq[k]->state + site,
                        mod->stepsize) == 1)
                        mat->P[j][k] += data->wght[site];
                }
            }
            if (len[j][k])
                mat-> P[j][k] /= lengthJK;
            else
                mat->P[j][k] = 1.0;
            mat->P[k][j] = mat->P[j][k];
            mat->dist[j][k] = mat->P[j][k];
            mat->dist[k][j] = mat->dist[j][k];
        }
    }
    for (i = 0; i < data->n_otu; i++)
        free(len[i]);
    free(len);
    return mat;
}

matrix *Kimura_p_Dist(allseq *data)
{
    int i, j;
    model jcmodel;
    phydbl x;
    
    jcmodel.stepsize = 1;
    jcmodel.datatype = 1;
    jcmodel.ns = 20;
    matrix *mat = Obs_Dist(data, &jcmodel);
    for (i = 0; i < data->n_otu; i++)
    {
        for (j = i + 1; j < data->n_otu; j++)
        {
            x = 1 - mat->dist[i][j] - 0.2 * mat->dist[i][j] * mat->dist[i][j];
            if (x <= 0)
                mat->dist[i][j] = NAN;
            else
                mat->dist[i][j] = -log(x);
            mat->dist[j][i] = mat->dist[i][j];
        }
    }
    return mat;
}

matrix *Gamma_p_Dist(allseq *data, double a)
{
    int i, j;
    model jcmodel;
    phydbl x;
    
    jcmodel.stepsize = 1;
    jcmodel.datatype = 1;
    jcmodel.ns = 20;
    matrix *mat = Obs_Dist(data, &jcmodel);
    for (i = 0; i < data->n_otu; i++)
    {
        for (j = i + 1; j < data->n_otu; j++)
        {
            x = a*(pow(1 - mat->dist[i][j], -1/a) - 1);
            mat->dist[i][j] = x;
            mat->dist[j][i] = mat->dist[i][j];
        }
    }
    return mat;
}

matrix *calc_dist_matrix(allseq *seqs, int distkind, int protein, double alpha,
    char **p_err_mess)
{
    matrix *phyml_mat = NULL;
    if (distkind == Obs_pdist)
    {//Observed
        model jcmodel;
        jcmodel.stepsize = 1;
        jcmodel.datatype = (protein ? 1 : 0);
        jcmodel.ns = (protein ? 20 : 4);
        phyml_mat = Obs_Dist(seqs, &jcmodel);
    }
    else if (distkind == Obs_gaps_pdist)
    {//Observed with gap as character
        model jcmodel;
        jcmodel.stepsize = 1;
        jcmodel.datatype = (protein ? 1 : 0);
        jcmodel.ns = (protein ? 20 : 4);
        phyml_mat = Obs_Dist_Gaps(seqs, &jcmodel);
    }
    else if (distkind == Kimura_pdist)
    {//Kimura
      phyml_mat = Kimura_p_Dist(seqs);
      if (phyml_mat == NULL)
          return NULL;
    }
    else if (distkind == Gamma_pdist)
    {//Gamma
      phyml_mat = Gamma_p_Dist(seqs, alpha);
      if (phyml_mat == NULL)
          return NULL;
    }
    else if (distkind == Poisson_pdist)
    {//Poisson
        model jcmodel;
        jcmodel.stepsize = 1;
        jcmodel.datatype = (protein ? 1 : 0);
        jcmodel.ns = (protein ? 20 : 4);
        phyml_mat = JC69_Dist(seqs, &jcmodel);
        if (phyml_mat == NULL)
            return NULL;
    }
    return phyml_mat;
}

void majuscules(char *name)
{
    name--;
    while(*(++name) != 0)
        *name = toupper(*name);
}

allseq *view_to_allseq(SEA_VIEW *view, int remove_all_gaps)
{
    int i, j, l;
    list_segments *ls;
    char *p;
    allseq *phyml_seqs = (allseq *)calloc(1, sizeof(allseq));
    phyml_seqs->n_otu = view->tot_sel_seqs == 0 ? view->tot_seqs : view->tot_sel_seqs;
    phyml_seqs->c_seq = (struct __Seq **)calloc(phyml_seqs->n_otu, sizeof(struct __Seq*));
    
    l = 0;
    if (view->active_region == NULL)
    {
        for (i = 0; i < view->tot_seqs; i++)
        {
            if (view->tot_sel_seqs != 0 && !view->sel_seqs[i])
                continue;
            if (view->each_length[i] > l)
                l = view->each_length[i];
        }
    }
    else
    {
        ls = view->active_region->list;
        while (ls != NULL)
        {
            l += ls->fin - ls->debut + 1;
            ls = ls->next;
        }
    }
    phyml_seqs->clean_len = l;
    phyml_seqs->wght = (int *)calloc(1, l * sizeof(int));
    for (i = 0; i < l; i++)
        phyml_seqs->wght[i] = 1;
    j = 0;
    for (i = 0; i < view->tot_seqs; i++)
    {
        if (view->tot_sel_seqs != 0 && !view->sel_seqs[i])
            continue;
        phyml_seqs->c_seq[j] = (struct __Seq *)calloc(1, sizeof(struct __Seq));
        phyml_seqs->c_seq[j]->name = view->seqname[i];
        phyml_seqs->c_seq[j]->len = l;
        phyml_seqs->c_seq[j]->state = (char *)malloc(l + 1);
        if (view->active_region == NULL)
        {
            memcpy(phyml_seqs->c_seq[j]->state, view->sequence[i], view->each_length[i] );
            if (l > view->each_length[i])
            {
                memset(phyml_seqs->c_seq[j]->state + view->each_length[i], '-', l -
                    view->each_length[i]);
            }
        }
        else
        {
            ls = view->active_region->list;
            p = phyml_seqs->c_seq[j]->state;
            while (ls != NULL)
            {
                if (ls->fin <= view->each_length[i])
                {
                    memcpy(p, view->sequence[i] + ls->debut - 1, ls->fin - ls->debut + 1 );
                }
                else
                {
                    int lrem = view->each_length[i] - ls->debut + 1;
                    if (lrem > 0)
                        memcpy(p, view->sequence[i] + ls->debut - 1, lrem );
                    if (lrem < 0)
                        lrem = 0;
                    memset(p + lrem, '-', ls->fin - ls->debut + 1 - lrem);
                }
                p += ls->fin - ls->debut + 1;
                ls = ls->next;
            }
        }
        phyml_seqs->c_seq[j]->state[l] = 0;
        majuscules(phyml_seqs->c_seq[j]->state);
        p = phyml_seqs->c_seq[j]->state;
        if (!view->protein)
        {
            while ((p = strchr(p, 'U')) != NULL)
                *p = 'T';
        }
        else
        {
            while ((p = strchr(p, '*')) != NULL)
                *p = '-';
        }
        if (!view->protein)
        {// remove non-nucleotide characters
            p = phyml_seqs->c_seq[j]->state;
            while (*p)
            {
                if (strchr("ABCDGHKMNRSTUVWXY-", *p) == NULL)
                    *p = 'N';
                p++;
            }
        }
        j++;
    }
    //remove gap-only or gap-with sites
    for (j = 0; j < phyml_seqs->clean_len; j++)
    {
        if (remove_all_gaps)
        {//remove any gap-containing site
            for (i = 0; i < phyml_seqs->n_otu; i++)
            {
                if (phyml_seqs->c_seq[i]->state[j] == '-')
                    break;
            }
            if (i == phyml_seqs->n_otu)
                continue;
        }
        else
        {//remove gap-only sites
            for (i = 0; i < phyml_seqs->n_otu; i++)
            {
                if (phyml_seqs->c_seq[i]->state[j] != '-')
                    break;
            }
            if (i != phyml_seqs->n_otu)
                continue;
        }
        for (i = 0; i < phyml_seqs->n_otu; i++)
        {
            memmove(phyml_seqs->c_seq[i]->state + j, phyml_seqs->c_seq[i]->state + j + 1,
                phyml_seqs->clean_len - j);
        }
        j--;
        phyml_seqs->clean_len--;
    }
    return phyml_seqs;
}

void check_symmetry(const MatrixXd& S)
{
    int i, j;
    double diff;
    char err_message[] = "Error: exchangeabilities matrix not symmetric";
    
    for (i = 0; i < 20; i++)
    {
        for (j = i + 1; j < 20; j++)
        {
            diff = S(i, j) - S(j, i);
            if (abs(diff) > 0.001)
                exit_error(err_message);
        }
    }
}

void set_exchanges(int model, MatrixXd& S)
{
    int i, j;
    
    if (model == LG_dist)
    {
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                S(i, j) = S_LG[i][j];
    }
    else if (model == JTT_dist)
    {
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                S(i, j) = S_JTT[i][j];
    }
    else if (model == WAG_dist)
    {
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                S(i, j) = S_WAG[i][j];
    }
    else if (model == PAM_dist)
    {
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                S(i, j) = S_PAM[i][j];
    }
    else if (model == Blosum_dist)
    {
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                S(i, j) = S_BLO[i][j];
    }
    //check_symmetry(S);
}

void set_equilibrium(int model, int freq, SEA_VIEW view, VectorXd& E)
{
    int i, ltot, n, m, test = FALSE;
    double sum = 0.0;
    char *ptr;
    char aa[] = "ARNDCQEGHILKMFPSTWYV";

    if (!freq)
    {
        if (model == PAM_dist)
        {
            for (i = 0; i < 20; i++)
                E(i) = E_PAM[i];
        }
        else if (model == JTT_dist)
        {
            for (i = 0; i < 20; i++)
                E(i) = E_JTT[i];
        }
        else if (model == WAG_dist)
        {
            for (i = 0; i < 20; i++)
                E(i) = E_WAG[i];
        }
        else if (model == LG_dist)
        {
            for (i = 0; i < 20; i++)
                E(i) = E_LG[i];
        }
        else if (model == Blosum_dist)
        {
            for (i = 0; i < 20; i++)
                E(i) = E_BLO[i];
        }
/*
        Rescale equilibrium frequencies to be sure that their sum = 1
*/
        for (i = 0; i < 20; i++)
            sum += E(i);
        for (i = 0; i < 20; i++)
            E(i) /= sum;
    }
/*
    If equilibrium frequencies are taken from dataset
*/    
    else
    {
        ltot = 0;
        for (n = 0; n < view.tot_seqs; n++)
        {
            for (m = 0; m < view.each_length[n]; m++)
            {
                ptr = strchr(aa, view.sequence[n][m]);
                if (ptr != NULL)
                {
                    i = (int)(ptr - aa);
                    E(i)++;
                    ltot++;
                }
            }
        }
        for (i = 0; i < 20; i++)
        {
            if (E(i) == 0)
                test = TRUE;
            E(i) /= double(ltot);
        }
        if (test)
            fprintf(stderr, "Warning: one or more equilibrium frequency is equal to zero\n");
    }
}

double calc_dist_model(double p, const VectorXd& K, const VectorXd& E, const MatrixXd& U,
    const MatrixXd& V)
{
    MatrixXd L = MatrixXd::Zero(20, 20);
    MatrixXd Pd(20, 20);
    int i, j, k;
    double a, b, d, phi, epsilon = 1e-6;
/*    
    if (p == 0)
    {
        return(0.0);
    }
*/
    for (k = int(p*100); k <= MAXDIST; k++)
    {
        for (i = 0; i < 20; i++)
            L(i, i) = pow(K(i), double(k));
        Pd = U*L*V;
        phi = 0.0;
        for (i = 0; i < 20; i++)
            phi += E(i)*Pd(i, i);
        phi = 1 - phi;
        if (phi > p)
            break;
    }
    a = double(k - 1);
    b = double(k);
    if (b >= MAXDIST)
        return(NAN);
    while (abs(phi - p) > epsilon)
    {
        d = (a+b)/2;
        for (i = 0; i < 20; i++)
            L(i, i) = pow(K(i), d);
        Pd = U*L*V;
        phi = 0.0;
        for (i = 0; i < 20; i++)
            phi += E(i)*Pd(i, i);
        phi = 1 - phi;
        if (abs(phi - p) <= epsilon)
            break;
        if (phi < p)
            a = d;
        else if (phi > p)
            b = d;
    }
    return(d/100);
}

void write_dst(FILE *fp, matrix *mat)
{
    int n, m;
    
    fprintf(fp, "%5d\n", mat->n_otu);
    for (n = 0; n < mat->n_otu; n++)
    {
        fprintf(fp, "%.10s  ", mat->name[n]);
        for (m = 0; m < mat->n_otu; m++)
            fprintf(fp, "%f ", mat->dist[n][m]);
        fprintf(fp, "\n");
    }
}

void test_duplicated_names(allseq *data)
{
    int i, j;
    char err_message[200];
    
    for (i = 0; i < data->n_otu - 1; i++)
    {
        for (j = i + 1; j < data->n_otu; j++)
        {
            if (strcmp(data->c_seq[i]->name, data->c_seq[j]->name) == 0)
            {
                sprintf(err_message, "Cannot run bootstrap because sequence name %s is used twice",
                    data->c_seq[i]->name);
                exit_error(err_message);
            }
        }
    }
}

void dist_miss_ultra(matrix *mat, double dmax)
{
    double d, Max, MinMax;
    int i, j, k;
    
    for (i = 0; i < mat->n_otu; i++)
    {
        for (j = i + 1; j < mat->n_otu; j++)
        {
            d = mat->dist[i][j];
            if (isnan(d) || isinf(d))
            {
                MinMax = dmax;
                for (k = 0; k < mat->n_otu; k++)
                {
                    if (k != i && k != j && !isnan(mat->dist[i][k]) &&
                    !isnan(mat->dist[j][k]) && !isinf(mat->dist[i][k]) &&
                    !isinf(mat->dist[j][k]))
                    {
                        Max = fmax(mat->dist[i][k], mat->dist[j][k]);
                        if (Max < MinMax)
                            MinMax = Max;
                    }
                }
                mat->dist[i][j] = mat->dist[j][i] = MinMax;
            }
        }
    }
}

void dist_miss_add(matrix *mat, double dmax)
{
    double d, Max, MinMax;
    int i, j, k, l;
    
    for (i = 0; i < mat->n_otu; i++)
    {
        for (j = i + 1; j < mat->n_otu; j++)
        {
            d = mat->dist[i][j];
            if (isnan(d) || isinf(d))
            {
                MinMax = dmax;
                for (k = 0; k < mat->n_otu; k++)
                {
                    for (l = k + 1; l < mat->n_otu; l++)
                    {
                        if (k != i && k != j && l != i && l != j &&
                        !isnan(mat->dist[i][k]) && !isnan(mat->dist[j][k]) &&
                        !isnan(mat->dist[i][l]) && !isnan(mat->dist[j][l]) &&
                        !isnan(mat->dist[k][l]) && !isinf(mat->dist[i][k]) &&
                        !isinf(mat->dist[j][k]) && !isinf(mat->dist[i][l]) &&
                        !isinf(mat->dist[j][l]) && !isinf(mat->dist[k][l]))
                        {
                            Max = fmax(mat->dist[i][k] + mat->dist[j][l], mat->dist[i][l] +
                                mat->dist[j][k] - mat->dist[k][l]);
                            if (Max < MinMax)
                                MinMax = Max;
                        }
                    }
                    
                }
                mat->dist[i][j] = mat->dist[j][i] = MinMax;
            }
        }
    }
}

char *get_name(char *str)
{
    char *ptr, *name;
    int i, c;
    
    ptr = strchr(str, '.');
    if (ptr == NULL)
        return(str);
    else
    {
        i = ptr - str;
        name = (char *)malloc(i);
        strncpy(name, str, i);
        return(name);
    }
}

void compute_dst(FILE *fp_dst, allseq *data, SEA_VIEW view, int estim, int boot,
    int model, double alpha, int freq, unsigned seed)
{
    MatrixXd L = MatrixXd::Zero(20, 20);
    MatrixXd F(20, 20), S(20, 20), Q(20, 20), U(20,20), V(20, 20), P(20, 20);
    VectorXd E(20), K(20);
    int i, n, m, repl, nnan;
    double d, dmax, s = 0, t = 0.01;
    char *err_message;
    char ret_message[200];
    matrix *mat;

    repl = 0;
/*
    Computation of distances
*/
    if (model == Kimura_pdist || model == Gamma_pdist || model == Poisson_pdist ||
        model == Obs_pdist || model == Obs_gaps_pdist)
    {
        do
        {
            mat = calc_dist_matrix(data, model, 1, alpha, &err_message);
            if (estim != 0)
            {
                nnan = 0;
                dmax = 0.0;
                for (n = 0; n < mat->n_otu; n++)
                {
                    for (m = n + 1; m < mat->n_otu; m++)
                    {
                        d = mat->dist[n][m];
                         if (isnan(d) || isinf(d))
                             nnan++;
                         else if (d > dmax)
                             dmax = d;
                    }
                }
                if (nnan != 0 && estim == 1 && mat->n_otu >= 4)
                    dist_miss_ultra(mat, dmax);
                else if (nnan != 0 && estim == 2 && mat->n_otu >= 5)
                    dist_miss_add(mat, dmax);
            }
            write_dst(fp_dst, mat);
            if (boot == 0)
                break;
            bootstrap_weights(data, seed);
            repl++;
        }
        while (repl <= boot);
    }
    else
    {
/*
        Construction of instantaneous rate matrix Q
*/
        set_exchanges(model, S);
        set_equilibrium(model, freq, view, E);
        F = E.asDiagonal();
        Q = S*F;
        for (i = 0; i < 20; i++)
            Q(i,i) = -Q.row(i).sum();
/*
        Renormalisation such as sum_{i}(pi_{i}*q_{ii}) = -1
*/
        for (i = 0; i < 20; i++)
            s -= E(i)*Q(i, i);
        Q = S*F/s;
        for (i = 0; i < 20; i++)
            Q(i,i) = -Q.row(i).sum();
/*
        Diagonalisation of Q and construction of P(t=0.01)
*/
        SelfAdjointEigenSolver < MatrixXd > es(Q);
        if (es.info() != Success)
        {
            strcpy(ret_message, "Error in instantaneous rate matrix diagonalisation");
            exit_error(ret_message);
        }
        U = es.eigenvectors();
        for (i = 0; i < 20; i++)
           L(i, i) = exp(es.eigenvalues()(i)*t);
        V = U.inverse();
        P = U*L*V;
/*
        Diagonalisation of P(t=0.01)
*/
        SelfAdjointEigenSolver < MatrixXd > es2(P);
        if (es2.info() != Success)
        {
            strcpy(ret_message, "Error in transition probabilities matrix diagonalisation");
            exit_error(ret_message);
        }
        for (i = 0; i < 20; i++)
            K(i) = es2.eigenvalues()(i);
        U = es2.eigenvectors();
        V = U.inverse();
/*
        For all p-distances in mat computes d
*/
        do
        {
            nnan = 0;
            dmax = 0.0;
            mat = calc_dist_matrix(data, Obs_pdist, 1, alpha, &err_message);
            #ifdef _OPENMP
            #pragma omp parallel private(n, d)
            #endif
            for (n = 0; n < mat->n_otu; n++)
            {
                #ifdef _OPENMP
                #pragma omp for
                #endif
                for (m = n + 1; m < mat->n_otu; m++)
                {
                     d = calc_dist_model(mat->dist[n][m], K, E, U, V);
                     if (isnan(d) || isinf(d))
                         nnan++;
                     else if (d > dmax)
                         dmax = d;
                     mat->dist[n][m] = d;
                }
            }
            for (n = 0; n < mat->n_otu - 1; n++)
                for (m = n + 1; m < mat->n_otu; m++)
                    mat->dist[m][n] = mat->dist[n][m];
            if (nnan != 0 && estim == 1 && mat->n_otu >= 4)
                dist_miss_ultra(mat, dmax);
            else if (nnan != 0 && estim == 2 && mat->n_otu >= 5)
                dist_miss_add(mat, dmax);
            write_dst(fp_dst, mat);
            if (boot == 0)
                break;
            bootstrap_weights(data, seed);
            repl++;
        }
        while (repl <= boot);
    }
    Free_Mat(mat);
}

void compute_trees(FILE *fp_dst, FILE *fp_tree)
{
    POINTERS *trees;   /* list of subtrees            */
    char *chain1;      /* stringized branch-length    */
    char *chain2;      /* idem                        */
    int *a, *b;        /* pair to be agglomerated     */
    float **delta;     /* delta matrix                */
    float la;          /* first taxon’s branch-length */
    float lb;          /* second taxon’s branch-length*/
    float vab;         /* variance of Dab             */
    float lamda;
    int i;
    int ok;
    int r;             /* number of subtrees          */
    int n;             /* number of taxa              */
    int x, y;
    double t;

    a = (int *)my_calloc(1, sizeof(int));
    b = (int *)my_calloc(1, sizeof(int));
    chain1 = (char *)my_calloc(LEN, sizeof(char));
    chain2 = (char *)my_calloc(LEN, sizeof(char));
    fscanf(fp_dst, "%d", &n);  
/*
    Create the delta matrix
*/
    delta = (float **)my_calloc(n + 1, sizeof(float*));
    for (i = 1; i <= n; i++)
        delta[i] = (float *)my_calloc(n + 1, sizeof(float));
    trees = (POINTERS *)my_calloc(n + 1, sizeof(POINTERS));
/*
    Initialise and symmetrize the running delta matrix
*/ 
    rewind(fp_dst);
    while (fscanf(fp_dst,"%d", &n) != EOF )
    {
        r = n;
        *a = 0;
        *b = 0;
        Initialize(delta, fp_dst, n, trees);
        ok = Symmetrize(delta, n);
        if (!ok)
            fprintf(stderr, "Warning: a distance matrix is not symmetric\n");
        while (r > 3)
        {
            Compute_sums_Sx(delta, n);               /* compute the sum Sx       */
            Best_pair(delta, r, a, b, n);            /* find the best pair by    */
            vab = Variance(*a, *b, delta);           /* minimizing (1)           */
            la = Branch_length(*a, *b, delta, r);    /* compute branch-lengths   */
            lb = Branch_length(*b, *a, delta, r);    /* using formula (2)        */
            lamda = Lamda(*a, *b, vab, delta, n, r); /* compute lambda* using (9)*/
            for (i = 1; i <= n; i++)
            {
                if (!Emptied(i, delta) && (i != *a) && (i != *b))
                {
                    if (*a > i)
                    {
                        x = *a;
                        y = i;
                    }
                    else
                    {
                        x = i;
                        y = *a;             /* apply reduction formulae */
                    }                       /* 4 and 10 to delta        */
                    delta[x][y] = Reduction4(*a, la, *b, lb, i, lamda, delta);
                    delta[y][x] = Reduction10(*a, *b, i, lamda, vab, delta);
                }
            }
            strcpy(chain1, "");                    /* agglomerate the subtrees */
            strcat(chain1, "(");                   /* a and b together with the*/
            Concatenate(chain1, *a, trees, 0);     /* branch-lengths according */
            strcpy(chain1, "");                    /* to the NEWSWICK format */
            strcat(chain1, ":");
            sprintf(chain1 + strlen(chain1), "%f", la);
            strcat(chain1, ",");
            Concatenate(chain1,*a, trees, 1);
            trees[*a].tail->suiv = trees[*b].head;
            trees[*a].tail = trees[*b].tail;
            strcpy(chain1, "");
            strcat(chain1, ":");
            sprintf(chain1 + strlen(chain1), "%f", lb);
            strcat(chain1, ")");
            Concatenate(chain1, *a, trees, 1);
            delta[*b][0] = 1.0;                     /* make the b line empty     */
            trees[*b].head = NULL;
            trees[*b].tail = NULL;
            r = r - 1;                                /* decrease r                */
        }
        Finish(delta, n, trees, fp_tree);    /* compute the branch-lengths */
        for (i = 1; i <= n; i++)             /* of the last three subtrees */
        {                                    /* and print the tree in the */
            delta[i][0] = 0.0;               /* output-file */
            trees[i].head = NULL;
            trees[i].tail = NULL;
        }
    }
    free(delta);
    free(trees);
}

int main(int argc, char **argv)
{
    SEA_VIEW view;
    FILE *fp_seq, *fp_dst, *fp_tree, *fp_root;
    char **unused2, *unused, *err_message;
    char ret_message[200], *file_name, *dst_name, *tree_name, *root_name, *name_opt;
    bool name_opt_bool;
    allseq *data;
    double a, alpha = 2.25;
    int i, j, model, gaps, freq, boot, estim, trds, root, split;
    unsigned seed;
    size_t lgth, lgth_2;

    model = PAM_dist;
    seed = estim = boot = 0;
    trds = 1;
    gaps = freq = root = split = FALSE;
/*
    Arguments parsing
*/

    if (argc < 2)
    {
        sprintf(ret_message, "Usage: %s <infile> [-n <name>] [-m <model>] [-a <alpha>] [-e <0|1|2>] [-b <number>] [-s <number>] [-t <number>] [-f] [-g] [-x]",
            argv[0]);
        exit_error(ret_message);
    }
    if ((fp_seq = fopen(argv[1], "r")) == NULL)
    {
        sprintf(ret_message, "Error: unable to open sequence file %s", argv[1]);
        exit_error(ret_message);
    }
    file_name = get_name(argv[1]); /* Sets the generic file name for output */
    lgth = strlen(file_name);
    
    for (i = 1; i < argc; i++)
    {
        if (strstr(argv[i], "-m") != NULL)
        {
            if (strstr(argv[i+1], "PAM") != NULL)
                model = PAM_dist;
            else if (strstr(argv[i+1], "JTT") != NULL)
                model = JTT_dist;
            else if (strstr(argv[i+1], "WAG") != NULL)
                model = WAG_dist;
            else if (strstr(argv[i+1], "LG") != NULL)
                model = LG_dist;
            else if (strstr(argv[i+1], "Blo") != NULL)
                model = Blosum_dist;
            else if (strstr(argv[i+1], "Kim") != NULL)
                model = Kimura_pdist;
            else if (strstr(argv[i+1], "Gam") != NULL)
                model = Gamma_pdist;
            else if (strstr(argv[i+1], "Poi") != NULL)
                model = Poisson_pdist;
            else if (strstr(argv[i+1], "Obs") != NULL)
                model = Obs_pdist;
            else if (strstr(argv[i+1], "Gap") != NULL)
                model = Obs_gaps_pdist;
            else
            {
                strcpy(ret_message, "Error: available substitution models: PAM, JTT, WAG, LG, Blo, Kim, Gam, Poi, Obs or Gap");
                exit_error(ret_message);
            }
        }
        if (strncmp(argv[i], "-n", 2) == 0)
        {
            free(file_name);
	    file_name = argv[i+1]; /* Sets the generic file name for output */
	    lgth = strlen(file_name);
        }
        if (strstr(argv[i], "-a") != NULL)
        {
            if (i < argc - 1)
                if (sscanf(argv[i+1], "%lf", &a) == 1)
                    alpha = a;
            if (alpha <= 0 || alpha > 200)
            {
                strcpy(ret_message, "Error: value for alpha parameter must be >0 and <=200");
                exit_error(ret_message);
            }
        }
        if (strstr(argv[i], "-b") != NULL)
        {
            boot = 100;
            if (i < argc - 1)
                if (sscanf(argv[i+1], "%i", &j) == 1)
                    boot = j;
        }
        if (strstr(argv[i], "-s") != NULL)
        {
            if (i < argc - 1)
                if (sscanf(argv[i+1], "%i", &j) == 1)
                    seed = abs(j);
        }
        else if (strstr(argv[i], "-t") != NULL)
        {
            if (i < argc - 1)
                if (sscanf(argv[i+1], "%i", &j) == 1)
                    trds = j;
        }
        if (strstr(argv[i], "-e") != NULL)
        {
            if (i < argc - 1)
                if (sscanf(argv[i+1], "%i", &j) == 1)
                    estim = j;
            if (estim < 0 || estim > 2)
            {
                strcpy(ret_message, "Value for missing distance estimation method must be 0, 1 or 2");
                exit_error(ret_message);
            }
        }
        if (strstr(argv[i], "-f") != NULL)
            freq = TRUE;
        if (strstr(argv[i], "-g") != NULL)
            gaps = TRUE;
        if (strstr(argv[i], "-x") != NULL)
            split = TRUE;
    }
    #ifdef _OPENMP
    omp_set_num_threads(trds);
    #endif
/*
    Fasta file opening and reading
*/
    memset(&view, 0, sizeof(view));
    view.tot_seqs = read_fasta_align(fp_seq, &view.sequence, &view.seqname, &unused2,
        &unused, &err_message, 0);
    if (view.tot_seqs == 0)
        exit_error(err_message);
    view.protein = 1;
    view.each_length = (int *)malloc(view.tot_seqs * sizeof(int));
    for (i = 0; i < view.tot_seqs; i++)
        view.each_length[i] = strlen(view.sequence[i]);
    data = view_to_allseq(&view, gaps);
    if (boot != 0)
        test_duplicated_names(data);
/*
    Distances computing
*/
    dst_name = (char *)malloc(lgth + 6);
    sprintf(dst_name, "%s.dst%c", file_name, '\0');
    if ((fp_dst = fopen(dst_name, "w+")) == NULL)
    {
        sprintf(ret_message, "Error: unable to open distances output file %s",
            dst_name);
        exit_error(ret_message);
    }
    compute_dst(fp_dst, data, view, estim, boot, model, alpha, freq, seed);
/*
    BioNJ trees computing
*/
    tree_name = (char *)malloc(lgth + 5);
    sprintf(tree_name, "%s.phb%c",file_name, '\0');
    if ((fp_tree = fopen(tree_name, "w+")) == NULL)
    {
        sprintf(ret_message, "Error: unable to open trees output file %s",
            tree_name);
        exit_error(ret_message);
    }
    rewind(fp_dst);
    compute_trees(fp_dst, fp_tree);
    fclose(fp_dst);
    rewind(fp_tree);
    if (split)
	{
	    split_trees(fp_tree, file_name, lgth);
	}
        
/*
    BioNJ trees rooting
    root_name = (char *)malloc(lgth + 10);
    sprintf(root_name, "%s_root.phb\0", file_name);
    if ((fp_root = fopen(root_name, "w")) == NULL)
    {
        sprintf(ret_message, "Error: unable to open rooted trees output file %s",
                root_name);
        exit_error(ret_message);
    }
    rewind(fp_tree);
    root_trees(fp_tree, fp_root);
    fclose(fp_root);
*/
    fclose(fp_tree);
}

void split_trees(FILE *fp_tree, char *file_name, int lgth)
{
    FILE *fp_out;
    char c, ret_message[200], *trees, *tree, *out_name;
    long l;
    int i;
    fseek(fp_tree, 0, SEEK_END);
    l = ftell(fp_tree);
    fseek(fp_tree, 0, SEEK_SET);
    trees = (char *)my_calloc(l + 1, sizeof(char));
    char *p = trees;
    while ((c = fgetc(fp_tree)) != EOF)
    {
        if (c != '\n' && c != '\r')
            *(p++) = c;
    }
    *p = 0;
    tree = strtok(trees, ";");
    tree = strtok(NULL, ";"); // the first tree is not a bootstrapped one
    i = 0;
    while (tree != NULL)
    {
        out_name = (char *)malloc(lgth + 10);
        sprintf(out_name, "%s_%d.phb%c",file_name, i, '\0');
        if ((fp_out = fopen(out_name, "w")) == NULL)
        {
            sprintf(ret_message, "Error: unable to open tree output file %s", out_name);
            exit_error(ret_message);
        }
        fprintf(fp_out, "%s;\n", tree);
        i++;
        tree = strtok(NULL, ";");
        fclose(fp_out);
        free(out_name);
    }
    free(tree);
    free(trees);
 }

void root_trees(FILE *fp_tree, FILE *fp_root)
{
    char c, ret_message[200], *trees, *tree;
    long l;
    int i;
    
    fseek(fp_tree, 0, SEEK_END);
    l = ftell(fp_tree);
    fseek(fp_tree, 0, SEEK_SET);
    trees = (char *)my_calloc(l + 1, sizeof(char));
    char *p = trees;
    while ((c = fgetc(fp_tree)) != EOF)
    {
        if (c != '\n' && c != '\r')
            *(p++) = c;
    }
    *p = 0;
    tree = strtok(trees, ";");
    i = 0;
    while (tree != NULL)
    {
        fprintf(stderr, "tree[%d]:\n%s;\n", i, tree);
        i++;
        tree = strtok(NULL, ";");
    }
 }

int calc_tree_count(char *tree)
{
    int count = 0;
    
    while (TRUE)
    {
        while (*tree && isspace(*tree))
            tree++;
        if (*tree == '[')
            do tree++; while (*tree && *tree != ']');
        while (*tree && *tree != '(')
            tree++;
        if (*tree == '(')
            tree = nextpar(tree);
        while (*tree && *tree != ';')
            tree++;
        if (*tree == 0)
            break;
        count++;
    }
    return count;
}

char *nextpar(char *pospar)
{
    char *pos;
    
    pos = pospar + 1;
    while (*pos != ')')
    {
        if (*pos == 0)
            return NULL;
        if (*pos == '(')
            pos = nextpar(pos);
        pos++;
    }
    return pos;
}


