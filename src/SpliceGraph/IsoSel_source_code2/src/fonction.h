#ifndef FONCTION_H_INCLUDED
#define FONCTION_H_INCLUDED

#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <list>  
#include <map>
#include <algorithm>
#ifdef _OPENMP
    #include <omp.h>
#endif

// Fonctions de manipulation de fichiers
void rename_for_fastdist(std::ifstream& base_MSA_file,const std::string& TmpDirectory,bool DS);
std::string StrReplace(std::string& str, const std::string& oldStr, const std::string& newStr);
// void replace_NAN(std::ifstream& matrices_file, const std::string& WorkDirectory,int nb_seq,int boot,bool DS);
int Isoforms_couting(std::string& transcript_File_Name,const std::string& TmpDir,int boot,const std::string& output_name);
void checkTranscriptFile(std::string& transcript_File_Name,std::string& InfilePath, const std::string& TmpDir,int boot,const std::string& output_name);

// Fonctions de création de tableau
std::vector<short> create_species_line(std::string& species_str,std::string& InfilePath, const std::string& TmpDirectory,int boot,const std::string& output_name);
std::vector<std::string> create_align_tab(std::ifstream& file);
void create_align_tab_pert(std::ifstream& file, std::string* tab);
void create_tab_num_base(float** tab_num, std::vector<std::string> tab, unsigned nb_seq, unsigned nb_site);
bool stats_alignment(float** tab_num, unsigned nb_seq, unsigned nb_site);
float MeanComputation(std::vector<int> tailles_seq, unsigned nb_seq);
void create_tab_num_pert(short** tab_num, std::string* tab, unsigned nb_seq, unsigned nb_site);
void create_pair_tab(short *** pair_tab, float** tab_num, unsigned nb_seq, unsigned nb_site, unsigned nb_paire);
void create_pair_tab_pert(short *** pair_tab, short** tab_num, unsigned nb_seq, unsigned nb_site, unsigned nb_paire);
void create_residu_score_for1msa(float** tab_num_base,short** tab_num_pert, unsigned nb_seq, unsigned nb_site, float** tab_residu_score_for1msa, std::vector<short> seq_lines, float* tab_nb_non_gap,bool GapPen);
void create_tab_nb_non_gap(float* tab_nb_non_gap, float** tab_num_base, unsigned nb_seq, unsigned nb_site_base);
void create_tab_score_seq(float* tab_score_par_seq, float** tab_score_res, unsigned nb_species, unsigned nb_seq, unsigned nb_site, std::vector<short> species_lines, float* tab_nb_non_gap,bool ShortPen);

// Maps initialisations
std::map<std::string, int > initMapCorNameNum(const std::string& TmpDir);
std::map<int, std::string > initMapCorNumName(const std::string& TmpDir);
std::map<std::string, std::vector<int> > initMapGeneTrans(const std::string& TmpDir,std::string& transcript_File_Name,std::map<std::string,int>& MapCorNameNum);
std::map<int, std::vector<int> > initMapNumAltTrans(std::map<std::string, std::vector<int> >& MapGeneTrans);
std::vector<int> ListTranscriptNumber(std::map<int, std::vector<int> >& MapNumAltTrans);
void initTransNumSize(const std::string& TmpDir,std::map<int,float>& map_name_num);
// Function option distances
void filtered_file_distances(const std::string& output_name,std::string& ScoreFilePath, std::string& infilePath,const std::string& WorkDirectory,const std::string& TmpDir,std::map< std::string, int>& CorrNameNum,std::map<int, std::string >& CorrNumName,std::map<std::string, std::vector<int> >& MapGeneTrans,std::vector<int>& TranscriptsList);
void supp_fichier_distances(const std::string& WorkDirectory,const std::string& output_name);
void PrintDistFile(std::map<int, std::string >& MapCorNumName,std::ofstream& DistFile,float& DistSum,int SeqNum,int SeqSize);
float MedianComputation(std::vector<int>& TransNumList,std::map<int, float>& TransNumSize, int nb_seq,bool only_uniq_prot);


// Fonction de création de fichier de resultats
void create_result_file(float* tab_score_per_res, unsigned nb_seq, std::string argv, double treshold, const std::string& WorkDirectory);
void create_seq_score_file_with_lines(float* tab_score_per_res, std::string& InfilePath, std::vector<short> species_lines, unsigned nb_species,std::string out_name, const std::string& WorkDirectory,const std::string& TmpDirectory);

void supp_fichier(const std::string& WorkDirectory,int boot,const std::string& output_name);
void check_nb_thread(unsigned int nb_thread);
void check_nb_boot(unsigned int nb_boot);
void MsaOrder(const std::string& TmpDirectory,int nb_tree);
void BaseOrder(const std::string& TmpDirectory,const std::string& WorkDirectory,const std::string& output_name);

int Seq_name_to_Number(std::string& infilePath, const std::string& TmpDirectory, const std::string& WorkDirectory,int boot,const std::string& output_name,bool DS);
void Seq_num_to_seq_names(const std::string& TmpDir,const std::string& WorkDirectory,const std::string& output_name,int boot);
void MafftTreeFormating(const std::string& TmpDir,int nb_tree);
std::string checkTree(std::string& tree);
std::string FindSisterLeafs(std::string& Tree,std::ofstream& RootedForMafft);
void filtered_file(const std::string& output_name,std::string& TranscriptFilePathName,std::string& ScoreFilePath,std::string& infilePath,const std::string& WorkDirectory,const std::string& TmpDir,int boot);
void find(std::string str);
#endif // FONCTION_H_INCLUDED
